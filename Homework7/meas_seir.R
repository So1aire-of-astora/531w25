library(tidyverse)
library(pomp)
library(foreach)
# library(doFuture)
library(parallel)
library(doParallel)
# library(iterators)

# plan(multisession)
cl <- makeCluster(detectCores() - 1) 
registerDoParallel(cl)

set.seed(1350254336)

# The code for the SEIR model is developed from https://kingaa.github.io/sbied/pfilter/model.R

read_csv(paste0("https://kingaa.github.io/sbied/stochsim/",
  "Measles_Consett_1948.csv")) |>
  select(week,reports=cases) -> meas


seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E, 1 - exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

seir_rinit <- Csnippet("
  S = nearbyint(eta*N);
  E = 0; 
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

seir_dmeas <- Csnippet("
  lik = dnbinom_mu(reports,k,rho*H,give_log);
")

seir_rmeas <- Csnippet("
  reports = rnbinom_mu(k,rho*H);
")


meas |> pomp(
    times="week",
    t0=0,
    rprocess=euler(seir_step,delta.t=1/7),
    rinit=seir_rinit,
    rmeasure=seir_rmeas,
    dmeasure=seir_dmeas,
    accumvars="H",
    statenames=c("S","E", "I","R","H"),
    paramnames=c("Beta", "mu_EI", "mu_IR","N","eta","rho","k")
  ) -> measSEIR

measSEIR |>
  simulate(
    params=c(Beta=7.5, mu_EI = .8, mu_IR=1.5, rho=0.5, k=10,
      eta=0.05, N=38000),
    nsim=20, format="data.frame", include.data=TRUE
  ) -> sims

sims |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color="none")


# measSEIR |>
#   pfilter(Np=1000) -> pf


fixed_params <- c(N=38000, mu_EI = .8, mu_IR=2, k=10)
coef(measSEIR,names(fixed_params)) <- fixed_params




tic <- Sys.time()
foreach(i=1:10,.combine=c,
#   .options.future=list(seed=TRUE)
    .packages = "pomp"
) %dopar% {
  measSEIR |> pfilter(Np=5000)
} -> pf

pf |> logLik() |> logmeanexp(se=TRUE) -> L_pf
L_pf
toc <- Sys.time()

pf[[1]] |> coef() |> bind_rows() |>
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2]) |>
  write_csv("measles_params.csv")



## LOCAL SEARCH

## What is this 'bake' function?
## See https://kingaa.github.io/sbied/pfilter/bake.html
## for an explanation.

coef(measSEIR) <- c(
  Beta = 7.5, rho = 0.5, eta = 0.05,
  N = 38000, mu_EI = 0.8, mu_IR = 2, k = 10
)

bake(file="local_search.rds",{
  foreach(i=1:20,.combine=c,
    .packages = "pomp"
  ) %dopar% {
    measSEIR |>
      mif2(
        Np=2000, Nmif=50,
        cooling.fraction.50=0.5,
        rw.sd=rw_sd(Beta=0.02, rho=0.02, eta=ivp(0.02)),
        partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
        paramnames=c("Beta","rho","eta")
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- getDoParWorkers()
  mifs_local
}) -> mifs_local

t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")

mifs_local |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")



bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.combine=rbind,
    .options.future=list(seed=900242057)
  ) %dopar% {
    evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf |> coef() |> bind_rows() |>
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- nbrOfWorkers()
  results
}) -> results
t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")

pairs(~loglik+Beta+eta+rho,data=results,pch=16)

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

# if (file.exists("CLUSTER.R")) {
#   source("CLUSTER.R")
# }


## GLOBAL SEARCH

runif_design(
  lower=c(Beta=5,rho=0.2,eta=0),
  upper=c(Beta=80,rho=0.9,eta=1),
  nseq=400
) -> guesses

mf1 <- mifs_local[[1]]


guesses = head(guesses) # test line

bake(file="global_search.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=1270401374)
    ) %dopar% {
      mf1 |>
        mif2(params=c(guess,fixed_params)) |>
        mif2(Nmif = 20) -> mf # Nmif default: 100
      replicate(
        2, #default: 10
        mf |> pfilter(Np = 200) |> logLik() # Np default: 5000
      ) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) |>
  filter(is.finite(loglik)) -> results

t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-50) |>
  bind_rows(guesses) |>
  mutate(type=if_else(is.na(loglik),"guess","result")) |>
  arrange(type) -> all

pairs(~loglik+Beta+eta+rho, data=all, pch=16, cex=0.3,
  col=ifelse(all$type=="guess",grey(0.5),"red"))

all |>
  filter(type=="result") |>
  filter(loglik>max(loglik)-10) |>
  ggplot(aes(x=eta,y=loglik))+
  geom_point()+
  labs(
    x=expression(eta),
    title="poor man's profile likelihood"
  )

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-20,loglik.se<2) |>
  sapply(range) -> box

box

freeze(seed=1196696958,
  profile_design(
    eta=seq(0.01,0.95,length=40),
    lower=box[1,c("Beta","rho")],
    upper=box[2,c("Beta","rho")],
    nprof=15, type="runif"
  )) -> guesses
plot(guesses)


