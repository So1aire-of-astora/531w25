library(doParallel)
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores)) cores <- detectCores()  
registerDoParallel(cores)

system.time(
 rnorm(10^8)
) -> time0

system.time(
  foreach(i=1:10) %dopar% rnorm(10^7)
) -> time1

system.time(
  foreach(i=1:10^2) %dopar% rnorm(10^6)
) -> time2

system.time(
  foreach(i=1:10^3) %dopar% rnorm(10^5)
) -> time3

system.time(
  foreach(i=1:10^4) %dopar% rnorm(10^4)
) -> time4

results <- rbind(
  time0,
  time1,
  time2,
  time3,
  time4
)

rownames(results) <- paste0("time", 0:4)

results_df <- as.data.frame(results)
results_df$cores <- cores

write.csv(results_df, file = "test.csv")



