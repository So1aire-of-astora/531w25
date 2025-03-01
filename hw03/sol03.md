---
title: "Solution to Homework 3"
author: "STATS/DATASCI 531, Winter 2025"
output:
  html_document:
    toc: no
bibliography: sol03.bib
csl: https://raw.githubusercontent.com/citation-style-language/styles/master/cell-numeric.csl
---



This analysis is developed from a previous homework submission [@rehnberg]. Homework reports are welcome to learn from previous solutions, as long as they are properly cited. However, your report is expected to demonstrate independent work that contributes beyond any sources. There are many different choices to make when carrying out data analysis, and this solution presents just one one set of choices. Also, this report was developed based on data only up to 2017. 

Most people noticed that it is hard to get a model fitting better than white noise for these data. Trying to do better can lead to unstable models with weakly identified parameters.






The data consist of the low temperature from each January, ranging from 1900 to 2024, in Ann Arbor, MI. By looking only at temperatures from the same month each year, we have simplified our problem by eliminating seasonal fluctuation. However, this also reduces our available data. Additionally, the low temperature is not available for January 1955 (seen in the plot below). There are various possibilities to deal with this missing value. Here, we follow [@rehnberg] by using the 1954 January temperature as a proxy for the missing 1955 value. 

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

This time series plot of the observed data shows wide variation in the January low temperatures across the years, ranging from below $-20^\circ F$ to almost $20^\circ F$. The data appear to fluctuate around the mean of -2.800$^\circ F$ without any obvious long-term patterns. Based on this, it seems reasonable to begin with a null hypothesis that a model with no trend is appropriate for the data. This analysis won't look at any models with trend, but that would be a logical next step.

From the time series plot, it also seems possible that the variation in the data is increasing over time, especially from about 1980 to present. Additional work could investigate whether a trend in variance is statistically supported.

We start the analysis by using maximum likelihood to fit an ARMA(p,q) model of the form:

$$ Y_n = \mu + \phi_1(Y_{n-1} - \mu) + \dots + \phi_p(Y_{n-p} - \mu) + \varepsilon_n + \psi_1\varepsilon_{n-1} + \dots + \psi_q\varepsilon_{n-q}$$

where ${\{\varepsilon_n}\}$ is a white noise process with distribution $\mathcal{N}(0,\sigma^2)$. The parameters for this model are $\theta = (\phi_1, \dots, \phi_p, \psi_1, \dots, \psi_q, \mu, \sigma^2)$, representing the coefficients for the autoregressive part of the model, the coefficients for the moving average part of the model, the population mean, and the error variance. In this model, $\mu$ does not depend on time because we are assuming a model without trend. To determine the best values of p and q for this data, we fit multiple ARMA(p,q) models with various values of p and q (shown below ranging from 0 to 4).

As an initial method to compare these various ARMA models, we will consider their Akaike information criteria (AIC) values following the approach in Chapter 5 of the notes [@notes531].
Models with low values of the AIC indicate higher prediction precision, and therefore, better models in terms of predictive power.
Though this is a somewhat informal method of model selection, it can be effective at eliminating models with very bad fits.


|            |    MA0|    MA1|    MA2|    MA3|    MA4|
|:-----------|------:|------:|------:|------:|------:|
|<b> AR0</b> | 861.02| 862.92| 864.75| 866.06| 867.61|
|<b> AR1</b> | 862.93| 864.61| 866.46| 867.09| 869.02|
|<b> AR2</b> | 864.78| 866.45| 868.35| 868.91| 871.09|
|<b> AR3</b> | 866.12| 867.39| 864.45| 870.82| 872.66|
|<b> AR4</b> | 867.91| 869.30| 871.38| 866.69| 868.76|

In the AIC table above, the lowest value is associated with an ARMA(0,0) model. This is a white noise model that assumes no temperature dependence between subsequent years. The ARMA(0,0) model is of the form $Y_n = \mu + \varepsilon_n$, where the ${\{\varepsilon_n}\}$ are as described above. Although the AIC table identifies this model as the most appropriate for the data, there is some climatological intuition indicating that there is dependence in temperature from one year to the next. Therefore, we will not restrict my analysis to the ARMA(0,0) model. In addition, we will look at some models that have a higher AIC, but that allow for the dependence we are interested in modeling, including the ARMA(0,1), ARMA(1,0), and ARMA(1,1).

We fit these four models to the data and the results are listed in the table below. The first thing to notice is that all four models give similar estimates for the intercept, around -2.8 but that their standard error estimates increase somewhat with model complexity, varying from $0.68$ for the ARMA(0,0) to $0.81$ for the ARMA(1,1).


|               |Intercept |SE(Intercept) |AR Coef. |MA Coef. |
|:--------------|:---------|:-------------|:--------|:--------|
|<b> ARMA(0, 0) |-2.800    |0.667         |--       |--       |
|<b> ARMA(0, 1) |-2.802    |0.687         |--       |0.030    |
|<b> ARMA(1, 0) |-2.802    |0.686         |0.028    |--       |
|<b> ARMA(1, 1) |-2.828    |0.787         |0.825    |-0.793   |

This may indicate that the ARMA(1,1) is more accurately capturing the dependence in the data than the other three models. Inadequately modeling dependence can result in artificially low standard errors for parameter estimates. These results indicate that the ARMA(1,1) is the better model to use, which is in opposition to the results of the AIC table above.

Due to the results of the AIC table and these fitted values, I will continue to consider the ARMA(0,0) model and the ARMA(1,1). The other two models, ARMA(0,1) and ARMA(1,0) have coefficient estimates very close to zero and don't seem to be doing anything significantly different from the ARMA(0,0). The ARMA(0,0) model can be written as $Y_n =  -2.800 + \varepsilon_n$, and the ARMA(1,1) model can be written as follows: 

$$\phi(B)(Y_n - (-2.828)) = \psi(B)\varepsilon_n$$

where $B$ is the backshift operator, and $\phi(x)$ and $\psi(x)$ are the AR and MA polynomials, respectively. For this fitted model, these polynomials are defined as follows:
$$\phi(x) = 1 - 0.825 x \hspace{3cm} \psi(x) = 1 -0.793 x$$



Something to consider with the ARMA(1,1) model are the roots of the AR and MA polynomials, which can be used to check for causality and invertibility.
The AR root is $1.212$ and the MA root is $1.262$. Both outside the unit circle, indicating that the fitted model is both causal and invertible, two attractive qualities for a time series model.
However, these roots are also relatively close in magnitude, which indicates the possibility of reducing the model to the ARMA(0,0).
It is difficult to tell if these roots are close enough to approximately cancel, but it definitely seems like a possibility.
This is another argument for using the ARMA(0,0) model over the ARMA(1,1).

A final test that we will do is a more formal hypothesis test using Wilks' approximation.
For this test, the null hypothesis corresponds to the ARMA(0,0) model, while the alternative corresponds to the ARMA(1,1). The approximation tells us:

$$\Lambda = 2(\mathcal{l}_1 - \mathcal{l}_0) \approx \chi^2_{D_1-D_0}$$



where $\mathcal{l}_i$ is the maximum log likelihood under hypothesis $H_i$ and $D_i$ is the number of parameters estimated under hypothesis $H_i$.
We will reject the null hypothesis if $\Lambda$ is larger than the $\chi^2$ cutoff.
When comparing ARMA(0,0) and ARMA(1,1), $\Lambda = 0.41$, which we can compare to the cutoff value of $5.99$ for a 95% significance level and 2 degrees of freedom.
Thus, this test does not provide evidence against our null hypothesis that the ARMA(0,0) model is adequate for the data.
Since this conclusion is supported both here, with the Wilks' approximate $\chi^2$ test, with the approximately canceling roots, and with the AIC, we will move forward with the white noise model.


Since we have identified the ARMA(0,0) as the best candidate model for the data, we should check that the model assumptions are valid.
First, we will look at the residuals of the fitted ARMA(0,0) model as a time series plot:

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

Any pattern in the residuals would be evidence against the model, but the time series plot shows no striking patterns.
Next, we can look at the autocorrelation plot of the residuals.
This will allow us to check our assumption that the errors $\{\varepsilon_n\}$ are uncorrelated.
There is only one lag with significant autocorrelation (lag 15), while the rest may be considered sufficiently close to zero.
While this may be the case, there are also some potentially non-negligible fluctuations in the autocorrelation that might be interesting to look into more carefully.
Perhaps this indicates that a model with trend could be appropriate for this data.

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

Finally, in fitting an ARMA model, we make the assumption that $\{\varepsilon_n\} \sim \mathcal{N}(0,\sigma^2)$ and we can check the normality assumption with a QQ-plot of the residuals. With the exception of a few points that deviate from the line, the residuals seem to be sufficiently normal to make this assumption valid.

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

Since the model fit seems to meet the assumptions, we can consider doing inference on the parameter estimate for $\mu$. The $\texttt{arima()}$ function in R uses the observed Fisher information to calculate standard errors for the coefficients. Those standard errors can then be used to construct approximate 95% confidence intervals for the parameter estimates. 



$$[-2.800 - (1.96)(0.667),
   -2.800 + (1.96)(0.667)]
   = [ -4.107 , -1.493 ] $$

The confidence interval for the mean does not contain zero, but this does not have much scientific meaning as a null hypothesis. Why?

As noted above, there is a possibility that the standard error from the ARMA(0,0) model ($0.667$) was artificially small. Therefore, I can check this confidence interval through simulations. Here is the distribution of the estimate for $\mu$ from 5,00 simulations:




![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

In this plot, the dashed vertical lines correspond to the upper and lower limits of the Fisher information confidence interval calculated above. From looking at this plot, the coverage of the confidence interval seems accurate, indicating that there are no problems with the Fisher information standard errors. I can further check the validity of the above confidence interval by looking at the profile log likelihood. Though not included here, this method also gives a confidence interval comparable to the one constructed using the Fisher information standard errors. This lends more credibility to the above analysis.

From this data exploration, it appears that the ARMA(0,0) model, a Gaussian white noise model, is most appropriate for the January low temperature data for Ann Arbor, MI. This is somewhat surprising, given the intuition that temperature might vary systematically from year to year. Further interesting work would be to consider models with trend to see if we can capture some gradual warming. It seems possible, however, that small changes (increases, fluctuations, etc.) could be difficult to detect with such little data on such a long time frame.


### Looking for a trend model


Since this time series is well modeled by white noise, we could fit a signal plus white noise model. This might be a more sensitive way to look for a trend.
First, we try some low-order polynomial trend specifications,

\[ Y_n=\sum_{k=0}^K \beta_k n^k + \epsilon_n\]

where $K$ is the order of the fitted trend model. We compare AIC for $K$ up to 5.


|           | K=0|   K=1|   K=2|   K=3|   K=4|   K=5|
|:----------|---:|-----:|-----:|-----:|-----:|-----:|
|<b>AIC</b> | 861| 862.2| 863.8| 865.4| 865.4| 866.9|

There is still no evidence suggesting anything other than a white noise model. Now, 
As one more attempt, we can compare the Michigan data to the global temperature series.


```
## Error in model.frame.default(formula = low ~ global_temp, data = y, drop.unused.levels = TRUE): invalid type (NULL) for variable 'global_temp'
```


![plot of chunk plot_jan_temp](figure/plot_jan_temp-1.png)

```
## Error in xy.coords(x, y): 'x' and 'y' lengths differ
```

```
## Error in xy.coords(x, y): 'x' and 'y' lengths differ
```

The red dashed line shows 10 times the global annual temperature anomaly (multiplied by 9/5 to move from Celcius to Fahrenheit) compared to the Michigan January low (in Fahrenheit).
The trends appear similar until 1970 or so, after which the global temperature increases while the Michigan January low does not.
However, caution is needed because of the relative scales: A scientifically plausible model probably can't have a coefficient much bigger than 1 degree centigrade in Michigan per degree of global mean temperature.
Given the size of the standard error resulting from year-to-year fluctuations, an effect of this size will be hard to see even if the model is good.
The blue dotted line shows the global climate anomaly in degrees Fahrenheit.
From this perspective, we can be skeptical about whether the apparent pattern (with warm January in the 1940s, cooler in the 1980, and currently reentering a cool phase) could be related to global temperature fluctations and trend.
Interpreting the Michigan data as trend may well be just reading patterns into noise.

------------

**<big>Sources</big>**.

A few people submitted homework assignments close to previously posted homework solutions, and some of these also had improper references to the source. Proper use of sources becomes increasingly important as we move toward midterm and final projects, for which we want to learn from past projects while acknowledging sources accurately. If you find yourself studying an online source closely for your solution, you should put thought into making sure that your report goes intellectually beyond the source - simple paraphrasing is not enough to demonstrate your own contribution. Developing your own original analysis usually involves a fair amount of time writing your own code to implement your own ideas.


**<big>References</big>**.


