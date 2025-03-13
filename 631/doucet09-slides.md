---
title: "Doucet & Johansen (2009)"
date: "Mar 13, 2025"
output:
  ioslides_presentation:
    smaller: no
    widescreen: true
    transition: "faster" 
---

## Impact

* Cited 2800 times

* Updates Doucet, de Freitas & Gordon, eds. (2001) "Sequential Monte Carlo methods in practice" (10400 cites)

* A practical introduction is: Arulampalam, Maskell, Gordon, & Clapp (2002) "A tutorial on particle filters for online nonlinear/non-Gaussian Bayesian tracking" (16300 cites)



## Plug-and-play property

* a.k.a. simulation based, likelihood-free, equation free

* The basic particle filter (sec. 4.1, with q=f) is p&p

* engineering applications may emphasize speed over model fidelity


## Particle depletion and the curse of dimensionality

* Why does the sample size needed for importance sampling increase exponentially with the dimension?

## Relationship between this review and the pomp package

* pomp has _de facto_ focused on plug-and-play

* Rao-Blackwellization can make things much quicker, e.g., for a self-driving car, for which filtering failures can [lead to tragic consequences](https://arxiv.org/abs/2409.17380)






