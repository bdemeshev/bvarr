bvarr
=====

The package `bvarr` may be useful for estimating bayesian VAR models.

* It contains a translation to R of the matlab code by
[Gary Koops and Dimitris Korobilis](http://personal.strath.ac.uk/gary.koop/bayes_matlab_code_by_koop_and_korobilis.html) with minor improvements.
* It has very versatile realisation of conjugate Normal-Inverse Wishart prior.

You may install the package usinge the commands:
```R
install.packages("devtools")
devtools::install_github("bdemeshev/bvarr")
```

Replication of KK code:
```R
library("bvarr")
data(Yraw)
model <- bvar(Yraw,prior = "independent")
bvar.imp.plot(model)
bvar.summary(model)
```

New conjugate Normal-Inverse Wishart function:
```R
library("bvarr")
data(Yraw)
priors <- Carriero_priors(Yraw, p = 4)
model <- bvar_conjugate0(priors = priors)
summary_conjugate(model) 
forecast_conjugate(model, h=2, output="wide")
forecast_conjugate(model, out_of_sample = FALSE, include="mean", level=NULL, type = "credible")
```



Goals of the package:

1. Good documentation

2. Versatile (many-many options)


You may also wish look at [BMR](http://bayes.squarespace.com/bmr/) and [MSBVAR](http://cran.r-project.org/web/packages/MSBVAR/) packages. They seem to be more professional, but they do not contain something I wish for :) Maybe `bvarr` will be redesigned completely in late 2015.


[![Travis-CI Build Status](https://travis-ci.org/bdemeshev/bvarr.svg?branch=master)](https://travis-ci.org/bdemeshev/bvarr)


