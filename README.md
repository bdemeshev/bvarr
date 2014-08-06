bvarr
=====

The package `bvarr` may be useful for estimating bayesian VAR models.
This is a translation to R of the matlab code by
[Gary Koops and Dimitris Korobilis](http://personal.strath.ac.uk/gary.koop/bayes_matlab_code_by_koop_and_korobilis.html)
with minor improvements.

You may install the package usinge the commands:
```R
library("devtools")
install_github("bdemeshev/bvarr")
```

Very basic example of usage:
```R
library("bvarr")
data(Yraw)
model <- bvar(Yraw,prior = "independent")
bvar.imp.plot(model)
```

Todo:
1. Describe!!! Good vignette is a first need.
2. Check case with no constant!
3. Check case with exogeneous variables!
4. Check border case with M=1!
5. Move hyperparameters in the definition of the function in the
unified fashion.
Something like: minnesota.hyper=c(...), ssvs.hyper=c(....), etc.
6. Make lovely function for prediction. Two types of prediction? 
7. Rewrite code using Rcpp.
8. Add steady-state VAR, see old BMR package vignette


