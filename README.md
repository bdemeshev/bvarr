bvarr
=====

The package `bvarr` may be useful for estimating bayesian VAR models.
This is a translation to R of the matlab code by [Gary Koops and Dimitris Korobilis](http://personal.strath.ac.uk/gary.koop/bayes_matlab_code_by_koop_and_korobilis.html).

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
