# torustest
Two-sample goodness-of-fit tests on the two-dimensional flat torus based on Wasserstein distance.

The goal of **torustest** is to provide some practical approaches to perform two-sample goodness-of-fit tests for measures supported on the two-dimensional flat torus based on Wasserstein distance. These techniques have been introduced in [1] and consist on:

1. Testing the equality of the one-dimensional projected distributions into $N_g$ closed geodesics [(twosample.geodesic.torus.test function)](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.geodesic.torus.test.R),
2. Upper bound p-values [(twosample.ubound.torus.test function)](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.ubound.torus.test.R).

Besides, a Wasserstein two-sample goodness-of-fit test for measures supported on the circle is also performed by [(twosample.test.s1.R)](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.test.s1.R).

Each file includes details on how to implement the test, as well as some minimal reproducible examples.

### Installing torustest

**torustest** can be installed using

```
devtools::install_github("https://github.com/gonzalez-delgado/torustest")
```

### References

[1] González-Delgado, J., González-Sanz, A., Cortés, J., & Neuvial, P. (2021). Two-sample goodness-of-fit tests on the flat torus based on Wasserstein distance and their relevance to structural biology. arXiv:2108.00165. [[arxiv]](https://arxiv.org/abs/2108.00165)[[HAL]](https://hal.archives-ouvertes.fr/hal-03369795v2).
