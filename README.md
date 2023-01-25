# torustest
Two-sample goodness-of-fit tests on the two-dimensional flat torus based on Wasserstein distance.

**torustest** is an $`\texttt{R}`$ package implementing the approaches introduced in [1] to perform two-sample goodness-of-fit tests for measures supported on the two-dimensional flat torus, based on Wasserstein distance.

### Installing torustest

**torustest** can be installed using

```
devtools::install_github("https://github.com/gonzalez-delgado/torustest")
```

### Two-sample goodness-of-fit tests

Given two measures $`P,Q\in\mathcal{P}(\mathbb{T}^2)`$, we consider the null hypothesis $`H_0:P=Q`$. The package **torustest** allows the assessment of $`H_0`$ through two different procedures, detailed below.

#### Projection to one-dimensional closed geodesics

The first approach tests the equality of the one-dimensional projected distributions of $`P`$ and $`Q`$ into $N_g$ closed geodesics of $`\mathbb{T}^2`$, which are isomorphic to the circle $`S^1`$. Given a pair of samples in the periodic $`[0,1]\times[0,1]`$, the function [twosample.geodesic.torus.test](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.geodesic.torus.test.R) returns a $`p`$-value for $`H_0`$, after specifying the number of projections $`N_g`$, which may be randomly chosen.





2. Upper bound p-values ([twosample.ubound.torus.test function](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.ubound.torus.test.R) function).

Besides, a Wasserstein two-sample goodness-of-fit test for measures supported on the circle is also performed by [twosample.test.s1.R](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.test.s1.R).


### References

[1] González-Delgado, J., González-Sanz, A., Cortés, J., & Neuvial, P. (2021). Two-sample goodness-of-fit tests on the flat torus based on Wasserstein distance and their relevance to structural biology. arXiv:2108.00165. [[arxiv]](https://arxiv.org/abs/2108.00165)[[HAL]](https://hal.archives-ouvertes.fr/hal-03369795v2).
