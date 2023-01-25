# torustest

**torustest** is an $\texttt{R}$ package implementing the approaches introduced in [1] to perform two-sample goodness-of-fit tests for measures supported on the two-dimensional flat torus, based on Wasserstein distance.

### Installing torustest

**torustest** can be installed using

```
devtools::install_github("https://github.com/gonzalez-delgado/torustest")
```

### Two-sample goodness-of-fit tests

Given two measures $P,Q\in\mathcal{P}(\mathbb{T}^2)$, we consider the null hypothesis $H_0:P=Q$. The package **torustest** allows the assessment of $H_0$ through two different procedures, detailed below.

#### Projection into one-dimensional closed geodesics

The first approach tests the equality of the one-dimensional projected distributions of $P$ and $Q$ into $N_g$ closed geodesics of $\mathbb{T}^2$, which are isomorphic to the circle $S^1$. Given a pair of samples in the periodic $[0,1]\times[0,1]$, the function [twosample.geodesic.torus.test](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.geodesic.torus.test.R) returns a $p$-value for $H_0$, after specifying the number of projections $N_g$, which may be randomly chosen.

#### Example

```
n <- 100 # Sample size
 
# Simulate the null distribution of the circle test statistic
sim_free_null <- sim.null.stat(500, NC = 2)

# Bivariate von Mises distributions
samp_1 <- BAMBI::rvmcos(n) / (2 * pi) 
samp_2 <- BAMBI::rvmcos(n) / (2 * pi)

#4 geodesics are chosen randomly
twosample.geodesic.torus.test(samp_1, samp_2, n_geodesics = 3, NC_geodesic = 2, sim_null = sim_free_null) 

#4 geodesics are chosen a priori
glist <- list(c(1, 0), c(0, 1), c(1, 1), c(2, 3))
twosample.geodesic.torus.test(samp_1, samp_2, geodesic_list = glist, NC_geodesic = 2, sim_null = sim_free_null) 

```




2. Upper bound p-values ([twosample.ubound.torus.test function](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.ubound.torus.test.R) function).

Besides, a Wasserstein two-sample goodness-of-fit test for measures supported on the circle is also performed by [twosample.test.s1.R](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.test.s1.R).


### References

[1] González-Delgado, J., González-Sanz, A., Cortés, J., & Neuvial, P. (2021). Two-sample goodness-of-fit tests on the flat torus based on Wasserstein distance and their relevance to structural biology. arXiv:2108.00165. [[arxiv]](https://arxiv.org/abs/2108.00165)[[HAL]](https://hal.archives-ouvertes.fr/hal-03369795v2).
