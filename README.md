# wgof_torus
Two-sample goodness-of-fit tests on the two-dimensional flat torus based on Wasserstein distance.

The goal of wgof_torus is to provide some practical approaches to perform two-sample goodness-of-fit tests for measures supported on the two-dimensional flat torus based on Wasserstein distance. These techniques have been introduced in https://arxiv.org/abs/2108.00165 and consist on:

1. Testing the equality of the one-dimensional projected distributions into Ng closed geodesics (file: twosample.geodesic.torus.test.R),
2. Upper bound p-values (file: twosample.ubound.torus.test.R).

Each file includes details on how to implement the test, as well as some minimal reproducible examples. A Wasserstein two-sample goodness-of-fit test for measures supported on the circle is also included, in the file twosample.test.s1.R.
