STAT547O - Bootstrap notes
================
Matias Salibian-Barrera
2019-11-21

#### LICENSE

These notes are released under the “Creative Commons
Attribution-ShareAlike 4.0 International” license. See the
**human-readable version**
[here](https://creativecommons.org/licenses/by-sa/4.0/) and the **real
thing**
[here](https://creativecommons.org/licenses/by-sa/4.0/legalcode).

# DRAFT (Read at your own risk)

## The Fast and Robust Bootstrap

Before we describe the fast and robust bootstrap (FRB) we briefly
illustrate the use of the non-parametric bootstrap (Efron 1979) with a
less-trivial example than the ones often used in tutorials.

Consider the problem of constructing a confidence region (or confidence
intervals) for the eigenvalues of the covariance matrix of a random
vector. To simplify the presentation we consider here a bi-variate case
`p = 2`). In order to construct these confidence sets we need to
estimate the (joint) distribution of the eigenvalue estimators. We will
use non-parametric boostrap to do this.

We first generate a random sample of `n = 200` observations from a
2-dimensional Gaussian random vector with mean 0. The true covariance
matrix has eigenvalues 12 and 3, and we randomly pick the corresponding
eigenvectors (from a uniform distribution on the unit ball). If `U` is
the matrix containing the eigenvectors as columns (so that `U'U = I`),
then the covariance matrix of the bi-variate random vector is `\Sigma =
U \Lambda U^T`, where `\Lambda = diag(12, 3)`.

``` r
p <- 2
# create a population covariance matrix for the example
set.seed(123456)
u <- qr.Q(qr(matrix(rnorm(p*p), p, p)))
la <- c(12, 3) # c(3, 1)
( sigma <- u %*% diag(la) %*% t(u) )
```

    ##           [,1]      [,2]
    ## [1,] 11.110840 -2.685487
    ## [2,] -2.685487  3.889160

We now take the random sample of `n = 200` observations:

``` r
n <- 200
set.seed(123456)
x <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma=sigma)
```

The estimated eigenvalues of the true covariance matrix `\Sigma` are

``` r
( la.hat <- svd(cov(x))$d )
```

    ## [1] 11.69325  3.32978

To construct 95% confidence intervals for the each of the true
eigenvalues, or build a 95% confidence region for both of them
simultaneously we need to estimate the distribution of the random vector
`(\hat{\lambda}_1, \hat{\lambda}_2)^T`. The non-parametric bootstrap
(Efron 1979) estimator for this distribution. is the “empirical” (Monte
Carlo) distribution obtained by re-computing the estimator on `B`
samples obtained from the original sample, with replacement. In this
example we will use `B = 500` bootstrap samples.

For each sample we will compute the eigenvalues of the corresponding
sample covariance matrix. The empirical distribution of these 500
bi-variate vectors constitutes the nonparametric bootstrap estimate of
the distribution of interest. We will save these 500 vectors in the
matrix `la.hat.b` below. The main bootstrap loop is below. In each step
we draw a sample from the original sample and compute the eigenvalues of
the sample covariance matrix.

``` r
B <- 500
la.hat.b <- array(0, dim=c(B, p))
set.seed(123456)
for(j in 1:B) {
  ii <- sample(n, repl=TRUE)
  la.hat.b[j, ] <- svd( cov(x[ii,]) )$d
}
```

The plot below displays the estimated distribution of the vector of
eigenvalues etimators. The red point is the center of the estimated
distribution, and the blue point is the original vector of estimated
eigenvalues.

``` r
plot(la.hat.b, pch=19, col='gray')
la.mu <- colMeans(la.hat.b)
points(la.mu[1], la.mu[2], col='red', pch=19, cex=1.2)
points(la.hat[1], la.hat[2], col='blue', pch=19, cex=1.2)
```

![](Bootstrap_files/figure-gfm/bootsam-1.png)<!-- -->

This distribution looks fairly elliptical (and Gaussian). We can draw a
95% confidence region based on this observation (the green ellipse
below). Note that the the true vector of eigenvalues (indicated with a
light blue point on the plot) falls inside the confidence region.

``` r
aa <- var(la.hat.b) 
plot(la.hat.b, pch=19, col='gray')
xx <- ellipse::ellipse(aa, centre=la.mu)
lines(xx, cex=.7, col='darkgreen', lwd=3)
points(la.mu[1], la.mu[2], col='red', pch=19, cex=1.2)
points(la[1], la[2], col='steelblue', pch=19, cex=1.5)
```

![](Bootstrap_files/figure-gfm/bootsam2-1.png)<!-- -->

To build marginal confidence intervals, we use quantiles of the
estimated marginal distributions. For the first eigenvalue:

``` r
qs <- as.numeric( quantile(la.hat.b[,1] - la.mu[1], c(.025, .975)) )
c(la.mu[1] - qs[2], la.mu[1] - qs[1])
```

    ## [1]  9.855476 13.550389

and for the second eigenvalue:

``` r
qs <- as.numeric( quantile(la.hat.b[,2] - la.mu[2], c(.025, .975)) )
c(la.mu[2] - qs[2], la.mu[2] - qs[1])
```

    ## [1] 2.577254 4.027706

<!-- Another way to look at it, using boxplots: -->

<!-- ```{r marg2} -->

<!-- boxplot(la.hat.b, col='gray80') -->

<!-- points(1, la[1], pch=19, cex=1.75, col='steelblue') -->

<!-- points(2, la[2], pch=19, cex=1.75, col='steelblue') -->

<!-- ``` -->

<!-- **Sanity check**: we do know the asymptotic distribution of each -->

<!-- eigenvalue estimator ( `N(\lambda_j, 2*\lambda_j^2 / n)` ).  -->

<!-- <!-- ```{r sanity} -->

–\> <!-- <!-- diag( var(la.hat.b) ) --> –\> <!-- <!-- 2*la.hat^2/n -->
–\> <!-- <!-- ``` --> –\>
<!-- The corresponding asymptotic CIs are: --> <!-- ```{r asym1} -->
<!-- c(la.hat[1] - qnorm(.975)*sqrt(2*la.hat[1]^2/n),  -->
<!--   la.hat[1] + qnorm(.975)*sqrt(2*la.hat[1]^2/n) ) -->

<!-- c(la.hat[2] - qnorm(.975)*sqrt(2*la.hat[2]^2/n),  -->

<!--   la.hat[2] + qnorm(.975)*sqrt(2*la.hat[2]^2/n) ) -->

<!-- ``` -->

### Fast and Robust Bootstrap for an M-estimator of location

When bootstrapping robust estimators with data that may contain ayptical
observations two main issues arise: (a) re-computing the robust
estimator several hundred times may be computationally too constly (or
unfeasible); and (b) some of the bootstrap samples may contain a larger
proportion of outliers than in the original sample, and thus the
bootstrapped estimators may be heavily influenced by them.

The Fast and Robust Bootstrap avoids both problems: it is very fast (for
each bootstrap sample we only need to compute a weighted average or
weighted least squares estimator), and outliers will be downweighted and
thus will not affect the bootstrapped estimators.

We described this method in class. Here is an application to the
simplest M-estimator of location (where the residual scale is assumed
known).

We first generate a simple data set, with 10% of oultiers:

``` r
set.seed(123456)
n <- 50
ep <- .1
n0 <- floor(n*(1-ep))
x <- c(rnorm(n0), rnorm(n-n0, mean=5, sd=.1))
```

The true location parameter is 0. The sample mean is `{r dr2 mean(x)` We
now compute a robust M-estimator using a bisquare score function (tuned
to achieve 90% efficiency if no outliers were present in the data).

``` r
si.hat <- mad(x) 
max.it <- 100
mu.hat <- median(x)
cc0 <- RobStatTM::lmrobdet.control(family='bisquare', efficiency=.95)$tuning.psi
rhoprime <- function(r, family='bisquare', cc=cc0) {
return(RobStatTM::rhoprime(r, family=family, cc=cc) )
}
tol <- 1e-7
mu.old <- mu.hat + 10*tol
j <- 0
while( (j  < max.it) & (abs(mu.old - mu.hat) > tol) ) {
  mu.old <- mu.hat
  re <- (x - mu.hat) / si.hat
  w <- rhoprime(re) / re
  mu.hat <- sum(w * x ) / sum(w)
  j <- j + 1
}
```

The M-estimator of location is 0.29. We now compute the correction
factor to account for the reduction in variability that results from the
one-step approximation.

``` r
re <- (x - mu.hat)/si.hat
w <- rhoprime(re) / re
tmp1 <- RobStatTM::rhoprime(re, family='bisquare', cc=cc0)
tmp2 <- RobStatTM::rhoprime2(re, family='bisquare', cc=cc0)
wprime <- (tmp1 - tmp2*re)/re^2
swp <- sum(wprime)
wwprime <- ( wprime * sum(w) - wprime * swp ) / ( sum(w) )^2
corr.f <- 1 / ( 1 - sum( wwprime * x ) ) 
```

We now compute the 500 bootstrapped M-estimators. Note that each of them
only involves computing a weighted average of the form
`sum(w*x)/sum(w)`.

``` r
B <- 500
mu.hat.b <- vector('numeric', B)
set.seed(123456)
for(j in 1:B) {
  ii <- sample(n, repl=TRUE)
  mu.hat.b[j] <- sum( w[ii] * x[ii] ) / sum( w[ii] )
}
```

Next we use the quantiles of the centered and corrected bootstrap
distribution to construct the 95% confidence interval for the true
location parameter:

``` r
# boxplot( corr.f*(mu.hat.b - mu.hat) )
uu <- quantile( corr.f*(mu.hat.b - mu.hat), c(.025, .975))
c( mu.hat - uu[2], mu.hat - uu[1])
```

    ##      97.5%       2.5% 
    ## -0.1417883  0.6879522

Note that the classical 95% confidence interval based on the sample mean
fails to contain the true value (zero) and that a 95% confidence
interval based on the median would require an estimate of the density of
the error distribution.

``` r
c(mean(x) - qnorm(.975) * sd(x)/sqrt(n), 
mean(x) + qnorm(.975) * sd(x)/sqrt(n))
```

    ## [1] 0.1593878 1.1368002

``` r
names( summary(quantreg::rq(x ~ 1)) )
```

    ## [1] "call"         "terms"        "coefficients" "residuals"   
    ## [5] "rdf"          "tau"
