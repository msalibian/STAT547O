
## R Scripts with some examples

We will generally use the packages `robustbase` and `RobStatTM`.

#### Location M-estimators

The implementation of location M-estimators in RobStatTM does not
include score functions in the Huber family. We can do it
with Tukey's bisquare family.

```
library(RobStatTM)
set.seed(123)
x <- rnorm(30, mean=3, sd=.5)
( tmp <- MLocDis(x, psi='bisquare', eff=.95, tol=1e-6) )
```

We check (sanity check) that the returned values solve
the desired equation. Note that `tmp$disper` is not
`mad(x)` (the residual scale estimator used to compute
the location estimator), so we compute it here:
```
si <- mad(x)
cc <- bisquare(e = .95)
mean( rhoprime(u=(x-tmp$mu)/si, family='bisquare', cc=cc) )
```
What is `tmp$disper` then? It is an M-scale estimator of the residuals
(with respect to the M-location estimator, computed with the mad):
```
cc2 <- lmrobdet.control(bb=.5, family='bisquare')$tuning.chi
mscale(x-tmp$mu, delta=.5, family='bisquare', tuning.chi=cc2)
```
However, this is still not equal to `tmp$disper`! The reason is 
that `MLocDis` does not use the correct tuning constant for a 
50% breakdown point M-scale estimator. Specifically, it
uses `1.56` instead of ``1.547...`:
```
tmp$disper
mscale(x-tmp$mu, delta=.5, family='bisquare', tuning.chi=1.56)
```

With the `robustbase` package we can use the function `huberM` which does
what we want:
```
library(robustbase)
( tmp <- huberM(x=x, k=1.5) )
```
Sanity check of the estimating equations
```
Hrhoprime <- function(a, cc) pmax(pmin(a, cc), -cc)
mean( Hrhoprime((x-tmp$mu)/tmp$s, cc=1.5) )
```
