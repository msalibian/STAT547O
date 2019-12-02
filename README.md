# A short course on modern Robust Statistics

This repository contains material (notes, reading suggestions and topics) for a short
graduate course on Robust Statistics (Fall 2019).


###### LICENSE
These notes are released under the
"Creative Commons Attribution-ShareAlike 4.0 International" license.
See the **human-readable version** [here](https://creativecommons.org/licenses/by-sa/4.0/)
and the **real thing** [here](https://creativecommons.org/licenses/by-sa/4.0/legalcode).

## Tentative schedule

- Week 1: Oct 22 / 24
    - Introduction, motivation, goals, general setting (atypical observations).
    - Location / scale (**brief**); M estimators; The problem of scale.
        - [Lecture 1 non-technical notes](Lecture1.md)
            - source available as [R Markdown](Lecture1.Rmd).
    - Contamination setting; Breakdown Point; "Bias", Influence Function.
    - Linear regression; M-estimators; The problem of scale.
    - Fixed designs;  Random explanatory variables.  
    - M-estimators of scale; S-estimators of regression.
        - [Lecture 2 non-technical notes](Lecture2.md)
             - source available as [R Markdown](Lecture2.Rmd).
    - Reference: Maronna *et al.* (2018), *Robust Statistics: Theory and Methods, 2nd Edition*, Wiley. [UBC Library link](http://tinyurl.com/wyfryxa)
- Week 2: Oct 29 / 31
    - M-estimators of scale; S-estimators of regression; MM-estimators.
        - References:
            - Maronna *et al.* (2018), *Robust Statistics: Theory and Methods, 2nd Edition*, Wiley. [UBC Library link](http://tinyurl.com/wyfryxa);
            - SB and Zamar (2004) Uniform asymptotics for robust location estimates when the scale is unknown, *The Annals of Statistics*, 32(4):1434-1447, [DOI](https://doi.org/10.1214/009053604000000544);
            - Omelka and SB (2010) Uniform asymptotics for S- and MM-regression estimators, *Annals of the Institute of Statistical Mathematics*, 62(5):897–927, [DOI](https://doi.org/10.1007/s10463-008-0189-x).
    - Non-parametric regression;  Kernel smoothers.
        - References:
             - Boente, Mart&iacute;nez and SB (2017) Robust estimators for
          additive models using backfitting. *Journal of Nonparametric
          Statistics*. 29:744-767:
          [DOI](https://doi.org/10.1080/10485252.2017.1369077);
              - Boente and Fraiman (1989), Robust nonparametric regression estimation,
          *Journal of Multivariate Analysis*, 29(2): 180-198,
          [DOI](https://doi.org/10.1016/0047-259X(89)90023-7);
             - H&auml;rdle.  (1984). Robust regression function estimation, *Journal of Multivariate
         Analysis*,  14(2): 169-180, [DOI](
         https://doi.org/10.1016/0047-259X(84)90003-4);
             <!-- - Welsh, A. (1996). Robust estimation of smooth regression and spread functions and their derivatives. *Statistica Sinica*, 6(2), 347-366. [JSTOR](https://www.jstor.org/stable/24306020); -->
             - H&auml;rdle, M&uuml;ller, Sperlich and Werwatz (2004) *Nonparametric and semiparametric models*, Springer, [UBC link](http://tinyurl.com/rn7z29x).
    - Brief intro to additive models / curse of dimensionality.
        - References:
             Hastie, T.J. and Tibshirani, R.J. (1990) *Generalized Additive Models*. Monographs on Statistics and Applied Probabilitiy, 43. Chapman & Hall, CRC.;
- Week 3: Nov 5 / 7
    - Algorithms for M-, S- and MM-regression estimators.
        - References:
            - SB & Yohai (2006) A Fast Algorithm for S-Regression Estimates, *Journal of Computational and Graphical Statistics*, 15:2, 414-427, [DOI](https://doi.org/10.1198/106186006X113629);
              - Maronna *et al.* (2018), *Robust Statistics: Theory and Methods, 2nd Edition*, Wiley. [UBC Library link](http://tinyurl.com/wyfryxa)
    - Practical / computational digression
        - [Computation by hand](Simple_examples.md)
            - source available as [R Markdown](Simple_examples.Rmd);
            - the corresponding sandbox as a [Jupyter Notebook](Simple_examples.ipynb).
    - Additive models  + Splines.
        - [Backfitting by hand](Example-backfitting.md)
            - source available as [R Markdown](Example-backfitting.Rmd)
            - the corresponding sandbox as a [Jupyter Notebook](Example-backfitting.ipynb).
        - References:
            - Boente, Mart&iacute;nez and SB (2017) Robust estimators for additive models using backfitting. *Journal of Nonparametric Statistics*. 29:744-767: [DOI](https://doi.org/10.1080/10485252.2017.1369077);
            - Tharmaratnam, Claeskens, Croux & SB (2010) S-Estimation for Penalized Regression Splines, *Journal of Computational and Graphical Statistics*, 19:3, 609-625, [DOI](https://doi.org/10.1198/jcgs.2010.08149);
               - Ruppert, Wand and Carroll (2003) *Semiparametric regression*, Cambridge University Press, [UBC link](http://tinyurl.com/qtdua46).
- Week 4: Nov 12 / 14
    - Multivariate analysis
         - Elliptical distributions
            - References:
                - Boente, SB and Tyler,  (2014). A characterization of elliptical distributions and some optimality properties of principal components for functional data,  *Journal of Multivariate Analysis*, 131:254-264, [DOI](https://doi.org/10.1016/j.jmva.2014.07.006);
                - Paindaveine, (2006), Elliptical Symmetry, *Encyclopedia of Environmetrics*, [DOI](https://doi.org/10.1002/9780470057339.vnn081).
         - Estimation; Outlier detection; PCA.
              - References: Maronna *et al.* (2018), *Robust Statistics: Theory and Methods, 2nd Edition*, Wiley. [UBC Library link](http://tinyurl.com/wyfryxa)
        - Functional Data Analysis:
            - References: Boente & SB (2015) S-Estimators for Functional Principal Component Analysis, *Journal of the American Statistical Association*, 110:511, 1100-1111, [DOI](https://doi.org/10.1080/01621459.2014.946991).
- Week 5: Nov 19 / 21
     - Principal Components Analysis (including brief mention to
         methods for Functional Data PCA).
        - References:
            - Boente, SB and Tyler,  (2014). A characterization of elliptical distributions and some optimality properties of principal components for functional data,  *Journal of Multivariate Analysis*, 131:254-264, [DOI](https://doi.org/10.1016/j.jmva.2014.07.006);
             - Boente & SB (2015) S-Estimators for Functional Principal Component Analysis, *Journal of the American Statistical Association*, 110:511, 1100-1111, [DOI](https://doi.org/10.1080/01621459.2014.946991).
    - Inference: asymptotics; bootstrap;
        - References:
            - SB & Zamar (2002) Bootrapping robust estimates of regression, *The Annals of Statistics*, 30(2), 556-582 [DOI](https://doi.org/10.1214/aos/1021379865);
            - SB, Van Aelst, & Willems (2008) Fast and robust bootstrap, *Statistical Methods and Applications*, 17(1): 41-71. [DOI](https://doi.org/10.1007/s10260-007-0048-6);
            - Christmann, SB, Van Aelst. (2013) Qualitative Robustness of Bootstrap Approximations for Kernel Based Methods. *Robustness and Complex Data Structures*, Berlin, Heidelberg: Springer Berlin Heidelberg; 263–278. [DOI](http://dx.doi.org/10.1007/978-3-642-35494-6_16).
    - [Bootstrap and Fast and Robust Bootstrap by hand](Bootstrap.md)
         - source available as [R Markdown](Bootstrap.Rmd);
         - the corresponding sandbox as a [Jupyter Notebook](Bootstrap.ipynb).
- Week 6: Nov 26 / 28
    - Robust Ensemble methods for regression
    - GAM + GAPLM
        - References:
             - Hastie, T.J. and Tibshirani, R.J. (1990) *Generalized Additive Models*. Monographs on Statistics and Applied Probabilitiy, 43. Chapman & Hall, CRC.;
             - Lee, T-Y. (2019), "Robust methods for generalized partial linear partial additive models with an application to detection of disease outbreaks", MSc Thesis, Department of Statistics, UBC, available on-line [here](https://dx.doi.org/10.14288/1.0380711);
            - Alimadad A, S-B, M. (2011), An Outlier-Robust Fit for Generalized Additive Models With Applications to Disease Outbreak Detection, *Journal of the American Statistical Association*, 106:719-731;
             - Boente, Mart&iacute;nez and SB (2017) Robust estimators for additive models using backfitting. *Journal of Nonparametric Statistics*. 29:744-767. [DOI](https://doi.org/10.1080/10485252.2017.1369077);           
- Other topics:
    - Random / Mixed effects models
    - Spatial statistics
    - State-space models
    - Non-standard random elements (phylogenetic trees, etc.)
