# R code implementation of GL Bretthorst (1993) "On the difference in means"

## Overview

Bretthorst (1993)  "On the difference in means" is an analytical solution to the Behrens-Fisher problem of same/ different means/standard deviations from a Bayesian perspective. It works without MCMC like R packages BESTmcmc or brms. Two implementations - very similar - are available.

## Main functions

- `DiM_Bretthorst_UMS.r` = R code as an adaptation of UM Studer (1998)
- `DiM_Bretthorst_PG.r` = R code as an adapatation of P Gregory (2005)

## Functions available

| Function | Description |
| --- | --- |
| `SucRatesIntBounds` | Calculate xbar/sd from successes/failures |
| `DiffinMeans` | calculate difference in means posterior probs after UMS |
| `UMSprint` | print results of `DiffinMeans` |
| `UMSplot` | plot successes/failures |
| | |
| `DiM.pg` | calculate difference in means posterior probs after PG |
| `DiM.print.pg` | print results of `DiM.pg` |
| `DiM.extract.limits` | extract and apply (new) limits for the plots |
| `DiM.plot.calc.pg` | calculate plot points of `DiM.pg` |
| `DiM.plot.pg` | plot results of `DiM.pg` |
| `ums2pg` | internal function to convert UMS successes/failures to PG style input |

## Output

The textual output is written as a hommage to the Bretthorst paper and mirrors the way results are reported there. Originally this was invented by UM Studer (1998) in his Mathematica script.

## Function calls and variants

### Qualitative version

The qualitative version is only available for the UMS implementation. It allows to compare e.g. "4 of 5 vs. 8 of 10" regarding same/different means/ standard deviations. The function `SucRatesIntBounds` calculates means and standard deviations from such qualitative successes/failures input values. Then the usual formulas are used to calculate posterior probabilities.

Example

Load scripts into R

```
source("DiM_Bretthorst_PG.r")
source("DiM_Bretthorst_UMS.r")
```
and call (successes/failures)
```
# calculate means/ standard deviations + lower/ upperbounds from successes/ failures
res.SIB.NRFmf <- SucRatesIntBounds(Si=20, Ni=47, Sii=13, Nii=28, smin=0, snames=c("male","female"))
# calculate posterior porbabilities of the difference in means (UMS version)
DiM.res.NRFmf <- DiffinMeans(inval=res.SIB.NRFmf, out=FALSE)
# output results
UMSprint(results=DiM.res.NRFmf)
```
The basic structure of the input values for DiffinMeans is:

```
# GENERAL call and preparation
# Studer 1998, p.47
# Bretthorst, 1993, p.189 (from Jaynes, 1976 + 1983)	
inval <- list(
		Si=NULL,		#UMS specific -> successes group1
		Ni = 4,			#N group1
		Sii=NULL,		#UMS specific -> successes group2
		Nii = 9,		#N group2
		smin = 0,		#UMS specific -> bounds on the mean -> only for SucRatesIntBounds()
		Di = 50,		#mean group1
		si = 6.48,		#sd group1
		Dii = 42,		#mean group2
		sii = 7.48,		#sd group2
		L = 34,			#mean lower bound
		H = 58,			#mean upper bound
		sL = 3,			#variance lower bound
		sH = 10,		#variance upper bound
		snames = c("Jaynes.1","Jaynes.2")
)
```
The UMS _specific_ values represent the qualitative version of UMS (comparison of two groups with respect to successes/ failures in each group). If those are `NULL` you can insert quantitative values ie. means and standard deviations plus lower/ upper bounds as priors for means and standard deviations (see example above by inspecting `res.SIB.NRFmf`).

### Quantitative version

The quantitative version is available for the UMS as well as for the PG version. One has to add lower and upper bounds for means and standard deviations as priors. Additionally, the actual empirical means and standard deviations of the two sets to be compared have to be announced as input values along of which type they are (`invtyp` can be UMS `ums` or PG `pg` scheme).

Example:

```
# call PG scheme
inputvalues <- list(snames = c("riverB.1","riverB.2"),
                   # sample 1
                   d1 = c(13.2,13.8,8.7,9,8.6,9.9,14.2,9.7,10.7,8.3,8.5,9.2),
                   # sample 2
                   d2 = c(8.9,9.1,8.3,6,7.7,9.9,9.9,8.9),
                    
                   # Input priors and no. of steps in evaluation of p(r|D_1,D_2,I) & p(delta|D_1,D_2,I)
                   # ndelta = number of steps in delta parameter (mean difference)
                   ndelta = 1000
                   # nr = number of steps in r parameter (ratio of the standard deviations)
                   nr = 1000
                    
                   # Set prior limits on mean and standard deviation
		   # assumed that the data set is the same
                   # upper mean
                   high = 12,
                   # lower mean
                   low = 7,
                   # upper sd
                   sigma.high = 4,
                   # lower sd
                   sigma.low = 1
		   )
# calculate posterior porbabilities of the difference in means (PG version)
dim.res <- DiM.pg(invtyp="pg", inputvalues, print.res=FALSE)
# output results
DiM.print.pg(dim.res)
# extract limits for plot
DiM.extract.limits(dim.res.newlimits, change=FALSE)
# cacöulate plot values
dim.res.calc <- DiM.plot.calc.pg(dim.res.newlimits, BROB=FALSE)
# plot results
DiM.plot.pg(dim.res.calc, filling=FALSE, BROB=FALSE)
```

One can convert input values to be used within the other script.

Example:

```
# calculate means/ standard deviations + lower/ upperbounds from successes/ failures 
res.SIB <- SucRatesIntBounds(Si=11, Ni=15, Sii=10, Nii=16, smin=0, snames=c("voluntary","non-voluntary"))
# use input values after UMS for PG version
DIM.pg.res <- DiM.pg(invtyp="ums", inputvalues=res.SIB, print.res=TRUE)
```

or

```
# insert successes/ failures manually and let the script calculate the rest
inputvalues.UMS <- list(snames=c("Jaynes.1","Jaynes.2"), si=6.48, Ni=4, sii=7.48, Nii=9, Di=50, Dii=42, L=34, H=58, sL=3, sH=10, ndelta=1000, nr=1000)
dim.res <- DiM.pg(invtyp="ums", inputvalues=inputvalues.UMS, print.res=TRUE)
```

## Usage of R package [Brobdingnag](https://github.com/RobinHankin/Brobdingnag)

The calculation of the posterior probabilities depend on Gamma function calls used for integration. With larger sample sizes (not the ones used in the original papers and books) this leads temporarily to numbers beyond normal computer capability. If that happens, the scripts contain a switch to use the R package [Brobdingnag](https://github.com/RobinHankin/Brobdingnag) for very large numbers. This works well but leads to a substantial loss in speed. A multi-threaded version, compiled code (just R code here) or a different programming language for this part could speed up the process.

Example output that fails due to very large numbers:

```
> res.SIB.NRFtotal <- SucRatesIntBounds(Si=(20+13), Ni=(47+28), Sii=(338 %/% 4), Nii=338, smin=0, snames=c("male","female"))
> DiM.res.NRF <- DiffinMeans(inval=res.SIB.NRFtotal, out=FALSE)

L - Mean_comb < 0 :  TRUE  [comparison L < DD]
'+'-sign between Gamma-factors is ok


Calculate PMV
Error in integrate(integpmv, lower = sL, upper = sH) :
non-finite function value
```

Look at the plot - it becomes obvious WHY large numbers occur, but it should work nevertheless.
```
UMSplot(inval=res.SIB.NRFtotal,pdfout=FALSE)
```

So we switch to `BROB=TRUE`:

```
> DiM.res.NRFmf.brob <- DiffinMeans(inval=res.SIB.NRFmf, out=TRUE, BROB=TRUE)

Using log() and package 'Brobdingnag'
All results are expressed as log(RESULT)

L - Mean_comb < 0 :  TRUE  [comparison L < DD]
'+'-sign between Gamma-factors is ok


Calculate PMV
Calculate PMbarV
Calculate PMVbar
Calculate PMbarVbar
Compile results
Create dataframes for output

Short output of 'the difference in means' based on BROBs:

 No. Standard Deviation Mean      Data set
 47  0.06998542         0.4285714 male    
 28  0.08960287         0.4666667 female  
 75  0.07948694         0.4427937 combined

 Numerical Example                    Value
 Prior Mean lower bound               0.020
 Prior Mean upper bound               0.980
 Prior Standard Deviation lower bound 0.019
 Prior Standard Deviation upper bound 0.090

 Hypothesis (Probabilities)                     abbrev    p(H|D1,D2,I) = exp(x) sign
 Same Mean,       same Standard Deviation       C&S       -0.6754089            TRUE
 Different Means, same Standard Deviation       Cbar,S    -1.5950066            TRUE
 Same Mean,       different Standard Deviations C,Sbar    -1.5786703            TRUE
 Different Means, different Standard Deviations Cbar&Sbar -2.5023128            TRUE
 The Means are the same                         C         -0.3351967            TRUE
 The Means are different                        Cbar      -1.2559591            TRUE
 The Standard Deviations are the same           S         -0.3398804            TRUE
 The Standard Deviations are different          Sbar      -1.2442930            TRUE
 The Data Sets are the same                     C&S       -0.6754089            TRUE
 The Data Sets are different                    Cbar|Sbar -0.7112058            TRUE

 Hypothesis in favor of ...                                             Odds Ratio (OR) = exp(x) sign
 A difference in Means                                                  -0.92076242              TRUE
 The same Means                                                          0.92076242              TRUE
 A difference in Standard Deviations                                    -0.90441264              TRUE
 The same Standard Deviations                                            0.90441264              TRUE
 A difference in the Sets (different Means and/ or Standard Deviations) -0.03579686              TRUE
 The same Sets (same Means and/ or Standard Deviations)                  0.03579686              TRUE
```
Be aware that the output contains now the sign of the results exp(x) like Brobdingnag stores it.

### Graphical output

There is a graphical output for the qualitative version as well as the quantitative one. The qualitative graphical output is based on the UMS version whereas the quantitative one is based on the PG version.

Example qualitative version:

```
UMSplot(inval=res.SIB.NRFmf,pdfout=FALSE)
```

Example quantitative version:

```
# ratio of SD requires input values scaleL (low) and scaleH (high), otherwise the script breaks
DiM.plotvalues.res.nonbrob <- DiM.plot.calc.pg(DIM.pg.res, scaleL=2, scaleH=8)
DiM.plot.pg(DiM.plotvalues.res.nonbrob, filling=TRUE)
```

It can happen that the graphical output for the quantitative version is not possible due to difficult ranges/ limits chosen. In such a case the function `DiM.extract.limits` with the switch `change=TRUE` allows to redefine the range which will allow to create a graphical output. This happens mostly for large numbers ie. larger sample sizes that require the usage of the R package Brobdingnag and esp. for the graphical output of the ratio of the standard deviations. This requires a little bit of experimentation to get the best result. The following example shows how to handle it.


```
# this example has very large numbers
DIM.pg.res.brob <- DiM.pg(invtyp="ums", inputvalues=res.SIB.NRFtotal, print.res=TRUE, BROB=TRUE)
# we can try out several limits before applying it
DiM.extract.limits(DIM.pg.res.brob, scaleL=5, scaleH=2, change=FALSE)
DiM.extract.limits(DIM.pg.res.brob, scaleL=10, scaleH=3, change=FALSE)
# apply new limits
DiM.newlimits <- DiM.extract.limits(DIM.pg.res.brob, scaleL=5, scaleH=2, change=TRUE)
# check
DiM.extract.limits(DiM.newlimits, change=FALSE)
# calculate plot values - this requires time our case here
DiM.newlimits.calc.plot <- DiM.plot.calc.pg(DiM.newlimits, BROB=TRUE)
# plot results as usual
DiM.plot.pg(DiM.newlimits.calc.plot, filling=FALSE, by1=TRUE, BROB=TRUE)
```

## TODO

Rewrite the integration function used to work multi-threaded and speed up the calculations. Us `logSum` and some other commands to speed up the calculation of the integrals. Another alternative is to use the package Rmpfr is.

## References

Bretthorst, G.L. (1993). [On the difference in means](https://bayes.wustl.edu/glb/diff.pdf). In _Physics & Probability Essays in honor of Edwin T. Jaynes_, W. T. Grandy and P. W. Milonni Eds., Cambridge University Press, England.

Gregory, P. (2005). [_Bayesian logical data analysis for the physical sciences_]( https://www.cambridge.org/nl/academic/subjects/statistics-probability/statistics-physical-sciences-and-engineering/bayesian-logical-data-analysis-physical-sciences-comparative-approach-mathematica-support?format=PB). A comparative approach with [Mathematica support](https://www.cambridge.org/nl/academic/subjects/statistics-probability/statistics-physical-sciences-and-engineering/bayesian-logical-data-analysis-physical-sciences-comparative-approach-mathematica-support?format=PB). Cambridge University Press.

Studer, U.M. (1998). [_Verlangen, Süchtigkeit und Tiefensystemik. Fallstudie des Suchttherapiezentrums für Drogensüchtige start again in Männedorf und Zürich von 1992 bis 1998_](https://www.bj.admin.ch/dam/data/bj/sicherheit/smv/modellversuche/evaluationsberichte/37.pdf). Bericht an das Bundesamt für Justiz (BAJ). Zürich.

## License and Credits

The R code is based on two Mathematica scripts as role models - very similar, but not identical implementations.

- (C) by UM Studer (1998). This Mathematica code is unpublished.
- (C) by P Gregory (2005). You can get his Mathematica notebook for free on the [book](https://www.cambridge.org/nl/academic/subjects/statistics-probability/statistics-physical-sciences-and-engineering/bayesian-logical-data-analysis-physical-sciences-comparative-approach-mathematica-support?format=PB)'s website.
- R code: GPL >= v3

## R version

All R scripts should work under R >=v3.

## Disclaimer

The R code was tested carefully, and cross-checked against various publication results to ensure proper results. However, it is provided "as is". Use common sense to compare results with expectations. NO WARRANTY of any kind is involved here. There is no guarantee that the software is free of error or consistent with any standards or even meets your requirements. Do not use the software or rely on it to solve problems if incorrect results may lead to hurting or injurying living beings of any kind or if it can lead to loss of property or any other possible damage to the world. If you use the software in such a manner, you are on your own and it is your own risk.
