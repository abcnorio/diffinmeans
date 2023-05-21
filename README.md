# R code implementation of GL Bretthorst (1993) "On the difference in means"

## Overview

Bretthorst paper from 1993 is an analytical solution to the Behrens-Fisher problem of same/ different means/standard deviations from a Bayesian perspective. It works without MCMC like R packages BESTmcmc or brms. Two implementations - very similar - are available.


```
> 10^(307:309)
[1] 1e+307 1e+308    Inf
> 10^(-323:-324)
[1] 9.881313e-324  0.000000e+00
```

## License

- Copyright of the original Mathematica code of UM Studer is with him.
- Copyright of the original Mathematica code of P Gregory is with him.
- R code adopted from Mathematica code: GPL >= v3

## Files

- DiM_Bretthorst_UMS.r = R code completely rewritten as an adaptation of a Mathematica script written by UM Studer (1998).
- DiM_Bretthorst_PG.r = R code completely rewritten as an adapatation of a Mathematica script originally created by P Gregory (2005).

## Functions used

- 

## Function calls

Functions are called identical to the original version. Use as normal but add .brob to the name of the function call and use a function that gives out a Brobdingnag object that acts as input for the integration function.

```
source(".r")
[...]
```

e.g.

```

```

## Output

The textual output is written as a hommage to the Bretthorst paper and mirrors the way results are reported there. Originally this was invented by UM Studer in his Mathematica script.

## Variants

### Qualitative version

The qualitative version is only available for the UMS implementation. It allows to compare "4 of 5 vs. 8 of 10" regarding same/different means/ standard deviations. The function XXX calculates means and standard deviations from such qualitative success/failure input values. Then the usual formula are used to calculate posterior probabilities.

Example

```

```

### Quantitative version

The quantitative version is available for the UMS as well as the PG version. One has to add lower and upper bounds for means and standard deviations as priors besides the actual empirical means and standard deviations that should be compared.

Example:

```

```

## Usage of R package Brobdingnag

The calculation of the posterior probabilities depend on Gamma function calls used for integration. With larger sample sizes (not the ones used in the original paper) this leads to numbers beyond normal computer capability. If that happens, the scripts contain a switch to use the R package Brobdingnag for very large numbers. This works but leads to a substantial loss in speed although a multi-threaded version of Simpson rule is used (but written in pure R).

Example:

```

```

### Graphical output

There is a graphical output for the qualitative version as well as the quantitative one. The qualitative graphical output is based on UMS whereas the quantitative one is based on PG.

Example qutalitative version:

```

```

Example quantitative version:

```

```

It can be that the graphical output for the quantitative version is not possible due to difficult ranges. In such a case the function XXX allows to redefine the range which will allow to create a graphical output. This happens mostly if large numbers require the usage of the R package Brobdingnag and for the ratio of the standard deviations. The following example shows how to handle it.


```

```

## TODO

Rewrite the integration routines to work with R package Rmpfr to allow for arbitrary large numbers and to gain more speed.

## References

Bretthorst, G.L. (1993). On the difference in means.
Gregory, P. (2005). Bayesian logical data analysis for the physical sciences. A comparative approach with Mathematica support. Cambridge University Press. https://www.cambridge.org/nl/academic/subjects/statistics-probability/statistics-physical-sciences-and-engineering/bayesian-logical-data-analysis-physical-sciences-comparative-approach-mathematica-support?format=PB
https://www.cambridge.org/nl/academic/subjects/statistics-probability/statistics-physical-sciences-and-engineering/bayesian-logical-data-analysis-physical-sciences-comparative-approach-mathematica-support?format=PB
Studer, U.M. (1998). Verlangen, Süchtigkeit und Tiefensystemik. Fallstudie des Suchttherapiezentrums für Drogensüchtige start again in Männedorf und Zürich von 1992 bis 1998. Bericht an das Bundesamt für Justiz. Zürich. URL: https://www.bj.admin.ch/dam/data/bj/sicherheit/smv/modellversuche/evaluationsberichte/37.pdf

## R version

All R scripts should work under R >=v3.

## Disclaimer

R code was tested but is provided "as is".


