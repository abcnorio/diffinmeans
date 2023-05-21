################################################################################
# ON THE DIFFERENCE OF MEANS
################################################################################

# original Mathematica code by Phil Gregory
# this was used as a role model to write the R code

# Ch. 9: Bayesian analysis of two independent samples
# Introduction
# This is a Mathematica implementation of the probability calculations discussed in the book in the section entitled,     
# "Bayesian Comparison of Two Samples?". 
#
# based on G.L. Bretthorst (1993) "On the difference of means"

#
# see also Mathematica code by UM Studer (1996 + 1998) on the same problem (paper)
# taken code from there to allow calculations based only on summary statistics
# and not on raw data (see also Bretthorst, 1993, for an example taken from Jaynes)

# R code by Leo G 2017
# first = 2017-04-19
# last = 2017-27-04
# last = 2024-05-21

# notes:
# *- introduce logs to integral calculations, but probably that won't help...
# *- very small numbers are slightly different from Mathematica -> e.g. e-230

# additional note:
# *- we use the sdame variable names and labels like P Gregory to allow comparison
#    with his Mathematica notebook. This is different from the UMS script where the
#    labels from UMS are used.

################################################################################
# pre-calculations taken from UMS 1998
ums2pg <- function(inputvalues)
#convert values from UMS implementation to PG implementation
{

  #mapping
  #n
  #NN <- Ni+Nii
  #
  #dbar
  #DD <- (Ni * Di + Nii * Dii) / NN
  #
  #d1squbar
  #Dsi <- (Ni-1) / Ni * si^2 + Di^2
  #
  #d2squbar
  #Dsii <- (Nii-1) / Nii * sii^2 + Dii^2
  #
  #dsqubar
  #DsD <- (Ni * Dsi + Nii * Dsii) / NN
  #
  #dstd
  #ss <- sqrt(NN * (DsD - DD^2) / (NN-1))

  d1bar <- inputvalues[["Di"]]
  d2bar <- inputvalues[["Dii"]]
  d1std <- inputvalues[["si"]]
  d2std <- inputvalues[["sii"]]
  n1 <- inputvalues[["Ni"]]
  n2 <- inputvalues[["Nii"]]
  low <- inputvalues[["L"]]
  high <- inputvalues[["H"]]
  sigma.low <- inputvalues[["sL"]]
  sigma.high <- inputvalues[["sH"]]
  snames <- inputvalues[["snames"]]

  if(is.null(inputvalues[["ndelta"]])) ndelta <- 1000 else ndelta <- inputvalues[["ndelta"]]
  if(is.null(inputvalues[["nr"]])) nr <- 1000 else nr <- inputvalues[["nr"]]

  # manual update of the data
  n <- n1 + n2
  dbar <- (d1bar*n1 + d2bar*n2)/n#(Ni+Nii)
  # Compute the mean square average of the 1st data set
  d1squbar <- (n1-1) / n1 * d1std^2 + d1bar^2
  # Compute the mean square average of the 2nd data set
  d2squbar <- (n2-1) / n2 * d2std^2 + d2bar^2
  # Compute the mean square average of the combined data set
  dsqubar <- (n1 * d1squbar + n2 * d2squbar) / n
  # Compute the standard deviation of the combined data set
  dstd <- sqrt(n * (dsqubar - dbar^2) / (n-1))

  res <- list(snames=snames,
              n1=n1, n2=n2, n=n,
              low=low, high=high, sigma.low=sigma.low, sigma.high=sigma.high,
              d1bar=d1bar, d2bar=d2bar, dbar=dbar,
              d1std=d1std, d2std=d2std, dstd=dstd,
              d1squbar=d1squbar, d2squbar=d2squbar, dsqubar=dsqubar,
              ndelta=ndelta, nr=nr)
  return(res)
}
# call:
#
################################################################################


################################################################################
DiM.pg <- function(invtyp=NULL, inputvalues=NULL, print.res=TRUE, dig=4, Nsteps=100, BROB=FALSE, ndelta=1000, nr=1000)
{
  #Bretthorst DiM
  #adapted from Mathematica code by Phil Gregory

  #for large numbers
  if(BROB)
  {
    library(Brobdingnag)
    cat("\nUsing log() and package 'Brobdingnag'\nAll results are expressed as log(RESULT)\n")
  }    

####################
# START INPUT VALUES
 if(invtyp=="ums")
 {
  d1 <- NULL
  d2 <- NULL
 
  ums2pg.res <- ums2pg(inputvalues)
  snames <- ums2pg.res[["snames"]]

  n1 <- ums2pg.res[["n1"]]
  n2 <- ums2pg.res[["n2"]]
  n <- ums2pg.res[["n"]]

  low <- ums2pg.res[["low"]]
  high <- ums2pg.res[["high"]]
  sigma.low <- ums2pg.res[["sigma.low"]]
  sigma.high <- ums2pg.res[["sigma.high"]]

  d1bar <- ums2pg.res[["d1bar"]]
  d2bar <- ums2pg.res[["d2bar"]]
  dbar <- ums2pg.res[["dbar"]]
 
  d1std <- ums2pg.res[["d1std"]]
  d2std <- ums2pg.res[["d2std"]]
  dstd <- ums2pg.res[["dstd"]]
 
  d1squbar <- ums2pg.res[["d1squbar"]]
  d2squbar <- ums2pg.res[["d2squbar"]]
  dsqubar <- ums2pg.res[["dsqubar"]]

  if(!is.null(inputvalues[["ndelta"]])) ndelta <- inputvalues[["ndelta"]]
  if(!is.null(inputvalues[["nr"]])) nr <- inputvalues[["nr"]]
  
 } else if(invtyp=="pg")
 {
 # input values PG scheme

  snames <- inputvalues[["snames"]]
  d1 <- inputvalues[["d1"]]
  d2 <- inputvalues[["d2"]]
  ndelta <- inputvalues[["ndelta"]]
  nr <- inputvalues[["nr"]]
  high <- inputvalues[["high"]]
  low <- inputvalues[["low"]]
  sigma.low <- inputvalues[["sigma.low"]]
  sigma.high <- inputvalues[["sigma.high"]]

  # sample sizes
  n1 <- length(d1)
  n2 <- length(d2)
  n <- n1 + n2

  # Compute the average of the 1st data set
  d1bar <- sum(d1)/n1 #mean
  # Compute the average of the 2nd data set
  d2bar <- sum(d2)/n2
  # Compute the average of the combined data set
  dbar <- sum(d1,d2)/n
  # Compute the mean square average of the 1st data set
  d1squbar <- sum(d1 * d1)/n1
  # Compute the mean square average of the 2nd data set
  d2squbar <- sum(d2 * d2)/n2
  # Compute the mean square average of the combined data set
  dsqubar <- (sum(d1 * d1) + sum(d2 * d2))/n
  # Neater way of computing the above using a vector dot product
  # (d1 %*% d1 + d2 %*% d2)/n
  # Compute the standard deviation of the 1st data set
  d1std <- sd(d1)
  # Compute the standard deviation of the 2nd data set
  d2std <- sd(d2)
  # Compute the standard deviation of the combined data set
  dstd <- sd(c(d1,d2))

 } else
 {
  stop(paste("\nNo valid input type for input values given.\n\n",sep=""))
 }

  
#### constants
# Set prior limits (assumed the same for each data set) on mean (low,high),
# and prior limits (assumed the same for each data set) on the
# standard deviation (sigmalow, sigmahigh).

  # range mean
  Rc <- high - low
  # range sd
  Rsigma <- sigma.high / sigma.low

# END INPUT VALUES


################################################################################
# pCSk
# Compute pCSk = p(C,S|D_1,D_2,I) * p(D_1,D_2|I)
# according to GL Bretthorst (1993) formulas and appendix C of Gregory (2005)

  z <- n * (dsqubar - dbar^2)
  uL <- sqrt(n/2) * (low - dbar)
  uH <- sqrt(n/2) * (high - dbar)

  (uH/uL) > 0

  # error function
  errf <- function(ERR) return( 2*pnorm(ERR*sqrt(2))-1 )

  # formula to calculate area below the curve / integral
  pCSk.hypo <- function(n, Rc, Rsigma, sigma1, z, uH, uL)
  {  
    fnc1 <- function(sigma1)
    {
      ( (2*pi)^(-n/2)*sqrt(pi/(2*n)) ) /
      ( 4*Rc*log(Rsigma) ) *
      sigma1^(-n) *
      exp(-z/(2*sigma1^2)) *
      ( errf(ERR=uH/sigma1) - errf(ERR=uL/sigma1) )
    }
    # integrate
    pCSk <- integrate(fnc1, lower=sigma.low, upper=sigma.high)$value
    return(pCSk)
  }
  #call:
  #pCSk <- pCSk.hypo(n=n, Rc=Rc, Rsigma=Rsigma, z=z, uH=uH, uL=uL)

# BROB
# error function
  errf.brob <- function(ERR) return( as.brob( 2*pnorm(ERR*sqrt(2))-1 ))

# formula to calculate area below the curve / integral
  pCSk.hypo.brob <- function(n,Rc,Rsigma,sigma1,z,uH,uL,Nsteps)
  {  
    fnc1.brob <- function(sigma1, Nsteps=Nsteps)
    {
      as.brob( (2*pi)^(-n/2)*sqrt(pi/(2*n)) ) /
        as.brob( 4*Rc*log(Rsigma) ) *
        as.brob( sigma1)^(-n) *
        exp(as.brob( -z/(2*sigma1^2)) ) *
        ( errf.brob(ERR=uH/sigma1) - errf.brob(ERR=uL/sigma1) )
    }
    # integrate
    pCSk.brob <- simpsonrule.brob(fx=fnc1.brob, sL=sigma.low, sH=sigma.high, Nsteps=Nsteps)
    return(pCSk.brob)
  }
  #call:
  #pCSk.brob <- pCSk.hypo.brob(n=n, Rc=Rc, Rsigma=Rsigma, z=z, uH=uH, uL=uL, Nsteps=Nsteps)
  #fnc1.brob(sigma.low)
  #fnc1.brob(sigma.high)

################################################################################
# pCSbark
# Compute pCSbark = p(C,Sbar|D_1,D_2,I) * p(D_1,D_2|I)

  u1A <- function(A) return( n1/2 * (d1squbar - 2 * A * d1bar + A^2) )
  u2A <- function(A) return( n2/2 * (d2squbar - 2 * A * d2bar + A^2) )

  # formula to calculate area below the curve / integral
  pCSbark.hypo <- function(n, n1, n2, Rc, Rsigma, sigma.high, sigma.log, low, high)
  {
    fnc2 <- function(A)
    {
      (2*pi)^(-n/2) * 
      gamma(n1/2) * 
      gamma(n2/2) /
      ( 16*Rc*log(Rsigma)^2 ) *
      u1A(A)^(-n1/2) *
      u2A(A)^(-n2/2) *
      ( pgamma(u1A(A)/sigma.high^2,n1/2) - pgamma(u1A(A)/sigma.low^2,n1/2) ) *
      ( pgamma(u2A(A)/sigma.high^2,n2/2) - pgamma(u2A(A)/sigma.low^2,n2/2) )
    }
    # integrate
    pCSbark <- integrate(fnc2, lower=low, upper=high)$value
    return(pCSbark)
  }
  #call:
  #pCSbark <- pCSbark.hypo(n=n, n1=n1, n2=n2, Rc=Rc, Rsigma=Rsigma, sigma.high=sigma.high, sigma.log=sigma.low, low=low, high=high)

# BROB
  u1A.brob <- function(A) return( as.brob( n1/2 * (d1squbar - 2 * A * d1bar + A^2) ))
  u2A.brob <- function(A) return( as.brob( n2/2 * (d2squbar - 2 * A * d2bar + A^2) ))

  # formula to calculate area below the curve / integral
  pCSbark.hypo.brob <- function(n, n1, n2, Rc, Rsigma, sigma.high, sigma.log, low, high, Nsteps)
  {
    fnc2.brob <- function(A)
    {
      as.brob( (2*pi))^(-n/2) *
      brob( lgamma(n1/2) ) *
      brob( lgamma(n2/2)) /
      as.brob( 16*Rc*log(Rsigma)^2 ) *
      as.brob( u1A(A))^(-n1/2) *
      as.brob( u2A(A))^(-n2/2) *
      ( as.brob( pgamma(u1A(A)/sigma.high^2,n1/2, log.p=FALSE) ) - as.brob( pgamma(u1A(A)/sigma.low^2,n1/2, log.p=FALSE) ) ) *
      ( as.brob( pgamma(u2A(A)/sigma.high^2,n2/2, log.p=FALSE) ) - as.brob( pgamma(u2A(A)/sigma.low^2,n2/2, log.p=FALSE) ) )
    }
    # integrate
    pCSbark.brob <- simpsonrule.brob(fx=fnc2.brob, sL=low, sH=high, Nsteps=Nsteps)
    return(pCSbark.brob)
  }
  #call:
  #pCSbark.brob <- pCSbark.hypo.brob(n=n, n1=n1, n2=n2, Rc=Rc, Rsigma=Rsigma, sigma.high=sigma.high, sigma.log=sigma.low, low=low, high=high, Nsteps=Nsteps)  
  
################################################################################
# pCbarSk
# Compute pCbarSk = p(Cbar,S|D_1,D_2,I) * p(D_1,D_2|I)

  z1 <- n1 * (d1squbar - d1bar^2)
  u1H <- sqrt(n1/2) * (high - d1bar)
  u1L <- sqrt(n1/2) * (low - d1bar)
  z2 <- n2 * (d2squbar - d2bar^2)
  u2H <- sqrt(n2/2) * (high - d2bar)
  u2L <- sqrt(n2/2) * (low - d2bar)

  (u1H/u1L) > 0
  (u2H/u2L) > 0

# formula to calculate area below the curve / integral
  pCbarSk.hypo <- function(n, Rc, Rsigma, n1, n2, z1, z2, u1H, u1L, u2H, u2L)
  {  
    fnc3 <- function(sigma1)
    {
      (2*pi)^(-n/2) * pi / ( 8 * Rc^2 * log(Rsigma) * sqrt(n1*n2) ) * sigma1^(-n+1) *
      exp(-(z1+z2)/(2*sigma1^2)) *
      (errf(ERR=u1H/sigma1) - errf(ERR=u1L/sigma1)) *
      (errf(ERR=u2H/sigma1) - errf(ERR=u2L/sigma1))
    }
    # integrate
    pCbarSk <- integrate(fnc3, lower=sigma.low, upper=sigma.high)$value
    return(pCbarSk)
  }
  #call:
  #pCbarSk <- pCbarSk.hypo(n=n, Rc=Rc, Rsigma=Rsigma, n1=n1, n2=n2, z1=z1, z2=z2, u1H=u1H, u1L=u1L, u2H=u2H, u2L=u2L)
  
# BROB
  pCbarSk.hypo.brob <- function(n, Rc, Rsigma, n1, n2, z1, z2, u1H, u1L, u2H, u2L, Nsteps)
  {  
    fnc3.brob <- function(sigma1)
    {
      as.brob(2*pi)^(-n/2) * pi /
      as.brob( ( 8 * Rc^2 * log(Rsigma) * sqrt(n1*n2) )) *
      as.brob(sigma1)^(-n+1) *
      exp(as.brob( -(z1+z2)/(2*sigma1^2)) ) *
      (errf.brob(ERR=u1H/sigma1) - errf.brob(ERR=u1L/sigma1)) *
      (errf.brob(ERR=u2H/sigma1) - errf.brob(ERR=u2L/sigma1))
    }
    # integrate
    pCbarSk.brob <- simpsonrule.brob(fx=fnc3.brob, sL=sigma.low, sH=sigma.high, Nsteps=Nsteps)
    return(pCbarSk.brob)
  }
  #call:
  #pCbarSk.brob <- pCbarSk.hypo.brob(n=n, Rc=Rc, Rsigma=Rsigma, n1=n1, n2=n2, z1=z1, z2=z2, u1H=u1H, u1L=u1L, u2H=u2H, u2L=u2L, Nsteps=Nsteps)
    
  
################################################################################
# pCbarSbark
# Compute pCbarSbark = p(Cbar,Sbar|D_1,D_2,I) * p(D_1,D_2|I)

  (u1H/u1L) > 0
  (u2H/u2L) > 0

# formula to calculate area below the curve / integral
  pCbarSbark.A.hypo <- function(n1, z1, u1H, u1L)
  {
    fnc4.A <- function(sigma1)
    {
      sigma1^(-n1) * exp(-z1/(2*sigma1^2)) *
      (errf(ERR=u1H/sigma1) - errf(ERR=u1L/sigma1))
    }
    # integrate
    pCbarSbark.A <- integrate(fnc4.A, lower=sigma.low, upper=sigma.high)$value
    return(pCbarSbark.A)
  }
  #call:
  #pCbarSbark.A <- pCbarSbark.A.hypo(n1=n1, z1=z1, u1H=u1H, u1L=u1L)

# BROB
  pCbarSbark.A.hypo.brob <- function(n1, z1, u1H, u1L, Nsteps)
  {
    fnc4.A.brob <- function(sigma1)
    {
      as.brob(sigma1) ^(-n1) *
      exp(as.brob(-z1/(2*sigma1^2)) ) *
      (errf.brob(ERR=u1H/sigma1) - errf.brob(ERR=u1L/sigma1))
    }
    # integrate
    pCbarSbark.A.brob <- simpsonrule.brob(fx=fnc4.A.brob, sL=sigma.low, sH=sigma.high, Nsteps=Nsteps)
    return(pCbarSbark.A.brob)
  }
  #call:
  #pCbarSbark.A.brob <- pCbarSbark.A.hypo.brob(n1=n1, z1=z1, u1H=u1H, u1L=u1L, Nteps=Nsteps)

  
  # formula to calculate area below the curve / integral
  pCbarSbark.B.hypo <- function(n1, z1, u2H, u2L)
  {
    fnc4.B <- function(sigma2)
    {
     sigma2^(-n2) * exp(-z2/(2*sigma2^2)) *
    (errf(ERR=u2H/sigma2) - errf(ERR=u2L/sigma2))
    }
    # integrate
    pCbarSbark.B <- integrate(fnc4.B, lower=sigma.low, upper=sigma.high)$value
    return(pCbarSbark.B)
  }
  #call:
  #pCbarSbark.B <- pCbarSbark.B.hypo(n1=n1, z1=z1, u2H=u2H, u2L=u2L)

# BROB
  pCbarSbark.B.hypo.brob <- function(n1, z1, u2H, u2L, Nsteps)
  {
    fnc4.B.brob <- function(sigma2)
    {
      as.brob( sigma2)^(-n2) *
      exp( as.brob(-z2/(2*sigma2^2)) ) *
      (errf.brob(ERR=u2H/sigma2) - errf.brob(ERR=u2L/sigma2))
    }
    # integrate
    pCbarSbark.B.brob <- simpsonrule.brob(fx=fnc4.B.brob, sL=sigma.low, sH=sigma.high, Nsteps=Nsteps)
    return(pCbarSbark.B.brob)
  }
  #call:
  #pCbarSbark.B.brob <- pCbarSbark.B.hypo.brob(n1=n1, z1=z1, u2H=u2H, u2L=u2L, Nsteps=Nsteps)

################################################################################
# compile results
  
# independent from BROB
# inputs
  # table/ dataframe with input values
  desc.df <- data.frame("Data set" = c(snames,"combined"),
                        "Sample Size" = c(n1,n2,n1+n2),
                        "Standard Deviation" = c(d1std,d2std,dstd),
                        "Variance" = c(d1std^2,d2std^2,dstd^2),
                        "Mean" = c(d1bar,d2bar,dbar),
                        check.names = FALSE)
  
  # table/ dataframe with prior values/ information
  prior.df <- data.frame("Type" = c("Prior Mean lower bound",
                                                 "Prior Mean upper bound",
                                                 "Prior Standard Deviation lower bound",
                                                 "Prior Standard Deviation upper bound",
                                                 "Number of steps for plotting p(delta | D_1, D_2, I)",
                                                 "Number of steps for plotting p(r | D_1, D_2, I)"),
                         "Value" = c(low,high,sigma.low,sigma.high,ndelta,nr),
                         check.names = FALSE)
  
  
# compute single hypotheses
if(!BROB)
{
  cat("\nCalculate pCSk")
  pCSk <- pCSk.hypo(n=n, Rc=Rc, Rsigma=Rsigma, z=z, uH=uH, uL=uL)
  cat("\nCalculate pCSbark")
  pCSbark <- pCSbark.hypo(n=n, n1=n1, n2=n2, Rc=Rc, Rsigma=Rsigma, sigma.high=sigma.high, sigma.log=sigma.low, low=low, high=high)
  cat("\nCalculate pCbarSk")
  pCbarSk <- pCbarSk.hypo(n=n, Rc=Rc, Rsigma=Rsigma, n1=n1, n2=n2, z1=z1, z2=z2, u1H=u1H, u1L=u1L, u2H=u2H, u2L=u2L)
  cat("\nCalculate pCbarSbark.A")
  pCbarSbark.A <- pCbarSbark.A.hypo(n1=n1, z1=z1, u1H=u1H, u1L=u1L)
  cat("\nCalculate pCbarSbark.B")
  pCbarSbark.B <- pCbarSbark.B.hypo(n1=n1, z1=z1, u2H=u2H, u2L=u2L)
  # compile pCbarSbark.A and pCbarSbark.B
  pCbarSbark <- ( (2*pi)^(-n/2) * pi * pCbarSbark.A * pCbarSbark.B ) / ( 8 * Rc^2 * log(Rsigma)^2 * sqrt(n1*n2) )
} else
{
# BROB
  cat("\nCalculate pCSk")
  pCSk <- pCSk.hypo.brob(n=n, Rc=Rc, Rsigma=Rsigma, z=z, uH=uH, uL=uL, Nsteps=Nsteps)
  cat("\nCalculate pCSbark")
  pCSbark <- pCSbark.hypo.brob(n=n, n1=n1, n2=n2, Rc=Rc, Rsigma=Rsigma, sigma.high=sigma.high, sigma.log=sigma.low, low=low, high=high, Nsteps=Nsteps)  
  cat("\nCalculate pCbarSk")
  pCbarSk <- pCbarSk.hypo.brob(n=n, Rc=Rc, Rsigma=Rsigma, n1=n1, n2=n2, z1=z1, z2=z2, u1H=u1H, u1L=u1L, u2H=u2H, u2L=u2L, Nsteps=Nsteps)
  cat("\nCalculate pCbarSbark.A")
  pCbarSbark.A <- pCbarSbark.A.hypo.brob(n1=n1, z1=z1, u1H=u1H, u1L=u1L, Nsteps=Nsteps)
  cat("\nCalculate pCbarSbark.B")
  pCbarSbark.B <- pCbarSbark.B.hypo.brob(n1=n1, z1=z1, u2H=u2H, u2L=u2L, Nsteps=Nsteps)
  # compile pCbarSbark.A and pCbarSbark.B
  pCbarSbark <- ( as.brob(2*pi)^(-n/2) * pi * pCbarSbark.A * pCbarSbark.B ) / ( 8 * Rc^2 * log(Rsigma)^2 * sqrt(n1*n2) )
}

  
# compute total probability
# Compute p(D_1,D_2|I) and the normalized probabilities
  cat("\nCompile results\n")
  # tot = pCSk + pCSbark + pCbarSk + pCbarSbark = p(D_1,D_2|I)^2
  tot <- pCSk + pCSbark + pCbarSk + pCbarSbark
  # pCS = probability that the means & standard deviations are the same
  pCS <- pCSk / tot
  # pCSbar = probability that the means are the same & standard deviations are different
  pCSbar <- pCSbark / tot
  # pCbarS = probability that the means are different & standard deviations are the same
  pCbarS <- pCbarSk / tot
  # pCbarSbar = probability that the means& standard deviations are different
  pCbarSbar <- pCbarSbark / tot
  # pC = probability that the means are the same independent of whether the standard deviations are the same or different
  pC <- pCS + pCSbar
  # pCbar = probability that the means are different independent of whether the standard deviations are the same or different
  pCbar <- pCbarS + pCbarSbar
  # Odds in favor of different means
  oddCbarC <- pCbar / pC
  # pS = probability that the standard deviations are the same independent of whether the means are the same or different
  pS <- pCS + pCbarS
  # pSbar = probability that the standard deviations are different independent of whether the means are the same or different
  pSbar <- pCSbar + pCbarSbar
  # Odds in favor of different standard deviations
  oddSbarS <- pSbar / pS
  # Odds in favor of different means and/or standard deviations
  odddiff <- (1 - pCS) / pCS


# create output tables
# table/ dataframe with resulting probabilities
  posterior.probs <- list(pCS,pCbarS,pCSbar,pCbarSbar,pC,pCbar,pS,pSbar,pCS,1-pCS)
  names(posterior.probs) <- c("pCS","pCbarS","pCSbar","pCbarSbar","pC","pCbar","pS","pSbar","pCS","One-minus-pCS")
  posterior.df <- data.frame("short" = c("C,S","Cbar,S","C,Sbar","Cbar,Sbar","C","Cbar","S","Sbar","C, S","Cbar, Sbar"),
                             "Hypothesis" = c(
                                           "Same Mean, same Standard Deviation",
                                           "Different Means, same Standard Deviation",
                                           "Same Mean, different Standard Deviations",
                                           "Different Means, different Standard Deviations",
                                           "The Means are the same",
                                           "The Means are different",                            
                                           "The Standard Deviations are the same",
                                           "The Standard Deviations are different",                            
                                           "Same Means and Standard Deviations",
                                           "One or Both are different"),                           
                                           "Probability" = NA,
                                           "Probability = exp(x)" = NA,
                                           "sign" = NA,
                                           check.names = FALSE)

# table/ dataframe with odds ratio results
  OR.post <- list(oddCbarC,oddSbarS,odddiff, 1/oddCbarC,1/oddSbarS,1/odddiff)
  names(OR.post) <- c("oddCbarC","oddSbarS","odddiff", "One-div-oddCbarC","One-div-oddSbarS","One-div-odddiff")
  OR.res.df <- data.frame("Hypothesis" = c("The odds Ratio in favour of a difference (means)",
                                          "The Odds Ratio in favour of a difference (standard deviations)",
                                          "The Odds Ratio in favour of a difference (means and standard deviations)",
                                          "The Odds Ratio in favour of the same (means)",
                                          "The Odds Ratio in favour of the same (standard deviations)",
                                          "The Odds Ratio in favour of the same (means and standard deviations)"),
                               "Odds Ratio" = NA,
                               "Odds Ratio = exp(x)" = NA,
                               "sign" = NA,
                               check.names = FALSE)
  
  if(!BROB)
  {
    posterior.df[,"Probability"] <- unlist(posterior.probs)
    OR.res.df["Odds Ratio"] <- unlist(OR.post)
    tot.out <- as.numeric(tot)
  } else
  {
    # posterior probs
    posterior.df[,"Probability"] <- unlist(lapply(posterior.probs,as.numeric))
    posterior.df[,"Probability = exp(x)"] <- unlist(lapply(posterior.probs, function(x) x@x))  
    posterior.df[,"sign"] <- unlist(lapply(posterior.probs, function(x) x@positive))  
    # ORs  
    OR.res.df[,"Odds Ratio"] <-  unlist(lapply(OR.post,as.numeric))
    OR.res.df[,"Odds Ratio = exp(x)"] <- unlist(lapply(OR.post, function(x) x@x))
    OR.res.df[,"sign"] <- unlist(lapply(OR.post, function(x) x@positive))
    tot.out <- paste(signif(tot@x,digits=dig+2)," = exp(x) [sign=",tot@positive,"]",sep="")
  }
  
  #results
  
 inputvalues[["low"]] <- low
 inputvalues[["high"]] <- high
 inputvalues[["ndelta"]] <- ndelta
 inputvalues[["n1"]] <- n1
 inputvalues[["n2"]] <- n2
 inputvalues[["n"]] <- n
 inputvalues[["d1bar"]] <- d1bar
 inputvalues[["d2bar"]] <- d2bar
 inputvalues[["Rc"]] <- Rc
 inputvalues[["Rsigma"]] <- Rsigma
 inputvalues[["dsqubar"]] <- dsqubar
 inputvalues[["dbar"]] <- dbar
 inputvalues[["sigma.high"]] <- sigma.high
 inputvalues[["sigma.low"]] <- sigma.low
 inputvalues[["d1squbar"]] <- d1squbar
 inputvalues[["d2squbar"]] <- d2squbar
 inputvalues[["d1std"]] <- d1std
 inputvalues[["d2std"]] <- d2std
 inputvalues[["nr"]] <- nr
  
  res <- list(invtyp = invtyp,
              inputvalues = inputvalues,
              prior.df = prior.df,
              desc.df = desc.df,
              tot = tot,
              tot.out = tot.out,
              posterior.probs = posterior.probs,
              posterior.df = posterior.df,
              OR.post = OR.post,
              OR.res.df = OR.res.df,
              BROB = BROB
              )
  
# print tables
  if(print.res)
  {
    DiM.print.pg(res)
  }

return(res)
}
# call:
# DiM.pg(invtyp="ums", inputvalues=res.SIB.NRFtotal, print.res=TRUE, BROB=TRUE)
################################################################################


################################################################################
DiM.extract.limits <- function(DiM.res, scaleL=NULL, scaleH=NULL, low=NULL, high=NULL, sigma.high=NULL, sigma.low=NULL, change=FALSE, dig=4, fac.brob=1)
{
  #prior
  if(is.null(low)) low <- DiM.res[["inputvalues"]]$low
  if(is.null(high)) high <- DiM.res[["inputvalues"]]$high
  ndelta <- DiM.res[["inputvalues"]]$ndelta
  nr <- DiM.res[["inputvalues"]]$nr
  if(is.null(sigma.low)) sigma.low <- DiM.res[["inputvalues"]]$sigma.low  
  if(is.null(sigma.high)) sigma.high <- DiM.res[["inputvalues"]]$sigma.high
  if(is.null(scaleL)) scaleL <- DiM.res[["inputvalues"]]$scaleL
  if(is.null(scaleH)) scaleH <- DiM.res[["inputvalues"]]$scaleH

  if(length(c(scaleL,scaleH, low, high)) !=  4)
  {
    cat("\nAny of ['scaleL', 'scaleH', 'low', 'high'] is NULL. Stopping.\n\n")
    stop()
  }
  
  #delta
  delta.low <- low - high
  delta.high <- high - low
  delta.delta <- (delta.high - delta.low) / (ndelta - 1)  
  delta.sek <- seq(delta.low,delta.high,delta.delta / fac.brob)
  delta.sek.l <- length(delta.sek)
  
  #ratio
  # upper and lower bounds of the plot can be tweaked with scaleL and scaleH.
  # use as input values to allow to change them:
  d1std <- DiM.res[["inputvalues"]]$d1std
  d2std <- DiM.res[["inputvalues"]]$d2std  
  r.low <- d1std / (scaleL * d2std)
  r.high <- scaleH * d1std / d2std
  r.delta <- (r.high - r.low) / (nr - 1)
  r.sek <- seq(r.low,r.high,r.delta / fac.brob)
  r.sek.l <- length(r.sek)
  
  #output
  res <- structure(list(low=low,
                 high=high,
                 ndelta=ndelta,
                 sigma.low=sigma.low,
                 sigma.high=sigma.high,
                 delta.low=delta.low,
                 delta.high=delta.high,
                 delta.delta=delta.delta,
                 delta.sek.l=delta.sek.l,
                 scaleL=scaleL,
                 scaleH=scaleH,
                 d1std=d1std,
                 d2std=d2std,
                 r.low=r.low,
                 r.high=r.high,
                 r.delta=r.delta,
                 r.sek.l=r.sek.l,
                 fac.brob=fac.brob
                 ))
  cat("\n\n")
  cat(paste(format(names(res), width = 15L, justify = "right"), 
            format(res, digits = dig), sep = " = "), sep = "\n")
  cat("\n\n")
  
  #output
  
  #change original values
  if(change)
  {
    DiM.res[["inputvalues"]]$low <- low
    DiM.res[["inputvalues"]]$high <- high
    DiM.res[["inputvalues"]]$ndelta <- ndelta
    DiM.res[["inputvalues"]]$sigma.low <- sigma.low
    DiM.res[["inputvalues"]]$sigma.high <- sigma.high
    DiM.res[["inputvalues"]]$delta.low <- delta.low
    DiM.res[["inputvalues"]]$delta.high <- delta.high
    DiM.res[["inputvalues"]]$delta.delta <- delta.delta
    DiM.res[["inputvalues"]]$scaleL <- scaleL
    DiM.res[["inputvalues"]]$scaleH <- scaleH
    DiM.res[["inputvalues"]]$d1std <- d1std
    DiM.res[["inputvalues"]]$d2std <- d2std
    DiM.res[["inputvalues"]]$r.low <- r.low
    DiM.res[["inputvalues"]]$r.high <- r.high
    DiM.res[["inputvalues"]]$r.delta <- r.delta
    DiM.res[["inputvalues"]]$fac.brob <- fac.brob
    return(DiM.res)
  }
}
# call:
# DiM.extract.limits(DiM.res, change=FALSE)
# Dim.res.newlimits <- DiM.extract.limits(DiM.res, scaleL=30, scaleH=10, change=TRUE)
# check:
# DiM.extract.limits(DiM.res.newlimits, change=FALSE)
################################################################################


################################################################################
DiM.plot.calc.pg <- function(DiM.res, BROB=FALSE, fac.brob=1, cMasses=c(0.89,0.96), low=NULL, high=NULL, sigma.low=NULL, sigma.high=NULL, scaleL=NULL, scaleH=NULL, Nsteps=100)
{
  
  cat("\n\nCalculate graphical output for the various probabilities...\n")
  library(progress)
  
  # for large numbers
  if(BROB)
  {
    library(Brobdingnag)
    fac.brob <- 0.1
    cat("\nUsing log() and package 'Brobdingnag'\nAll results are expressed as log(RESULT)\n")
  }  
  
  # error function
  errf <- function(ERR) return( 2*pnorm(ERR*sqrt(2))-1 )

# BROB
  errf.brob <- function(ERR) return( as.brob( 2*pnorm(ERR*sqrt(2))-1 ))
  
  # map values
  if(is.null(low)) low <- DiM.res[["inputvalues"]]$low
  if(is.null(high)) high <- DiM.res[["inputvalues"]]$high
  ndelta <- DiM.res[["inputvalues"]]$ndelta
  n1 <- DiM.res[["inputvalues"]]$n1
  n2 <- DiM.res[["inputvalues"]]$n2
  n <- DiM.res[["inputvalues"]]$n
  d1bar <- DiM.res[["inputvalues"]]$d1bar
  d2bar <- DiM.res[["inputvalues"]]$d2bar
  Rc <- DiM.res[["inputvalues"]]$Rc
  Rsigma <- DiM.res[["inputvalues"]]$Rsigma
  dsqubar <- DiM.res[["inputvalues"]]$dsqubar
  dbar <- DiM.res[["inputvalues"]]$dbar
  if(is.null(sigma.low)) sigma.low <- DiM.res[["inputvalues"]]$sigma.low  
  if(is.null(sigma.high)) sigma.high <- DiM.res[["inputvalues"]]$sigma.high
  d1squbar <- DiM.res[["inputvalues"]]$d1squbar
  d2squbar <- DiM.res[["inputvalues"]]$d2squbar
  d1std <- DiM.res[["inputvalues"]]$d1std
  d2std <- DiM.res[["inputvalues"]]$d2std
  nr <- DiM.res[["inputvalues"]]$nr

  if(is.null(scaleL)) scaleL <- DiM.res[["inputvalues"]]$scaleL
  if(is.null(scaleH)) scaleH <- DiM.res[["inputvalues"]]$scaleH
  if(is.null(fac.brob)) fac.brob <- DiM.res[["inputvalues"]]$fac.brob
  
  # range mean: Rc <- high - low
  # range sd: Rsigma <- sigma.high / sigma.low
  
################################################################################
# p(delta|S,D_1,D_2,I) problem
# plot The difference in the means (delta) if sd is the same

# create sequence for values to calculate the posterior probability density(ies)
  delta.low <- low - high
  delta.high <- high - low
  delta.delta <- (delta.high - delta.low) / (ndelta - 1)
  
  beta.low <- 2 * low
  beta.high <- 2 * high
  delta.N <- (n1 - n2) / n
  b <- (n1 * d1bar - n2 * d2bar) / (2 * n)

  # functions
  V <- function(delta, BETA)
  {
   (n/2) * (dsqubar - (2 * delta * b) - (BETA * dbar) + (BETA^2/4) + (delta^2/4) + (delta * BETA * delta.N/ 2) )
  }

  fndbV <- function(BETA)
  {
    gamma(n/2) / (8 * Rc^2 * log(Rsigma)) * V(delta=delta,BETA=BETA)^(-n/2) *
    ( pgamma(V(delta=delta,BETA=BETA)/sigma.high^2,n/2) - pgamma(V(delta=delta,BETA=BETA)/sigma.low^2,n/2))
  }

# BROB
  V.brob <- function(delta, BETA)
  {
    as.brob( (n/2) * (dsqubar - (2 * delta * b) - (BETA * dbar) + (BETA^2/4) + (delta^2/4) + (delta * BETA * delta.N/ 2) ) )
  }

  fndbV.brob <- function(BETA)
  {
    brob( lgamma(n/2) ) /
    as.brob( (8 * Rc^2 * log(Rsigma)) ) * 
    as.brob( V(delta=delta,BETA=BETA))^(-n/2) *
    (
      brob( pgamma(V(delta=delta,BETA=BETA) / sigma.high^2,n/2,log.p=TRUE,lower=FALSE)) -
      brob( pgamma(V(delta=delta,BETA=BETA) / sigma.low^2,n/2,log.p=TRUE,lower=FALSE)) 
    )  
  }

  # integrate for each element of delta.sek
  # fac.brob ~ 10% of non-brob (fac = 1) -> calculation time!
  cat("Calculate pdel.SD1D2I and pdelA...\n")
  delta.sek <- seq(delta.low,delta.high,delta.delta / fac.brob) 
  delta.sek.l <- length(delta.sek)
  pdel.SD1D2I <- vector(mode="list",length=delta.sek.l)
  pdel.SD1D2I.sum <- 0
  pb <- progress_bar$new(total = delta.sek.l)
  for(i in 1:delta.sek.l)
  {
    #print(i)
    pb$tick()
    delta <- delta.sek[i]
    if(!BROB)
    {
      pdel.SD1D2I[[i]] <- integrate(fndbV, lower=beta.low, upper=beta.high)$value
    } else
    {
      pdel.SD1D2I[[i]] <- simpsonrule.brob(fx=fndbV.brob, sL=beta.low, sH=beta.high, Nsteps=Nsteps)
    }
    pdel.SD1D2I.sum <- pdel.SD1D2I.sum + pdel.SD1D2I[[i]]
  }
  
  pdel.SD1D2I[[50]]
  pdel.SD1D2I.sum
  
  # compile
  pdelA <- lapply(pdel.SD1D2I, function(x) x/ (delta.delta * pdel.SD1D2I.sum))
  #pdelA <- pdel.SD1D2I / (delta.delta * sum(pdel.SD1D2I))
  # create dataframe with input and results
  if(!BROB)
  {
    pdelta.df <- data.frame(delta.sek = delta.sek,
                            pdel.SD1D2I = unlist(pdel.SD1D2I),
                            pdel.SD1D2I.sign = NA,
                            pdelA = unlist(pdelA),
                            pdelA.sign = NA
                            )
  } else
  {
    pdelta.df <- data.frame(delta.sek = delta.sek,
                            pdel.SD1D2I = unlist(lapply(pdel.SD1D2I,function(x) x@x)),
                            pdel.SD1D2I.sign = unlist(lapply(pdel.SD1D2I,function(x) x@positive)),
                            pdelA = unlist(lapply(pdelA,function(x) x@x)),
                            pdelA.sign = unlist(lapply(pdelA,function(x) x@positive))
                            )
  }    
  head(pdelta.df)
 
################################################################################
# p(delta|Sbar,D_1,D_2,I) problem
# plot the difference in means assuming sds are different.

# TODO tweak range for each plot!
  #redundant, not necessary
  #delta.low <- low - high
  #delta.high <- high - low
  #delta.delta <- (delta.high - delta.low) / (ndelta - 1)
  #beta.low <- 2 * low
  #beta.high <- 2 * high

  w1 <- function(delta, BETA)
  {
    (n1/2) * (d1squbar - (BETA + delta) * d1bar + ((BETA + delta)^2/4) )
  }

  w2 <- function(delta, BETA)
  {
    (n2/2) * (d2squbar - (BETA - delta) * d2bar + ((BETA - delta)^2/4) )
  } 

  fndbVbar <- function(BETA)
  {
    gamma(n1/2) * gamma(n2/2) / (16 * Rc^2 * log(Rsigma)^2) *
    w1(delta=delta,BETA=BETA)^(-n1/2) * w2(delta=delta,BETA=BETA)^(-n2/2) *
    ( pgamma(w1(delta=delta,BETA=BETA)/sigma.high^2,n1/2) - pgamma(w1(delta=delta,BETA=BETA)/sigma.low^2,n1/2) ) *
    ( pgamma(w2(delta=delta,BETA=BETA)/sigma.high^2,n2/2) - pgamma(w2(delta=delta,BETA=BETA)/sigma.low^2,n2/2) )
  }

# BROB
  w1.brob <- function(delta, BETA)
  {
    as.brob( (n1/2) * (d1squbar - (BETA + delta) * d1bar + ((BETA + delta)^2/4) ) )
  }

  w2.brob <- function(delta, BETA)
  {
    as.brob( (n2/2) * (d2squbar - (BETA - delta) * d2bar + ((BETA - delta)^2/4) ) )
  } 

  fndbVbar.brob <- function(BETA)
  {
    brob( lgamma(n1/2) ) * brob( lgamma(n2/2) ) / as.brob( (16 * Rc^2 * log(Rsigma)^2) ) *
    as.brob( w1(delta=delta,BETA=BETA)^(-n1/2) ) *
    as.brob( w2(delta=delta,BETA=BETA)^(-n2/2) ) *
    ( brob( pgamma(w1(delta=delta,BETA=BETA)/sigma.high^2,n1/2, log.p=TRUE,lower=FALSE) ) -
      brob( pgamma(w1(delta=delta,BETA=BETA)/sigma.low^2,n1/2, log.p=TRUE,lower=FALSE) ) ) *
    ( brob( pgamma(w2(delta=delta,BETA=BETA)/sigma.high^2,n2/2, log.p=TRUE,lower=FALSE) ) -
      brob( pgamma(w2(delta=delta,BETA=BETA)/sigma.low^2,n2/2, log.p=TRUE,lower=FALSE) ) )
  }

# integrate for each element of delta.sek
# TODO: do it with sapply
  cat("Calculate pdel.SbarD1D2I and pdelB...\n")
  pdel.SbarD1D2I <- vector(mode="list",length=delta.sek.l)
  pdel.SbarD1D2I.sum <- 0
  pb <- progress_bar$new(total = delta.sek.l)
  for(i in 1:delta.sek.l)
  {
    #print(i)
    pb$tick()
    delta <- delta.sek[i]
    if(!BROB)
    {
      pdel.SbarD1D2I[[i]] <- integrate(fndbVbar, lower=beta.low, upper=beta.high)$value
    } else
    {
      pdel.SbarD1D2I[[i]] <- simpsonrule.brob(fx=fndbVbar.brob, sL=beta.low, sH=beta.high, Nsteps=Nsteps)
    }
    pdel.SbarD1D2I.sum <- pdel.SbarD1D2I.sum + pdel.SbarD1D2I[[i]]
  }

  pdel.SbarD1D2I[[50]]
  pdel.SbarD1D2I.sum
  #compile
  #pdelB <- pdel.SbarD1D2I / (delta.delta * sum(pdel.SbarD1D2I))
  pdelB <- lapply(pdel.SbarD1D2I, function(x) x/ (delta.delta * pdel.SbarD1D2I.sum))
  # add to dataframe
  if(!BROB)
  {
    pdelta.df[,c("pdel.SbarD1D2I","pdel.SbarD1D2I.sign","pdelB","pdelB.sign")] <- data.frame(unlist(pdel.SbarD1D2I), NA, unlist(pdelB),NA)
  } else
  {
    pdelta.df[,c("pdel.SbarD1D2I","pdel.SbarD1D2I.sign","pdelB","pdelB.sign")] <- data.frame(
                            unlist(lapply(pdel.SbarD1D2I,function(x) x@x)),
                            unlist(lapply(pdel.SbarD1D2I,function(x) x@positive)),
                            unlist(lapply(pdelB,function(x) x@x)),
                            unlist(lapply(pdelB,function(x) x@positive))
    )
  }    
  head(pdelta.df)

################################################################################
# p(delta|D_1,D_2,I) problem
# plot the difference in means independent of the sds are the same or not
# This is a weighted average of p(delta|S,D_1,D_2,I) and p(delta|Sbar,D_1,D_2,I)
# with weights are pS and pSbar ie. probs of sds are same/ different.
# It also calcs the CIs.

  # to remember - s.a.:
  # p(delta|S,D_1,D_2,I)
  # xdel <- delta.sek                         #vector
  # xda <- pdelta.df[,c("delta.sek","pdelA")] #table/ data.frame with 2 cols
  #
  # p(delta|Sbar,D_1,D_2,I)
  # xdel <- delta.sek                         #vector
  # xdB <- pdelta.df[,c("delta.sek","pdelB")] #table/ data.frame with 2 cols

  # taken from results
  pS <- DiM.res[["posterior.probs"]]$pS
  pSbar <- DiM.res[["posterior.probs"]]$pSbar

  #pdelave = average of pdelA and pdelB
  cat("Calculate average of pdelA and pdelB...\n")
  if(!BROB)
  {
    pdelave <- (pS * pdelta.df[,"pdelA"] + pSbar * pdelta.df[,"pdelB"])
    pdelta.df[,c("pdelave","pdelave.sign","pdelave.delta.delta","pdelave.delta.delta.sign")] <- data.frame(as.numeric(pdelave / sum(pdelave)),NA,as.numeric(pdelave / sum(pdelave))/delta.delta,NA)
  } else
  {
    pS.pdelA <- lapply(pdelA, function(x) x*pS)
    pSbar.pdelB <- lapply(pdelB, function(x) x*pSbar)
    pS.pdelA.pSbar.pdelB <- sapply(seq_along(1:length(pS.pdelA)), function(x) pS.pdelA[[x]] + pSbar.pdelB[[x]] )
    pS.pdelA.pSbar.pdelB.sum <- 0
    for(i in 1:length(pS.pdelA.pSbar.pdelB)) pS.pdelA.pSbar.pdelB.sum <- pS.pdelA.pSbar.pdelB[[i]] + pS.pdelA.pSbar.pdelB.sum
    pS.pdelA.pSbar.pdelB <- lapply(pS.pdelA.pSbar.pdelB, function(x) x/pS.pdelA.pSbar.pdelB.sum)
    
    pdelta.df[,"pdelave"] <- unlist(lapply(pS.pdelA.pSbar.pdelB, function(x) x@x))
    pdelta.df[,"pdelave.sign"] <- unlist(lapply(pS.pdelA.pSbar.pdelB, function(x) x@positive))
    pdelta.df[,"pdelave.delta.delta"] <- unlist(lapply(pS.pdelA.pSbar.pdelB, function(x) (x/delta.delta)@x ))
    pdelta.df[,"pdelave.delta.delta.sign"] <- unlist(lapply(pS.pdelA.pSbar.pdelB, function(x) (x/delta.delta)@positive ))
  } 
  head(pdelta.df)
  
######## FUNCTION
  # function to calculate and extract HDI infos from densities
  extractcMass.values <- function(dframe, cMass=0.95, BROB=FALSE)
  {
    if(BROB) return(NA)
    library(HDInterval)
    orig.names <- colnames(dframe)
    dframe.d <- dim(dframe)
    stopifnot(dframe.d[2] == 2)
    colnames(dframe) <- c("x","y")
    dframe.cMass <- HDInterval:::hdi.density(object=dframe, credMass=cMass)
    res <- as.list(dframe.cMass)
    res[["cMass"]] <- cMass
    res[["height"]] <- attr(dframe.cMass,"height")
    dframe[,"cMass.incl"] <- dframe[,"y"] >= res[["height"]]
    colnames(dframe)[1:2] <- orig.names
    res[["dframe"]] <- dframe
    #x values for cMass
    #dframe[,"delta.sek"][dframe[,"cMass.incl"] == TRUE]
    return(res)
  } 
  #call:
  #extractcMass.values(dframe, cMass=0.95)

  #create 89% and 96% HDI infos    
  cat("Calculate HDIs for the average of pdelA and pdelB...\n")
  pdelta.df.dframe <- pdelta.df[,c("delta.sek","pdelave.delta.delta")]
  pdelave.delta.delta.cMass1 <- extractcMass.values(pdelta.df.dframe, cMass=cMasses[1], BROB=BROB)
  pdelave.delta.delta.cMass2 <- extractcMass.values(pdelta.df.dframe, cMass=cMasses[2], BROB=BROB)
  
  #plot al three lines next to each other
  #plot(pdelta.df[,"delta.sek"],pdelta.df[,"pdelave"]/delta.delta, type="l", col="darkred")
  #abline(h=attr(hdiobj.cMass1,"height"))
  #lines(pdelta.df[,"delta.sek"],pdelta.df[,"pdelA"], col="blue")
  #lines(pdelta.df[,"delta.sek"],pdelta.df[,"pdelB"], col="green")

  # Plot p(delta|v,D_1,D_2,I), p(delta|vbar,D_1,D_2,I) and the weighted average of them
  # plus CIs.

################################################################################
# p(r|C,D_1,D_2,I) problem
# Plot the ratio of sds if the means are the same

  #use as input values to allow to change them:
  #scaleL <- 30
  #scaleH <- 4

  r.low <- d1std / (scaleL * d2std)
  r.high <- scaleH * d1std / d2std
  r.delta <- (r.high - r.low) / (nr - 1)
  AL <- low
  AH <- high

  u <- function(r) ( (n1 * d1squbar) / r^2 ) + n2 * d2squbar
  v <- function(r) ( (n1 * d1bar) / r^2 ) + n2 * d2bar
  w <- function(r) n1 / r^2 + n2
  z <- function(r) u(r=r) - v(r=r)^2 / w(r=r)
  XH <- function(sigma) sqrt(w(r=r) / (2 * sigma^2)) * (AH - v(r=r)/w(r=r))
  XL <- function(sigma) sqrt(w(r=r) / (2 * sigma^2)) * (AL - v(r=r)/w(r=r))

  fnrS <- function(sigma)
  {
    (2*pi)^(-n/2) * sqrt(pi/ (8*w(r=r))) / (Rc * log(Rsigma)^2) *
    r^(-(n1+1)) * sigma^(-n) *
    exp(-z(r=r)/ (2*sigma^2)) * ( errf(ERR=XH(sigma=sigma)) - errf(ERR=XL(sigma=sigma)) )
  }

  # BROB
  u.brob <- function(r) ( (n1 * d1squbar) / r^2 ) + n2 * d2squbar
  v.brob <- function(r) ( (n1 * d1bar) / r^2 ) + n2 * d2bar
  w.brob <- function(r) n1 / r^2 + n2
  z.brob <- function(r) u(r=r) - v(r=r)^2 / w(r=r)
  XH.brob <- function(sigma) sqrt(w(r=r) / (2 * sigma^2)) * (AH - v(r=r)/w(r=r))
  XL.brob <- function(sigma) sqrt(w(r=r) / (2 * sigma^2)) * (AL - v(r=r)/w(r=r))
  
  fnrS.brob <- function(sigma)
  {
    as.brob(2*pi)^(-n/2) *
    as.brob( sqrt(pi/ (8*w(r=r))) / (Rc * log(Rsigma)^2) ) *
    as.brob(r)^(-(n1+1)) *
    as.brob(sigma)^(-n) *
    exp(-z(r=r)/ (2*sigma^2)) *
    ( errf.brob(ERR=XH(sigma=sigma)) - errf.brob(ERR=XL(sigma=sigma)) )
  }  
  
  # integrate for each element of r.sek
  r.sek <- seq(r.low,r.high,r.delta / fac.brob)
  r.sek.l <- length(r.sek)
  
  # do it later with sapply
  cat("Calculate pr.CD1D2I and prA...\n")
  pr.CD1D2I <- vector(mode="list",length=r.sek.l)
  pr.CD1D2I.sum <- 0
  pb <- progress_bar$new(total = delta.sek.l)
  for(i in 1:r.sek.l)
  {
    #print(i)
    pb$tick()
    r <- r.sek[i]
    if(!BROB)
    {
      pr.CD1D2I[[i]] <- integrate(fnrS, lower=sigma.low, upper=sigma.high)$value
    } else
    {
      pr.CD1D2I[[i]] <- simpsonrule.brob(fx=fnrS.brob, sL=sigma.low, sH=sigma.high, Nsteps=Nsteps)
    }
    pr.CD1D2I.sum <- pr.CD1D2I.sum + pr.CD1D2I[[i]]
  }
  
  pr.CD1D2I[[50]]
  pr.CD1D2I.sum
  #note by LG:
  #first value different to MA:
  #MA: 3.15108e-230
  #R: 3.206374e-230
  #and consequently prB is different, too:
  #MA: 1.0686e-210
  #R: 1.087353e-210
  #all other values are 100% identical...
  #assumed is that this is due to the very large numbers used and therefor caused by "accuracy"
  
  #compile
  prA <- lapply(pr.CD1D2I, function(x) x/ (r.delta * pr.CD1D2I.sum))
  #prA
  # create dataframe with input and results
  if(!BROB)
  {
    pr.df <- data.frame(r.sek, pr.CD1D2I = unlist(pr.CD1D2I), pr.CD1D2I.sign=NA, prA=unlist(prA), prA.sign=NA)
  } else
  {
    pr.df <- data.frame(r.sek = r.sek,
                        pr.CD1D2I = unlist(lapply(pr.CD1D2I,function(x) x@x)),
                        pr.CD1D2I.sign = unlist(lapply(pr.CD1D2I,function(x) x@positive)),
                        prA = unlist(lapply(prA,function(x) x@x)),
                        prA.sign = unlist(lapply(prA,function(x) x@positive))    
    )
  }  
  head(pr.df)
  
################################################################################
# p(r|Cbar,D_1,D_2,I) problem  
# This section computes and plots the posterior probability density
# for the ratio of the standard deviations assuming the means are the different. 

  #redundant, not necessary
  #scaleL <- 30
  #scaleH <- 4
  #r.low <- d1std / (scaleL * d2std)
  #r.high <- scaleH * d1std / d2std
  #r.delta <- (r.high - r.low) / (nr - 1)
  #AL <- low
  #AH <- high

  z1 <- n1 * (d1squbar - d1bar^2)
  z2 <- n2 * (d2squbar - d2bar^2)

  X1H <- function(r, sigma) sqrt(n1 / (2 * r^2 * sigma^2)) * (AH - d1bar)
  X1L <- function(r, sigma) sqrt(n1 / (2 * r^2 * sigma^2)) * (AL - d1bar)
  X2H <- function(sigma) sqrt(n2 / (2 * sigma^2)) * (AH - d2bar)
  X2L <- function(sigma) sqrt(n2 / (2 * sigma^2)) * (AL - d2bar)

  fnrSb <- function(sigma)
  {
    ( (2*pi)^(-n/2) * pi ) / ( 4 * Rc^2 * log(Rsigma)^2 * sqrt(n1*n2) ) *
    r^(-n1) * sigma^(1-n) *
    exp( -(z1 / (2 * r^2 * sigma^2)) - (z2 / (2*sigma^2)) ) *
    ( errf(ERR=X1H(r=r, sigma=sigma)) - errf(ERR=X1L(r=r, sigma=sigma)) ) *
    ( errf(ERR=X2H(sigma=sigma)) - errf(ERR=X2L(sigma=sigma)) )
  }

  # BROB
  X1H.brob <- function(r, sigma) sqrt(n1 / (2 * r^2 * sigma^2)) * (AH - d1bar)
  X1L.brob <- function(r, sigma) sqrt(n1 / (2 * r^2 * sigma^2)) * (AL - d1bar)
  X2H.brob <- function(sigma) sqrt(n2 / (2 * sigma^2)) * (AH - d2bar)
  X2L.brob <- function(sigma) sqrt(n2 / (2 * sigma^2)) * (AL - d2bar)
  
  fnrSb.brob <- function(sigma)
  {
    ( as.brob(2*pi)^(-n/2) * pi ) /
    ( 4 * Rc^2 * log(Rsigma)^2 * sqrt(n1*n2) ) *
    as.brob(r)^(-n1) *
    as.brob(sigma)^(1-n) *
    exp( as.brob( -(z1 / (2 * r^2 * sigma^2)) - (z2 / (2*sigma^2)) ) ) *
    ( errf.brob(ERR=X1H(r=r, sigma=sigma)) - errf.brob(ERR=X1L(r=r, sigma=sigma)) ) *
    ( errf.brob(ERR=X2H(sigma=sigma)) - errf.brob(ERR=X2L(sigma=sigma)) )
  }  
  
  # integrate for each element of r.sek
  cat("Calculate pr.CbarD1D2I and prB...\n")
  pr.CbarD1D2I <- vector(mode="list",length=r.sek.l)
  pr.CbarD1D2I.sum <- 0
  pb <- progress_bar$new(total = delta.sek.l)
  for(i in 1:r.sek.l)
  {
    #print(i)
    pb$tick()
    r <- r.sek[i]
    if(!BROB)
    {
      pr.CbarD1D2I[[i]] <- integrate(fnrSb, lower=sigma.low, upper=sigma.high)$value
    } else
    {
      pr.CbarD1D2I[[i]] <- simpsonrule.brob(fx=fnrSb.brob, sL=sigma.low, sH=sigma.high, Nsteps=Nsteps)
    }
    pr.CbarD1D2I.sum <- pr.CbarD1D2I.sum + pr.CbarD1D2I[[i]]
  }

  pr.CbarD1D2I[[50]]
  pr.CbarD1D2I.sum  

  #note by LG:
  #first value different to MA:
  #MA: 4.06446e-230
  #R: 4.135770e-230
  #and consequently prB is different, too:
  #MA: 4.31101e-211
  #R: 4.386639e-211
  #all other values are 100% identical...
  #assumed is that this is due to the very large numbers used and therefor caused by "accuracy"
  
  #compile
  prB <- lapply(pr.CbarD1D2I, function(x) x/ (r.delta * pr.CbarD1D2I.sum))
  #prB
  # create dataframe with input and results
  if(!BROB)
  {
    pr.df[,c("pr.CbarD1D2I","pr.CbarD1D2I.sign","prB","prB.sign")] <- data.frame(unlist(pr.CbarD1D2I), NA, prB = unlist(prB), NA)
  } else
  {
    pr.df[,c("pr.CbarD1D2I","pr.CbarD1D2I.sign","prB","prB.sign")] <- data.frame(
                        unlist(lapply(pr.CbarD1D2I,function(x) x@x)),
                        unlist(lapply(pr.CbarD1D2I,function(x) x@positive)),
                        unlist(lapply(prB,function(x) x@x)),
                        unlist(lapply(prB,function(x) x@positive))
    )
  }    
  head(pr.df)

################################################################################
# p(r|D_1,D_2,I)
# This section computes and plots the prob of the ratio of the sds independent
# of same means or not.
# This is a weighted average p(r|C,D1,D2,I) and p(r|Cbar,D1,D2,I)
# with weights pC and pCbar,

  #note by LG:
  #first (smallest values) are slightly different due to following errors from
  #above -> prA and prB, not pC and pCbar -> probably due to precision on the
  #not-log scale (ie. not using logs to integrate)

  pC <- DiM.res[["posterior.probs"]]$pC
  pCbar <- DiM.res[["posterior.probs"]]$pCbar

  #pdelave = average of pdelA and pdelB
  cat("Calculate average of prA and prB...\n")
  if(!BROB)
  {
    prave <- (pC * pr.df[,"prA"] + pCbar * pr.df[,"prB"])
    pr.df[,c("prave","prave.sign","prave.r.delta","prave.r.sign")] <- data.frame( as.numeric(prave / sum(prave)), NA, as.numeric(prave / sum(prave))/r.delta, NA)
  } else
  {
    pC.prA <- lapply(prA, function(x) x*pC)
    pCbar.prB <- lapply(prB, function(x) x*pCbar)
    pC.prA.pCbar.prB <- sapply(seq_along(1:length(pC.prA)), function(x) pC.prA[[x]] + pCbar.prB[[x]] )
    pC.prA.pCbar.prB.sum <- 0
    for(i in 1:length(pC.prA.pCbar.prB)) pC.prA.pCbar.prB.sum <- pC.prA.pCbar.prB[[i]] + pC.prA.pCbar.prB.sum
    pC.prA.pCbar.prB <- lapply(pC.prA.pCbar.prB, function(x) x/pC.prA.pCbar.prB.sum)
    
    pr.df[,c("prave","prave.sign","prave.r.delta","prave.r.sign")] <- data.frame( unlist(lapply(pC.prA.pCbar.prB, function(x) x@x)),
                                                                                  unlist(lapply(pC.prA.pCbar.prB, function(x) x@positive)),
                                                                                  unlist(lapply(pC.prA.pCbar.prB, function(x) (x/r.delta)@x )),
                                                                                  unlist(lapply(pC.prA.pCbar.prB, function(x) (x/r.delta)@positive ))
    )
  }  
  #head(pr.df)
  
  #create 89% and 96% HDI infos  
  cat("Calculate HDIs for the average of prA and prB...\n")
  prave.df.dframe <- pr.df[,c("r.sek","prave.r.delta")]
  prave.r.delta.cMass1 <- extractcMass.values(prave.df.dframe, cMass=cMasses[1], BROB=BROB)
  prave.r.delta.cMass2 <- extractcMass.values(prave.df.dframe, cMass=cMasses[2], BROB=BROB)
  
  # Here we plot  p(r|s,D1,D2,I), p(r|sbar,D1,D2,I) and their weighted average
  head(pdelta.df)
  head(pr.df)
  
  cat("Compile results...\n\n")
  DiM.plotvalues.res <- list(pdelta.df=pdelta.df,
                             pr.df=pr.df, 
                             cMasses=cMasses,
                             cMass=list(pdelave.delta.delta.cMass1=pdelave.delta.delta.cMass1,
                                        pdelave.delta.delta.cMass2=pdelave.delta.delta.cMass2,
                                        prave.r.delta.cMass1=prave.r.delta.cMass1,
                                        prave.r.delta.cMass2=prave.r.delta.cMass2)
                             )
  #str(DiM.plotvalues.res)
return(DiM.plotvalues.res)
}
#call:
#DiM.plot.calc.pg(DIM.pg.res, BROB=FALSE)
################################################################################


################################################################################
# PGplot
DiM.plot.pg <- function(DiM.plotvalues.res, cols=NULL, fac=1.05, lwd=2, alphav=0.2,
                        filling=FALSE, by1=TRUE, colb="red", dig=2, BROB=FALSE)
{

  pdelta.df <- DiM.plotvalues.res[["pdelta.df"]]
  pr.df <- DiM.plotvalues.res[["pr.df"]]
  cMasses <- DiM.plotvalues.res[["cMasses"]]
  cMass <- DiM.plotvalues.res[["cMass"]]
  
  # define colors if not done
  if(is.null(cols)) cols <- c("red","darkgreen","magenta","orange","lightgreen","violet")
  cols.rgb <- col2rgb(cols[4:6])/255
  
  # dirty hack for graphical output
  if(BROB)
  {
    pr.df.minusInf <- list(names(pr.df))
    head(pr.df)
    for(i in 1:dim(pr.df)[2])
    {
      ID <- which(pr.df[,i] == -Inf)
      MIN <- min(pr.df[pr.df[,i] > -Inf,i])
      pr.df[ID,i] <- MIN
      pr.df.minusInf[[i]] <- ID
    }  
    head(pr.df)
    pr.df.minusInf
    cat("\nIn table 'pr.df' (ratios) are values of '-Inf' - those will be replaced for the graphical output by the minimum value of those vectors\n")
    cat("Remember: exp(-Inf) = 0\n\n")
  }
  #if(BROB) pr.df[INFid] <- -Inf
  head(pr.df)
  
################################################################################
# p(delta|S,D_1,D_2,I) problem
# plot p(delta|S,D_1,D_2,I)
  if(!BROB)
  {
    y.axistextA <- expression(paste("p(",delta," | S, D" ["1"] ,", D" ["2"] ,", I)",sep=""))
  } else y.axistextA <- expression(paste("p(",delta," | S, D" ["1"] ,", D" ["2"] ,", I) on log() scale",sep=""))
  x.limsA <- range(pdelta.df[,"delta.sek"])
  y.limsA <- range(pdelta.df[,"pdelA"])*c(1,fac)
  if(by1) par(ask=TRUE)
  par(mar=c(5,6,5,5))
  par("cex.axis" = 0.8)
  plot(x.limsA[2], y.limsA[2],
      col="white",
      main="",
      xlab="",
      ylab=y.axistextA,
      bty="l",
      pre.plot=grid(),
      xlim=x.limsA, ylim=y.limsA)
  if(filling) polygon( c(pdelta.df[1,"delta.sek"],pdelta.df[,"delta.sek"], pdelta.df[dim(pdelta.df)[1],"delta.sek"],pdelta.df[1,"delta.sek"]),
                       c(y.limsA[1],pdelta.df[,"pdelA"],y.limsA[1],y.limsA[1]),
                       col=rgb(cols.rgb[1,1],cols.rgb[2,1],cols.rgb[3,1],0.5), border=NA)
  points(pdelta.df[,"delta.sek"], pdelta.df[,"pdelA"], col=cols[1], type="l", lty=1, lwd=lwd)
  mtext(expression(paste(delta," Mean Difference",sep="")), 3, line=2, cex=1.5)
  mtext("assuming the standard deviations are the same", 3, line=0.8)
  mtext(expression(paste(delta," = ", mu["1"]," - ",mu["2"])), 1, line=3, cex=1.5, col=colb)

################################################################################
# p(delta|Sbar,D_1,D_2,I) problem
# plot p(delta|Sbar,D_1,D_2,I)
  if(!BROB)
  {
    y.axistextB <- expression(paste("p(",delta," | ",bar(S),", D" ["1"] ,", D" ["2"] ,", I)",sep=""))
  } else y.axistextB <- expression(paste("p(",delta," | ",bar(S),", D" ["1"] ,", D" ["2"] ,", I) on log() scale",sep=""))
  x.limsB <- range(pdelta.df[,"delta.sek"])
  y.limsB <- range(pdelta.df[,"pdelB"])*c(1,fac)
  if(by1) par(ask=TRUE)
  par(mar=c(5,6,5,5))
  par("cex.axis" = 0.8)
  plot(x.limsB[2], y.limsB[2],
      col="white",
      main="",
      xlab="",
      ylab=y.axistextB,
      pre.plot=grid(),
      bty="l",
      xlim=x.limsB, ylim=y.limsB)
  if(filling) polygon( c(pdelta.df[1,"delta.sek"],pdelta.df[,"delta.sek"], pdelta.df[dim(pdelta.df)[1],"delta.sek"],pdelta.df[1,"delta.sek"]),
                       c(y.limsB[1],pdelta.df[,"pdelB"],y.limsB[1],y.limsB[1]),
                       col=rgb(cols.rgb[1,2],cols.rgb[2,2],cols.rgb[3,2],0.5), border=NA)
  points(pdelta.df[,"delta.sek"], pdelta.df[,"pdelB"], col=cols[2], type="l", lty=1, lwd=lwd)
  mtext(expression(paste(delta," Mean Difference",sep="")), 3, line=2, cex=1.5)
  mtext("assuming the standard deviations are different", 3, line=0.8)
  mtext(expression(paste(delta," = ", mu["1"]," - ",mu["2"])), 1, line=3, cex=1.5, , col=colb)

################################################################################
# p(delta|D_1,D_2,I) problem
# plot p(delta|v,D_1,D_2,I), p(delta|vbar,D_1,D_2,I)
  if(!BROB)
  {
    y.axistextAB <- c("Probability Density")
  } else y.axistextAB <- c("Probability Density on log() scale")
  x.limsAB <- range(pdelta.df[,"delta.sek"])
  y.limsAB <- range(pdelta.df[,c("pdelA","pdelB","pdelave.delta.delta")])*c(1,fac)
  if(by1) par(ask=TRUE)
  par(mar=c(5,6,5,5))
  par(oma=c(2,1,1,1))
  par("cex.axis" = 0.8)
  plot(x.limsAB[2], y.limsAB[2], col="white",
       main="",
       xlab="",
       ylab=y.axistextAB,
       pre.plot=grid(),
      bty="l", xlim=x.limsAB, ylim=y.limsAB)
  if(filling)
  {
    polygon( c(pdelta.df[1,"delta.sek"],pdelta.df[,"delta.sek"], pdelta.df[dim(pdelta.df)[1],"delta.sek"],pdelta.df[1,"delta.sek"]),
             c(y.limsAB[1],pdelta.df[,"pdelA"],y.limsAB[1],y.limsAB[1]),
             col=rgb(cols.rgb[1,1],cols.rgb[2,1],cols.rgb[3,1],0.5), border=NA)
    polygon( c(pdelta.df[1,"delta.sek"],pdelta.df[,"delta.sek"], pdelta.df[dim(pdelta.df)[1],"delta.sek"],pdelta.df[1,"delta.sek"]),
             c(y.limsAB[1],pdelta.df[,"pdelB"],y.limsAB[1],y.limsAB[1]),
             col=rgb(cols.rgb[1,2],cols.rgb[2,2],cols.rgb[3,2],0.5), border=NA)
    polygon( c(pdelta.df[1,"delta.sek"],pdelta.df[,"delta.sek"], pdelta.df[dim(pdelta.df)[1],"delta.sek"],pdelta.df[1,"delta.sek"]),
             c(y.limsAB[1],pdelta.df[,"pdelave.delta.delta"],y.limsAB[1],y.limsAB[1]),
             col=rgb(cols.rgb[1,3],cols.rgb[2,3],cols.rgb[3,3],0.5), border=NA)
  }    
  points(pdelta.df[,"delta.sek"], pdelta.df[,"pdelA"], col=cols[1], type="l", lty=1, lwd=lwd)
  points(pdelta.df[,"delta.sek"], pdelta.df[,"pdelB"], col=cols[2], type="l", lty=2, lwd=lwd)
  points(pdelta.df[,"delta.sek"], pdelta.df[,"pdelave.delta.delta"], col=cols[3], type="l", lty=3, lwd=lwd)

  mtext(expression(paste(delta," Mean Difference",sep="")), 3, line=3, cex=1.5)
  mtext("independent of the standard deviations being different/ the same", 3, line=1.5)
  mtext(expression(paste(delta," = ", mu["1"]," - ",mu["2"])), 1, line=3, cex=1.5, col=colb)

  # add a nice legend with information
  delta.peak <- with( pdelta.df, delta.sek[which(pdelave == max(pdelave))])
  if(!BROB)
  {
    delta.mean <- with( pdelta.df, pdelave %*% delta.sek )
    delta.mean.out <- round(delta.mean,dig)
    mtext(eval(substitute(expression(paste(delta.cMass,"% lower ",delta," = ", delta.min, "  |  ",
                                           delta.cMass,"% upper ",delta," = ", delta.max)),
                          list(delta.cMass=cMass$pdelave.delta.delta.cMass1$cMass,
                               delta.min=round(cMass$pdelave.delta.delta.cMass1$lower,dig),
                               delta.max=round(cMass$pdelave.delta.delta.cMass1$upper,dig) )
    )) , 4, line=2, cex=0.9, col="black", adj=0)
  
    mtext(eval(substitute(expression(paste(delta.cMass,"% lower ",delta," = ", delta.min, "  |  ",
                                           delta.cMass,"% upper ",delta," = ", delta.max)),
                          list(delta.cMass=cMass$pdelave.delta.delta.cMass2$cMass,
                               delta.min=round(cMass$pdelave.delta.delta.cMass2$lower,dig),
                               delta.max=round(cMass$pdelave.delta.delta.cMass2$upper,dig) )
    )) , 4, line=3, cex=0.9, col="black", adj=0)    
  } else
  {
    delta.mean <- sum(brob(pdelta.df[,"pdelave"] + as.brob(pdelta.df[,"delta.sek"])@x))
    delta.mean.out <- paste("exp(",signif(delta.mean@x,digits=dig+2),")",sep="")
  }
  mtext(eval(substitute(expression(paste("mode (peak) ",delta," = ",delta.peak, "  |  mean ",delta," = ",delta.mean.out)),
                        list(delta.peak=round(delta.peak,dig), delta.mean.out=delta.mean.out))) , 4, line=1, cex=0.9, col="black", adj=0)

  # add a nice legend to explain curves and associated probability densities
  par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
  legend("bottom", legend=c(expression(paste("p(",delta," | S, D" ["1"] ,", D" ["2"] ,", I)")),
                            expression(paste("p(",delta," | ",bar(S),", D" ["1"] ,", D" ["2"] ,", I)")),
                            expression(paste("p(",delta," | D" ["1"] ,", D" ["2"] ,", I)"))),
         xpd=TRUE, horiz=TRUE, inset=c(0,0), y.intersp=2.4,
         col=cols[1:3], lty=c(2,3,1), lwd=lwd, bty="n", cex=0.9)               

################################################################################
# p(r|C,D_1,D_2,I) problem
# plot p(r|C,D_1,D_2,I)
  if(!BROB)
  {
    y.axistextSA <- expression(paste("p(r | C, D" ["1"] ,", D" ["2"] ,", I)",sep=""))
  } else y.axistextSA <- expression(paste("p(r | C, D" ["1"] ,", D" ["2"] ,", I) on log() scale",sep=""))
  x.limsSA <- range(pr.df[,"r.sek"])
  y.limsSA <- range(pr.df[,"prA"])*c(1,fac)
  if(by1) par(ask=TRUE)
  par(mar=c(5,6,5,5))
  plot(x.limsSA[2],y.limsSA[2], type="l", lty=1, col="white",
      main="",
      xlab="",
      ylab=y.axistextSA,
      pre.plot=grid(),
      bty="l", xlim=x.limsSA, ylim=y.limsSA)
  if(filling) polygon( c(pr.df[1,"r.sek"],pr.df[,"r.sek"], pr.df[dim(pr.df)[1],"r.sek"],pr.df[1,"r.sek"]),
                       c(y.limsSA[1],pr.df[,"prA"],y.limsSA[1],y.limsSA[1]),
                       col=rgb(cols.rgb[1,1],cols.rgb[2,1],cols.rgb[3,1],0.5), border=NA)
  points(pr.df[,"r.sek"], pr.df[,"prA"], col=cols[1], type="l", lty=1, lwd=lwd)
  mtext(expression(paste(sigma," Standard Deviation Ratio (r)",sep="")), 3, line=2, cex=1.5)
  mtext("assuming the means are the same", 3, line=0.8)
  mtext(expression(paste("r = ",sigma[1],"/",sigma[2])), 1, line=3, cex=1.5, col=colb)

################################################################################
# p(r|Cbar,D_1,D_2,I) problem
# plot p(r|Cbar,D_1,D_2,I)
  if(!BROB)
  {
    y.axistextSB <- expression(paste("p(r | ",bar(C),", D" ["1"] ,", D" ["2"] ,", I)",sep=""))
  } else y.axistextSB <- expression(paste("p(r | ",bar(C),", D" ["1"] ,", D" ["2"] ,", I) on log() scale",sep=""))
  x.limsSB <- range(pr.df[,"r.sek"])
  y.limsSB <- range(pr.df[,"prB"])*c(1,fac)
  if(by1) par(ask=TRUE)
  par(mar=c(5,6,5,5))
  plot(x.limsSB[2], y.limsSB[2], col="white",
       main="",
       xlab="",
       ylab=y.axistextSB,
       pre.plot=grid(),
      bty="l", xlim=x.limsSB, ylim=y.limsSB)
  if(filling) polygon( c(pr.df[1,"r.sek"],pr.df[,"r.sek"], pr.df[dim(pr.df)[1],"r.sek"],pr.df[1,"r.sek"]),
                       c(y.limsSB[1],pr.df[,"prB"],y.limsSB[1],y.limsSB[1]),
                       col=rgb(cols.rgb[1,1],cols.rgb[2,1],cols.rgb[3,1],0.5), border=NA)
  points(pr.df[,"r.sek"], pr.df[,"prB"], col=cols[1], type="l", lty=1, lwd=lwd)

  mtext(expression(paste(sigma," Standard Deviation Ratio (r)",sep="")), 3, line=2, cex=1.5)
  mtext("assuming the means are different", 3, line=0.8)
  mtext(expression(paste("r = ",sigma[1],"/",sigma[2])), 1, line=3, cex=1.5, col=colb)

################################################################################
# p(r|D_1,D_2,I)
# plot p(r|C,D_1,D_2,I), p(r|Cbar,D_1,D_2,I)
  if(!BROB)
  {
    y.axistextSAB <- c("Probability Density")
  } else y.axistextSAB <- c("Probability Density on log() scale")
  x.limsSAB <- range(pr.df[,"r.sek"])
  y.limsSAB <- range(pr.df[,c("prA","prB","prave.r.delta")])*c(1,fac)
  if(by1) par(ask=TRUE)
  par(mar=c(5,6,5,5))
  par(oma=c(2,1,1,1))
  par("cex.axis" = 0.8)
  plot(x.limsSAB[2], y.limsSAB[2], col="white",
       main="",
       xlab="",
       ylab=y.axistextSAB,
       pre.plot=grid(),
       bty="l",
       xlim=x.limsSAB, ylim=y.limsSAB)
  if(filling)
  {
    polygon( c(pr.df[1,"r.sek"],pr.df[,"r.sek"], pr.df[dim(pr.df)[1],"r.sek"],pr.df[1,"r.sek"]),
             c(y.limsSAB[1],pr.df[,"prA"],y.limsSAB[1],y.limsSAB[1]),
             col=rgb(cols.rgb[1,1],cols.rgb[2,1],cols.rgb[3,1],0.5), border=NA)
    polygon( c(pr.df[1,"r.sek"],pr.df[,"r.sek"], pr.df[dim(pr.df)[1],"r.sek"],pr.df[1,"r.sek"]),
             c(y.limsSAB[1],pr.df[,"prB"],y.limsSAB[1],y.limsSAB[1]),
             col=rgb(cols.rgb[1,2],cols.rgb[2,2],cols.rgb[3,2],0.5), border=NA)
    polygon( c(pr.df[1,"r.sek"],pr.df[,"r.sek"], pr.df[dim(pr.df)[1],"r.sek"],pr.df[1,"r.sek"]),
             c(y.limsSAB[1],pr.df[,"prave.r.delta"],y.limsSAB[1],y.limsSAB[1]),
             col=rgb(cols.rgb[1,3],cols.rgb[2,3],cols.rgb[3,3],0.5), border=NA)
  }    
  points(pr.df[,"r.sek"], pr.df[,"prA"], col=cols[1], type="l", lty=1, lwd=lwd)
  points(pr.df[,"r.sek"], pr.df[,"prB"], col=cols[2], type="l", lty=2, lwd=lwd)
  points(pr.df[,"r.sek"], pr.df[,"prave.r.delta"], col=cols[3], type="l", lty=3, lwd=lwd)

  mtext(expression(paste("r Standard Deviation Ratio",sep="")), 3, line=3, cex=1.5)
  mtext("independent of the means being different/ the same", 3, line=1.5)
  mtext(expression(paste("r = ",sigma[1],"/",sigma[2])), 1, line=3, cex=1.5, col=colb)

  # add a nice legend with information
  r.peak <- with( pr.df, r.sek[which(prave == max(prave))])
  if(!BROB)
  {
    r.mean <- with( pr.df, prave %*% r.sek)
    r.mean.out <- round(r.mean,dig)
    mtext(eval(substitute(expression(paste(r.cMass,"% lower r = ", r.min, "  |  ",
                                           r.cMass,"% upper r = ", r.max)),
                          list(r.cMass=cMass$prave.r.delta.cMass1$cMass,
                               r.min=round(cMass$prave.r.delta.cMass1$lower,dig),
                               r.max=round(cMass$prave.r.delta.cMass1$upper,dig) )
    )) , 4, line=2, cex=0.9, col="black", adj=0)
    mtext(eval(substitute(expression(paste(r.cMass,"% lower r = ", r.min, "  |  ",
                                           r.cMass,"% upper r = ", r.max)),
                          list(r.cMass=cMass$prave.r.delta.cMass2$cMass,
                               r.min=round(cMass$prave.r.delta.cMass2$lower,dig),
                               r.max=round(cMass$prave.r.delta.cMass2$upper,dig) )
    )) , 4, line=3, cex=0.9, col="black", adj=0)
  } else
  {
    r.mean <- sum(brob(pr.df[,"prave"] + as.brob(pr.df[,"r.sek"])@x))
    r.mean.out <- paste("exp(",signif(r.mean@x,digits=dig+2),")",sep="")    
  }    
  mtext(eval(substitute(expression(paste("mode (peak) r = ",r.peak, "  |  mean r = ",r.mean.out)),
                        list(r.peak=round(r.peak,dig), r.mean.out=r.mean.out ))) , 4, line=1, cex=0.9, col="black", adj=0)
  # add a nice legend to explain curves and associated probability densities
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend=c(expression(paste("p(r | C, D" ["1"] ,", D" ["2"] ,", I)")),
                            expression(paste("p(r | ",bar(C),", D" ["1"] ,", D" ["2"] ,", I)")),
                            expression(paste("p(r | D" ["1"] ,", D" ["2"] ,", I)"))),
         xpd=TRUE, horiz=TRUE, inset=c(0,0), y.intersp=2.4,
         col=cols[1:3], lty=c(2,3,1), lwd=lwd, bty="n", cex=0.9)
}
# call:
# DiM.plot.pg(DiM.plotvalues.res.nonbrob.brob, filling=TRUE, BROB=TRUE)
################################################################################


################################################################################
# print results / tables
DiM.print.pg <- function(DIM.pg.res, dig=4)
{
   cat("\n##########################################################################")
   cat("\n#\n# On the 'Difference in means'")
   cat("\n\n# original: G.L. Bretthorst (1993)")
   cat("\n# R code after Mathematica code by Phil Gregory")
   cat("\n\n# Descriptive values of samples\n\n")
   print(DIM.pg.res$desc.df, right=FALSE, row.names=FALSE)
   cat("\n# Prior information\n\n")
   print(DIM.pg.res$prior.df, right=FALSE, row.names=FALSE)
   cat("\n# Posterior probabilities\n\n")
   cat("\nTotal Probability:\t",DIM.pg.res$tot.out,"\n\n",sep="")
   print(DIM.pg.res$posterior.df, digits=dig, right=FALSE, row.names=FALSE)
   cat("\n# Odds Ratios\n\n")
   print(DIM.pg.res$OR.res.df, digits=dig, right=FALSE, row.names=FALSE)
   cat("\n##########################################################################")
   cat("\n")
}
# call:
# DiM.print.pg(DIM.pg.res)
################################################################################
 
 
################################################################################
# END OF 'ON THE DIFFERENCE OF MEANS'
################################################################################

