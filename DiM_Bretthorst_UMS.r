
################################################################################
# On the difference in means
# paper: G.L. Bretthorst "On the difference of means" (1993)

#
# R code based on Mathematica code by Urban Studer (90's, Zürich/ CH)
# R code by Leo G
# first: 12-06-05
# latest: 21-06-06, 20-04-17, 07-09-20

# additional note:
# *- we use the sdame variable names and labels like UM Studer to allow comparison
#    with his Mathematica code. This is different from the PG script where the
#    labels from P Gregory are used.

################################################################################
# success rates and integration bounds based on successes and total N
SucRatesIntBounds <- function(Si, Ni, Sii, Nii, smin, snames=c("sample1","sample2"))
{
  # --------------------------------------------------
  # Success rates and integration bounds in the case
  # of the (conservative) Bayes-Laplace prior
  # --------------------------------------------------
  # necessary variables:
  # N, S(uccesses) for sample 1 (i) and 2 (ii)
  
  ########################
  # defintion of constants
  
  # based on frequencies (successes + failures)
  # Di = mean group 1
  # si = standard deviation group 1
  # Dii = mean group 2
  # sii = standard deviation group 2
  
  Di <- (Si + 1) / (Ni + 2)
  si <- sqrt(Di * (1 - Di) / (Ni + 3))
  Dii <- (Sii + 1) / (Nii + 2)
  sii <- sqrt(Dii * (1 - Dii) / (Nii + 3))
  
  ########################
  # calculation boundaries
  Nmax <- max(Ni, Nii)
  Nmin <- min(Ni, Nii)
  
  # boundaries standard deviations
  sL <- floor(1000 * sqrt((Nmax+1) / (Nmax+3)) / (Nmax + 2)) / 1000
  sH <- ceiling(1000 / (2 * sqrt(Nmin + 3))) / 1000
  
  # boundaries means
  L <- floor(100 * (smin + 1) / (Nmax + 2)) / 100
  H <- 1 - L
  
  res <- list(Si=Si, Ni=Ni, Sii=Sii, Nii=Nii, smin=smin,
              Di=Di, si=si, Dii=Dii, sii=sii,
              L=L, H=H, sL=sL, sH=sH,
              snames=snames)
  attr(res,"typ") <- c("SRIB")
  return(res)
}
#call:
#res.SIB <- SucRatesIntBounds(Si=11, Ni=15, Sii=10, Nii=16, smin=0, snames=c("voluntary","non-voluntary"))


################################################################################  
#originally taken from Bolstad2 and tweaked to work with Brobdingnag
#DEV BROB
simpsonrule.brob <- function(fx, sL, sH, Nsteps=100)
{
  # taken from sintegral from Bolstad2
  sek <- seq(sL,sH,length=Nsteps)
  l.intv <- 2*Nsteps+1
  intv.x <- approx(sek,sek,n=l.intv)$x
  h <- diff(intv.x)[1]
  inity <- as.list(sapply(seq_along(1:l.intv), function(x) fx(intv.x[x])))
  inity1 <- inity[2 * (1:Nsteps) - 1]
  inity2 <- inity[2 * (1:Nsteps)]
  inity3 <- inity[2 * (1:Nsteps) + 1]
  sum <- 0
  for(i in inity1) sum <- sum + unlist(i)
  for(i in inity2) sum <- sum + 4*unlist(i)
  for(i in inity3) sum <- sum + unlist(i)
  sum <- sum*(h/3)
return(sum)
}

################################################################################
#DEV
#Gamma funcs
#(generalized) incomplete gamma function
gamma_inc_gen <- function(a,z0,z1=NA, log=FALSE)
{
  # incomplete
  # equivalent to Mathematica N[Gamma[a,z0]]  
  
  # generalized incomplete  
  # equivalent to Mathematica N[Gamma[a,z0,z1]] = G[a, z0] - G[a, z1]
  
  # comparison to Mathematica
  
  # G[z0] = integral_0^Inf t^(z0-1) exp(-t) dt
  # G[z0] is gamma(z0)
  
  # G[a, z0] = 1/Gamma(a) integral_z0^Inf t^(a-1) exp(-t) dt
  # G[a, z0] is pgamma(z0, a, lower=FALSE)
  
  # 1 - G[a, z0] = 1/Gamma(a) integral_0^z0 t^(a-1) exp(-t) dt
  # 1 - G[a, z0] is pgamma(z0, a, lower=TRUE) #default
  
  # a = real part
  # z0 = upper limit
  # z1 = lower limit   
  
  if(!log)
  {
    if(is.na(z1)) return( pgamma(z0,a, lower=FALSE) )
    else return( pgamma(z0,a, lower=FALSE) - pgamma(z1,a, lower=FALSE) )
  }  
  else
  {
    library(Brobdingnag)
    if(is.na(z1)) return( pgamma(z0,a, lower=FALSE, log=TRUE))
    else return( brob(pgamma(z0,a, lower=FALSE,log=TRUE)) - brob(pgamma(z1,a, lower=FALSE,log=TRUE)) )
  }  
}


################################################################################
# --------------------------------------------------
# ON THE DIFFERENCE IN MEANS
# --------------------------------------------------
DiffinMeans <- function(inval=NULL, out=FALSE, smin=0, BROB=FALSE, Nsteps=100)
{
  ################################################################################
  # necessary variables in the function
  #{ NN, DD, Dsi, Dsii, DsD, ss,
  #  dd, lownum, upnum, low, up, pmv, PMV,
  #  zz, lowinum, upinum, lowiinum, upiinum, pmbarv, PMbarV,
  #  pmvbar, PMVbar, pmbarvbar, PMbarVbar, cc,
  #  mv, mbarv, mvbar, mbarvbar,
  #  samemeans, diffmeans, samevars, diffvars, diffsets },
  
  ################################################################################
  # extract values for calculations from inputv object (class "UMS")
  
  #for large numbers
  if(BROB)
  {
    library(Brobdingnag)
    cat("\nUsing log() and package 'Brobdingnag'\nAll results are expressed as log(RESULT)\n")
  }  
  
  Di <- inval[["Di"]]
  Dii <- inval[["Dii"]]
  si <- inval[["si"]]
  sii <- inval[["sii"]]
  Ni <- inval[["Ni"]]
  Nii <- inval[["Nii"]]
  L <- inval[["L"]]
  H <- inval[["H"]]
  sL <- inval[["sL"]]
  sH <- inval[["sH"]]
  snames <- inval[["snames"]]
  smin <- inval[["smin"]]
  
  ################################################################################
  # defintion of constants
  NN <- Ni+Nii
  DD <- (Ni * Di + Nii * Dii) / NN
  Dsi <- (Ni-1) / Ni * si^2 + Di^2
  Dsii <- (Nii-1) / Nii * sii^2 + Dii^2
  DsD <- (Ni * Dsi + Nii * Dsii) / NN
  ss <- sqrt(NN * (DsD - DD^2) / (NN-1))
  
  cat("\nL - Mean_comb < 0 : ",L < DD," [comparison L < DD]\n")
  if(L < DD)
  {
    cat("'+'-sign between Gamma-factors is ok\n\n")
  }  else cat("'+'-sign between Gamma-factors is false!\n\n")
    
################################################################################
# p_mv
  dd <- NN * (DsD - DD^2)
  lownum <- NN * (L-DD)^2
  upnum  <- NN * (H-DD)^2
  
  #low = lownum / (2*s^2)
  #up  = upnum / (2*s^2)
  PMV.hypo <- function(NN,dd,upnum,lownum,sL,sH)
  {
    integpmv <- function(s)
    { 
      1 / (s^NN) * exp(-dd / (2 * s^2)) *
      ( pgamma(upnum/(2*s^2),1/2)*gamma(1/2) + pgamma(lownum/(2*s^2),1/2)*gamma(1/2) )
    }
    pmv <- integrate(integpmv, lower=sL, upper=sH)$value
    #cat("\npmv =\t",pmv,"\n")
    PMV <- pmv / sqrt(2*NN)
    return(PMV)  
  }  
  #call:
  #PMV <- PMV.hypo(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH)

# BROB
  PMV.hypo.brob <- function(NN,dd,upnum,lownum,sL,sH,Nsteps=100)
  {
    integpmv.brob <- function(s)
    {
      1 / (as.brob(s)^NN) * exp(as.brob(-dd) / (2 * s^2)) *
        ( gamma_inc_gen(1/2,0,upnum/(2*s^2),log=TRUE)*sqrt(pi) + gamma_inc_gen(1/2,0,lownum/(2*s^2),log=TRUE)*sqrt(pi) )
    }
    pmv.brob <- simpsonrule.brob(fx=integpmv.brob, sL=sL, sH=sH, Nsteps=Nsteps)
    PMV.brob <- pmv.brob / sqrt(2*NN)
    return(PMV.brob)
  }
  #call:
  #PMV.brob <- PMV.hypo.brob(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH, Nsteps=100)

################################################################################
# p_mbarv
  zz <- Ni * (Dsi - Di^2) + Nii * (Dsii - Dii^2)
  lowinum <- Ni * (L - Di)^2
  upinum <- Ni * (H - Di)^2
  lowiinum <- Nii * (L - Dii)^2
  upiinum <- Nii * (H - Dii)^2
  
  #Bretthorst, p.6f.
  #sign in inner brackets between the two pgamma*gamma's
  #sign + in case if U1H, H1L are of different sign
  #sign - in case if U1H, U1L are of same sign
  #same for U2H, U2L
  PMbarV.hypo <- function(NN,zz,dd,upinum,lowinum,sL,sH,H,L,Ni,Nii)
  {
    integpmbarv <- function(s)
    {
      #(* + if (L-DD) < 0 *) 
      1/(s^(NN-1)) * exp(-zz/(2*s^2)) *
      ( pgamma(upinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowinum/(2*s^2),1/2)*gamma(1/2) ) *
      ( pgamma(upiinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowiinum/(2*s^2),1/2)* gamma(1/2) )
    }
    pmbarv <- integrate(integpmbarv, lower=sL, upper=sH)$value
    PMbarV <- pmbarv / (2*(H-L) * sqrt(Ni*Nii))
    return(PMbarV)  
  }  
  #call:
  #PMbarV <- PMbarV.hypo(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)

# BROB
  PMbarV.hypo.brob <- function(NN,zz,dd,upinum,lowinum,sL,sH,H,L,Ni,Nii,Nsteps=100)
  {
    integpmbarv.brob <- function(s)
    {
      1/(as.brob(s)^(NN-1)) * exp(as.brob(-zz)/(2*s^2)) *
      ( pgamma(upinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowinum/(2*s^2),1/2)*gamma(1/2) ) *
      ( pgamma(upiinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowiinum/(2*s^2),1/2)*gamma(1/2) )
    }
    pmbarv.brob <- simpsonrule.brob(fx=integpmbarv.brob, sL=sL, sH=sH, Nsteps=Nsteps)
    PMbarV.brob <- pmbarv.brob / (2*(H-L) * sqrt(Ni*Nii))
  return(PMbarV.brob)  
  }  
  #call:
  #PMbarV.brob <- PMbarV.hypo.brob(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=100)

################################################################################
# p_mvbar

# note Mathematica: gamma[a,z0,z1] = gamma[a,z1] - gamma[a,z0]
  UiA <- function(A) Ni*(Dsi-2*Di*A+A^2)/2
  UiiA <- function(A) Nii*(Dsii-2*Dii*A+A^2)/2

  PMVbar.hypo <- function(Dsi,Dsii,sL,sH,H,L,Ni,Nii)
  {
    integpmvbar <- function(A)
    {
      1/UiA(A)^(Ni/2) *
      1/UiiA(A)^(Nii/2) *
      ( pgamma(UiA(A)/(sH^2),Ni/2)*gamma(Ni/2) - pgamma(UiA(A)/(sL^2),Ni/2)*gamma(Ni/2) ) *
      ( pgamma(UiiA(A)/(sH^2),Nii/2)*gamma(Nii/2) - pgamma(UiiA(A)/(sL^2),Nii/2)*gamma(Nii/2) )
    }
    pmvbar <- integrate(integpmvbar, lower=L, upper=H)$value
    PMVbar <- pmvbar / (4*log(sH/sL))
    return(PMVbar)  
  }  
  #call:
  #PMVbar <- PMVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)

# not possible for pgamma due to its definition without brob numbers  
  UiA.brob <- function(A) as.brob( Ni*(Dsi-2*Di*A+A^2)/2 )
  UiiA.brob <- function(A) as.brob( Nii*(Dsii-2*Dii*A+A^2)/2 )
# BROB
  PMVbar.hypo.brob <- function(Dsi,Dsii,sL,sH,H,L,Ni,Nii,Nsteps=100)
  {
    integpmvbar.brob <- function(A)
    {
      1/as.brob(UiA.brob(as.brob(A))^(as.brob(Ni)/2)) *
      1/as.brob(UiiA.brob(as.brob(A))^(as.brob(Nii)/2)) *
      ( brob( (pgamma(UiA(A)/(sH^2),Ni/2,log.p=TRUE,lower=FALSE)+lgamma(Ni/2)) ) -
        brob( (pgamma(UiA(A)/(sL^2),Ni/2,log.p=TRUE,lower=FALSE)+lgamma(Ni/2)) ) ) *
      ( brob( (pgamma(UiiA(A)/(sH^2),Nii/2,log.p=TRUE,lower=FALSE)+lgamma(Nii/2)) ) - 
        brob( (pgamma(UiiA(A)/(sL^2),Nii/2,log.p=TRUE,lower=FALSE)+lgamma(Nii/2)) )
      )
    }
    pmvbar.brob <- simpsonrule.brob(fx=integpmvbar.brob, sL=L, sH=H, Nsteps=Nsteps)
    PMVbar.brob <- pmvbar.brob / (4*log(sH/sL))    
    return(PMVbar.brob)
  }
  #call:
  #PMVbar.brob <- PMVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=100)

################################################################################
# p_mbarvbar
  
  #Bretthorst, p.11
  #sign: minus if U1H, U1L have same sign
  #sign: plus if U1H, U1L have different sign
  #same true for U2H, U2L
  
  UiiB <- function(B) Nii*(Dsii-2*Dii*B+B^2)/2
  PMbarVbar.hypo <- function(Dsi,Dsii,sL,sH,H,L,Ni,Nii)
  {
    integpmbarvbar1 <- function(A)
    {
      1/(UiA(A)^(Ni/2)) *
        #lower=FALSE for integral z^inf
        #the absolute value is ok, only the sign would differ
        #the multiplication of both integrals if handled identical
        #will lead to positive values
        #and with prob distributions we do not have negative probs!
      ( pgamma(UiA(A)/(sH^2),Ni/2,lower=FALSE)*gamma(Ni/2) -
        pgamma(UiA(A)/(sL^2),Ni/2,lower=FALSE)*gamma(Ni/2) )
    }
    
    integpmbarvbar2 <- function(B)
    {
      1/(UiiB(B)^(Nii/2)) *
        #lower=FALSE for integral z^inf
      ( pgamma(UiiB(B)/(sH^2),Nii/2,lower=FALSE)*gamma(Nii/2) -
        pgamma(UiiB(B)/(sL^2),Nii/2,lower=FALSE)*gamma(Nii/2) )
    }

    pmbarvbar1 <- integrate(integpmbarvbar1, lower=L, upper=H)$value
    pmbarvbar2 <- integrate(integpmbarvbar2, lower=L, upper=H)$value
    pmbarvbar <- pmbarvbar1 * pmbarvbar2
    PMbarVbar <- pmbarvbar / (4*(H-L) * log(sH/sL))
    return(PMbarVbar)  
  }  
  #call:
  #PMbarVbar <- PMbarVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)

# BROB
  UiiB.brob <- function(B) as.brob( Nii*(Dsii-2*Dii*B+B^2)/2 )
  PMbarVbar.hypo.brob <- function(Dsi,Dsii,sL,sH,H,L,Ni,Nii,Nsteps=100)
  {
    #first integral
    integpmbarvbar1.brob <- function(A)
    {
      1/( UiA(as.brob(A))^as.brob(Ni/2) )  *
      ( brob( pgamma(UiA(A)/(sH^2),Ni/2,log.p=TRUE,lower=TRUE)+lgamma(Ni/2) ) -
        brob( pgamma(UiA(A)/(sL^2),Ni/2,log.p=TRUE,lower=TRUE)+lgamma(Ni/2) ) )  
    }  
    #second integral
    integpmbarvbar2.brob <- function(B)
    {
      1/( UiiB.brob(B)^(Nii/2)) *
      ( brob(pgamma(UiiB(B)/(sH^2),Nii/2,log.p=TRUE,lower=TRUE)+lgamma(Nii/2)) -
        brob(pgamma(UiiB(B)/(sL^2),Nii/2,log.p=TRUE,lower=TRUE)+lgamma(Nii/2)) )
    }
    pmbarvbar1.brob <- simpsonrule.brob(fx=integpmbarvbar1.brob, sL=L, sH=H, Nsteps=Nsteps)
    pmbarvbar2.brob <- simpsonrule.brob(fx=integpmbarvbar2.brob, sL=L, sH=H, Nsteps=Nsteps)
    pmbarvbar.brob <- pmbarvbar1.brob * pmbarvbar2.brob
    PMbarVbar.brob <- pmbarvbar.brob / (4*(H-L) * log(sH/sL))
    return(PMbarVbar.brob)
  }
  #call:
  #PMbarVbar.brob <- PMbarVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=100)
  
#compile results independent from BROB
#inputs
  constants <- data.frame(sample1=snames[1],sample2=snames[2],
                          Di,Dii,Ni,si,Nii,sii,
                          NN,DD,ss,L,H,sL,sH,Dsi,Dsii,
                          DsD,dd,
                          lownum,upnum,lowinum,upinum,lowiinum,upiinum,
                          smin)
  #resulting table/ dataframe with input values
  input.df <- data.frame("No." = c(Ni,Nii,NN),
                         "Standard Deviation" = c(si,sii,ss),
                         "Mean" = c(Di,Dii,DD),
                         "Data set" = c(snames,"combined"),
                         check.names = FALSE)
  #table/ dataframe with prior values/ information
  prior.df <- data.frame("Numerical Example" = c("Prior Mean lower bound","Prior Mean upper bound",
                                                 "Prior Standard Deviation lower bound","Prior Standard Deviation upper bound"),
                         "Value" = c(L,H,sL,sH), check.names = FALSE) 

  if(!BROB)
  {
    cat("\nCalculate PMV\n")
    PMV <- PMV.hypo(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH)
    cat("Calculate PMbarV\n")
    PMbarV <- PMbarV.hypo(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
    cat("Calculate PMVbar\n")
    PMVbar <- PMVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
    cat("Calculate PMbarVbar\n")
    PMbarVbar <- PMbarVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
  } else
  {
    cat("\nCalculate PMV\n")
    PMV <- PMV.hypo.brob(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH, Nsteps=Nsteps)
    cat("Calculate PMbarV\n")
    PMbarV <- PMbarV.hypo.brob(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)
    cat("Calculate PMVbar\n")
    PMVbar <- PMVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)
    cat("Calculate PMbarVbar\n")
    PMbarVbar <- PMbarVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)
    # calculate total probability (denominator of Bayes' Theorem)
    # sum of all four hypotheses (probabilities)
  }
    
    cat("Compile results\n")
    cc <- 1 / (PMV + PMbarV + PMVbar + PMbarVbar)
    
    # no *100%, just *1
    percfac <- 1
    mv <- percfac * cc * PMV
    mbarv <- percfac * cc * PMbarV
    mvbar <- percfac * cc * PMVbar
    mbarvbar <- percfac * cc * PMbarVbar
    
    samemeans <- mv + mvbar
    diffmeans <- mbarv + mbarvbar
    samevars <- mv + mbarv
    diffvars <- mvbar + mbarvbar
    diffsets <- mvbar + mbarv + mbarvbar
    
    OR.diffmeans <- diffmeans/samemeans
    OR.samemeans <- samemeans/diffmeans
    OR.diffvars <- diffvars/samevars
    OR.samevars <- samevars/diffvars
    OR.diffsets <- diffsets/mv
    OR.samesets <- mv/diffsets    
    
    #unnormalized + total prob
    posteriorprobs.un <- list(PMV,PMbarV,PMVbar,PMbarVbar,cc)
    names(posteriorprobs.un) <- c("PMV","PMbarV","PMVbar","PMbarVbar","cc")
    
    #normalized
    posteriorprobs.no <- list(mv,mbarv,mvbar,mbarvbar,
                              samemeans,diffmeans,samevars,diffvars,
                              samesets=1-diffsets,diffsets)
    names(posteriorprobs.no) <- c("mv","mbarv","mvbar","mbarvbar",
                                  "samemeans","diffmeans","samevars","diffvars",
                                  "samesets","diffsets")
    
    #odds ratios
    ORs <- list(OR.diffmeans,OR.samemeans,
                OR.diffvars,OR.samevars,
                OR.diffsets,OR.samesets)
    names(ORs) <- c("OR.diffmeans","OR.samemeans",
                    "OR.diffvars","OR.samevars",
                    "OR.diffsets","OR.samesets")
    
    cat("Create dataframes for output\n")
    
    #table/ dataframe with resulting probabilities
    prob.res.df <- data.frame("Hypothesis (Probabilities)" = c(
        "Same Mean,       same Standard Deviation",
        "Different Means, same Standard Deviation",
        "Same Mean,       different Standard Deviations",
        "Different Means, different Standard Deviations",
        "The Means are the same",
        "The Means are different",
        "The Standard Deviations are the same",
        "The Standard Deviations are different",
        "The Data Sets are the same",
        "The Data Sets are different"),
        "abbrev" = c("C&S","Cbar,S","C,Sbar","Cbar&Sbar","C","Cbar","S","Sbar","C&S","Cbar|Sbar"),
        "p(H|D1,D2,I)"=NA,
        sign=NA,
        check.names=FALSE)
      
    #table/ dataframe with odds ratio results
    OR.res.df <- data.frame("Hypothesis in favor of ..." = c(
        "A difference in Means",
        "The same Means",
        "A difference in Standard Deviations",
        "The same Standard Deviations",
        "A difference in the Sets (different Means and/ or Standard Deviations)",
        "The same Sets (same Means and/ or Standard Deviations)"),
        "Odds Ratio (OR)"=NA,
        sign=NA,
        check.names = FALSE)
    
    #ifelse(BROB==TRUE, data.frame( unlist(lapply(ORs,function(x) x@x)) ), unlist(ORs))
    if(!BROB)
    {
      prob.res.df[,3] <- unlist(posteriorprobs.no)
      OR.res.df[,2] <- unlist(ORs)
    } else
    {
      colnames(prob.res.df)[3] <- paste(colnames(prob.res.df)[3]," = exp(x)", sep="")
      prob.res.df[,3] <- unlist(lapply(posteriorprobs.no,function(x) x@x))
      #sign
      prob.res.df[,4] <- unlist(lapply(posteriorprobs.no,function(x) x@positive))
      colnames(OR.res.df)[2] <- paste(colnames(OR.res.df)[2]," = exp(x)", sep="")
      OR.res.df[,2] <- unlist(lapply(ORs,function(x) x@x))
      #sign
      OR.res.df[,3] <- unlist(lapply(ORs,function(x) x@positive))
    }
    
    #results
    res <- list(consts = constants,
                pp.un = posteriorprobs.un,
                pp.no = posteriorprobs.no,
                OR = ORs,
                iv.df = input.df,
                prior.df = prior.df,
                prob.df = prob.res.df,
                OR.df = OR.res.df,
                BROB=BROB)

# print tables
  if(out)
  {
    if(!BROB)
    {
      cat(paste("\nShort output of 'the difference in means':\n\n"))
    } else cat(paste("\nShort output of 'the difference in means' based on BROBs:\n\n"))
    
    print(input.df, right=FALSE, row.names=FALSE)
    cat(paste("\n"))
    print(prior.df, right=FALSE, row.names=FALSE)
    cat(paste("\n"))
    print(prob.res.df, right=FALSE, row.names=FALSE)
    cat(paste("\n"))
    print(OR.res.df, right=FALSE, row.names=FALSE)
    cat(paste("\n"))    
  }

return(res)
}
#call:
#DiM.results <- DiffinMeans(inval=res.SIB, out=FALSE, BROB=FALSE)


################################################################################
# print nice output
UMSprint <- function(results, SRIB=NA, dig=6)
{
  
  cat("\n\n#########################################################################################\n")
  cat("###\n")
  cat("### ON THE DIFFERENCE IN MEANS\n")
  cat("### G.L. Bretthorst (1993)\n")
  cat("###\n### original Mathematica code by U.M. Studer (90s, Switzerland)\n")
  cat("\n\nNote:\nIf any probability is printed as '1' (= one) or '0' (= zero), it means that the\nprobability is practically that value by giving respect to limited computer precision.\n")
  BROB <- results[["BROB"]]
  consts <- results[["consts"]]
  
  if(!any(is.na(SRIB)) && attr(SRIB,"typ") == "SRIB") srib.base <- TRUE else srib.base <- FALSE

  pp.un.temp <- results[["pp.un"]]
  pp.no.temp <- results[["pp.no"]]
  
  term.diffmeansvssame <- c("different")
  term.diffvarsvssame <- c("different")
  term.diffsetsvssame <- c("different means and/ or")
  
  if(!BROB)
  {
    pp.un.temp <- lapply(pp.un.temp, function(x) signif(x,digits=dig))
    pp.no.temp <- lapply(pp.no.temp, function(x) signif(x,digits=dig))
    
    diffmeansvssame <- signif(pp.no.temp$diffmeans/pp.no.temp$samemeans, digits=dig)
    if(pp.no.temp$diffmeans < pp.no.temp$samemeans)
    {
      diffmeansvssame <- 1/diffmeansvssame
      term.diffmeansvssame <- c("the same")
    }  
    
    diffvarsvssame <- signif(pp.no.temp$diffvars/pp.no.temp$samevars, digits=dig)
    if(pp.no.temp$diffvars < pp.no.temp$samevars)
    {
      diffvarsvssame <- 1/diffvarsvssame
      term.diffvarsvssame <- c("the same")
    }  
    
    diffsetsvssame <- signif(pp.no.temp$diffsets/pp.no.temp$mv, digits=dig)
    if(pp.no.temp$diffsets < pp.no.temp$mv)
    {
      diffsetsvssame <- 1/diffsetsvssame
      term.diffsetsvssame <- c("the same means and")
    }  
    
  } else
  {
    
    #convert back
    pp.un.expcheck <- unlist(lapply(pp.un.temp, function(x) exp(x[]@x)))
    pp.no.expcheck <- unlist(lapply(pp.no.temp, function(x) exp(x[]@x)))
    #check for non-infinite and exact-zero values that should replace existent values (BROBs)
    pp.un.IDs <- as.numeric( which(!is.infinite(pp.un.expcheck) & pp.un.expcheck != 0 ))
    pp.no.IDs <- as.numeric( which(!is.infinite(pp.no.expcheck) & pp.no.expcheck != 0 ))
    
    #replace BROBs if non-infinite-converted values are present
    if(length(pp.un.IDs) > 0)
    {
      pp.un.temp[pp.un.IDs] <- signif( pp.un.expcheck[pp.un.IDs], digits=dig)
    }
    pp.un.ID2 <- which(!(1:length(pp.un.expcheck)) %in% pp.un.IDs == TRUE)
    if(length(pp.un.ID2) > 0) pp.un.temp[pp.un.ID2] <- unlist( lapply(seq_along(pp.un.ID2), function(x) paste( ifelse(pp.un.temp[pp.un.ID2[x]][[1]]@positive,"+","-"), "exp(",signif( pp.un.temp[pp.un.ID2[x]][[1]]@x, digits=dig),")",sep="") ) )
    
    if(length(pp.no.IDs) > 0)
    {
      pp.no.temp[pp.no.IDs] <- signif( pp.no.expcheck[pp.no.IDs], digits=dig)
    }
    pp.no.ID2 <- which(!(1:length(pp.no.expcheck)) %in% pp.no.IDs == TRUE)
    if(length(pp.no.ID2) > 0) pp.no.temp[pp.no.ID2] <- unlist( lapply(seq_along(pp.no.ID2), function(x) paste( ifelse(pp.no.temp[pp.no.ID2[x]][[1]]@positive,"+","-"), "exp(",signif( pp.no.temp[pp.no.ID2[x]][[1]]@x, digits=dig),")",sep="") ) )
    
    pp.un.temp
    pp.no.temp
    
    #diffmeans vs. samemeans
    if(results[["pp.no"]]$diffmeans > results[["pp.no"]]$samemeans)
    {
      diffmeansvssame <- results[["pp.no"]]$diffmeans / results[["pp.no"]]$samemeans
    } else
    {
      diffmeansvssame <- results[["pp.no"]]$samemeans / results[["pp.no"]]$diffmeans
      term.diffmeansvssame <- c("the same")
    }  
      
    if(is.infinite(exp(diffmeansvssame@x)))
    {  
      diffmeansvssame <- paste("exp(",signif(diffmeansvssame@x,digits=dig),")",sep="")
    } else diffmeansvssame <- round(exp(diffmeansvssame@x),digits=dig)
    
    #diffvars vs. samevars
    if(results[["pp.no"]]$diffvars > results[["pp.no"]]$samevars)
    {
      diffvarsvssame <- results[["pp.no"]]$diffvars / results[["pp.no"]]$samevars
    } else
    {
      diffvarsvssame <- results[["pp.no"]]$samevars / results[["pp.no"]]$diffvars
      term.diffvarsvssame <- c("the same")
    }

    if(is.infinite(exp(diffvarsvssame@x)))
    {  
      diffvarsvssame <- paste("exp(",signif(diffvarsvssame@x,digits=dig),")",sep="")
    } else diffvarsvssame <- round(exp(diffvarsvssame@x),digits=dig)
    
    #diffsets vs. samesets
    if(results[["pp.no"]]$diffsets > results[["pp.no"]]$mv)
    {
      diffsetsvssame <- results[["pp.no"]]$diffsets / results[["pp.no"]]$mv
    } else
    {
      diffsetsvssame <- results[["pp.no"]]$mv / results[["pp.no"]]$diffsets
      term.diffsetsvssame <- c("the same means and")
    }

    if(is.infinite(exp(diffsetsvssame@x)))
    {  
      diffsetsvssame <- paste("exp(",signif(diffsetsvssame@x,digits=dig),")",sep="")
    } else diffsetsvssame <- round(exp(diffsetsvssame@x),digits=dig)
    
  }    
  
  # data (input) = descriptive statistics
  cat("\n-------------------------------- Data (Input) ------------------------------------------\n\n")
  if(!srib.base)
  {
    cat(paste("N_1 = ",consts$Ni ," :\t\t\tMean_1 ± SD_1\t\t= ", round(consts$Di,digits=dig)," ± ",round(consts$si,digits=dig), "\n", sep=""))
    cat(paste("N_2 = ",consts$Nii," :\t\t\tMean_2 ± SD_2\t\t= ", round(consts$Dii,digits=dig)," ± ",round(consts$sii,digits=dig), "\n", sep=""))
  } else
  {
    cat(paste("s_1/N_1 = ",SRIB$Si,"/",consts$Ni," :\t\t\tMean_1 ± SD_1\t\t= ", round(consts$Di,digits=dig)," ± ",round(consts$si,digits=dig), "\n", sep=""))
    cat(paste("s_2/N_2 = ",SRIB$Sii,"/",consts$Nii," :\t\t\tMean_2 ± SD_2\t\t= ", round(consts$Dii,digits=dig)," ± ",round(consts$sii,digits=dig), "\n", sep=""))
  }
  cat(paste("N_total = N_1 + N_2  = ", consts$NN ," :\t\tMean_comb ± SD_comb\t= ", round(consts$DD,digits=dig)," ± ",round(consts$ss,digits=dig), "\n", sep=""))
  cat("\n")
  cat(paste("Bounds on the Mean (s_min = ",consts$smin,"):\t\tMean_L = ",consts$L, ",\tMean_H = ",consts$H,"\n",sep=""))
  cat(paste("Bounds on the Standard Deviation:\t  SD_L = ",consts$sL,",\t  SD_H = ",consts$sH,"\n",sep=""))
  if(consts$L < consts$DD)
  {
    cat(paste("\nMean_L - Mean_comb < 0 = ", (consts$L<consts$DD), "\t(-> '+'-sign between Gamma-fcts o.k.)","\n", sep=""))
  } else cat(paste("\nMean_L - Mean_comb < 0 = ", (consts$L<consts$DD), "\t(-> '+'-sign between Gamma-fcts false!)","\n", sep=""))
  
  # results
  cat("\n-------------------------------- Results ------------------------------------------------\n\n")
  
  if(BROB) cat("...based on BROBs - try to convert back to non-BROB numbers...\n\n")
  
  cat(paste("p(mv | D_1, D_2, I)\t\t= const. ", pp.un.temp$PMV, "\n", sep=""))
  cat(paste("p(mbarv | D_1, D_2, I)\t\t= const. ", pp.un.temp$PMbarV, "\n", sep=""))
  cat(paste("p(mvbar | D_1, D_2, I)\t\t= const. ", pp.un.temp$PMVbar, "\n", sep=""))
  cat(paste("p(mbarvbar | D_1, D_2, I)\t= const. ", pp.un.temp$PMbarVbar,"\n\n", sep=""))
  cc.x <- signif(1/(4*(consts$H-consts$L) * log(consts$sH/consts$sL) * (2*pi)^(consts$NN/2)),digits=dig)
  cat(paste("where\t\tconst.\t= ", cc.x," / p(D_1,D_2|I)\n\t\t\t= ",pp.un.temp$cc, "\n",sep=""))
  
  # model probabilities
  cat(paste("\n--------------- Model --------------------------------- Probability ---------------------\n\n", sep=""))
  cat(paste("mv:\t\tSame Mean,      Same Variance:\t\t", pp.no.temp$mv, "\n", sep=""))
  cat(paste("mbarv:\t\tDifferent Mean, Same Variance:\t\t", pp.no.temp$mbarv, "\n", sep=""))
  cat(paste("mvbar:\t\tSame Mean,      Different Variance:\t", pp.no.temp$mvbar, "\n", sep=""))
  cat(paste("mbarvbar:\tDifferent Mean, Different Variance:\t", pp.no.temp$mbarvbar, "\n", sep=""))
  
  # odds ratios
  cat("\n------------------------------ Odds Ratios ----------------------------------------------\n\n")
  cat(paste("The probability the means are the same is:  ",pp.no.temp$samemeans , "\n", sep=""))
  cat(paste("The probability the means are different is: ",pp.no.temp$diffmeans , "\n", sep=""))
  cat(paste("The odds ratio is ",signif(diffmeansvssame,dig)," to 1 in favor of ",term.diffmeansvssame," means.\n", sep=""))
  cat("\n")
  cat(paste("The probability the variances are the same is:  ",pp.no.temp$samevars,
            "\nThe probability the variances are different is: ",pp.no.temp$diffvars, "\n", sep=""))
  cat(paste("The odds ratio is ",signif(diffvarsvssame,dig)," to 1 in favor of ",term.diffvarsvssame," variances\n", sep=""))
  cat("\n")
  cat(paste("The probability the data sets are the same is:  ",pp.no.temp$mv , "\n", sep=""))
  cat(paste("The probability the data sets are different is: ",pp.no.temp$diffsets , "\n", sep=""))
  cat(paste("The odds ratio is ",signif(diffsetsvssame,dig)," to 1 in favor of ",term.diffsetsvssame," variances.\n", sep=""))
  cat("\n-------------------------------- End ----------------------------------------------------\n\n")
  cat("\n#########################################################################################\n\n")
  
}
################################################################################
#call:
#UMSprint(results=DiM.results)
#UMSprint(results=DiM.results, SRIB=SRIB.results)


################################################################################
# plot graphical comparison
# only qualitative Bretthorst
UMSplot <- function(inval, dig=4, pdfout=FALSE, fname="UMSplot.pdf", loga=TRUE, fac=1.1, Nsteps=100)
{
  
  Si <- inval[["Si"]] 
  Sii <- inval[["Sii"]]
  Di <- inval[["Di"]]
  Dii <- inval[["Dii"]]
  si <- inval[["si"]]
  sii <- inval[["sii"]]
  Ni <- inval[["Ni"]]
  Nii <- inval[["Nii"]]
  L <- inval[["L"]]
  H <- inval[["H"]]
  sL <- inval[["sL"]]
  sH <- inval[["sH"]]
  snames <- inval[["snames"]]
  smin <- inval[["smin"]]
  
  # determining plot range
  sigma <- min(si, sii)
  plr <- round(10 / (sqrt(2 * pi) * sigma) + 5) / 10
  sek <- seq(from=0, to=1, by=1/Nsteps)
  
  # check whether to use the log() to plot (factorials!) = default
  if(loga)
  {
    consti.l <- lfactorial(Ni + 1) - ( lfactorial(Si) + lfactorial(Ni - Si) )
    constii.l <- lfactorial(Nii + 1) - ( lfactorial(Sii) + lfactorial(Nii - Sii) )
    pBLi.l <- function(x) { x^Si * (1 - x)^(Ni - Si) }
    probs.i <- exp(consti.l + log(pBLi.l(sek)))
    pBLii.l <- function(x) { x^Sii * (1 - x)^(Nii - Sii) }
    probs.ii <- exp(constii.l + log(pBLii.l(sek)))
    probs.i
    probs.ii
  } else {
    # defining functions for later plotting
    consti <- factorial(Ni + 1) / ( factorial(Si) * factorial(Ni - Si) )
    pBLi <- function(x) { consti * x^Si * (1 - x)^(Ni - Si) }
    constii <- factorial(Nii + 1) / ( factorial(Sii) * factorial(Nii - Sii) )
    pBLii <- function(x) { constii * x^Sii * (1 - x)^(Nii - Sii) }
    # calculating probability densities
    probs.i <- pBLi(sek)
    probs.ii <- pBLii(sek)
  }
  
  # calculate vertical range
  ylim.max <- max(c(probs.i,probs.ii)) * fac
  
  ### plot graphical comparison p(H[delta]|S_1/N_1,I) versus p(H[delta]|S_2/N_2,I)
  if(pdfout) pdf(fname,width=9,height=6,paper="A4r")
  par(mar=c(5,6,5,5))
  plot(sek, probs.i,
       xlab="", ylab="probabilities",
       main="",
       type="l", lty=1, lwd=1.75, col="red", bty="l",
       ylim=c(0,ylim.max))
  points(sek, probs.ii, type="l", lty="dashed", lwd=1.75, col="blue")
  mtext(expression(paste(delta," Mean Difference",sep="")), 3, line=2, cex=1.5)
  mtext(eval(substitute(expression(paste("qualitative comparison of ",Si,"/",Ni," (red) vs. ",Sii,"/",Nii," (blue)")),
                        list(Si=Si,Ni=Ni,Sii=Sii,Nii=Nii))), 3, line=0.8)
  mtext(expression(paste(delta," = ", mu["1"]," - ",mu["2"])), 1, line=3, cex=1.5)
  mtext(
    eval(substitute(expression(paste(bar(x)["1"] ," = ",Di," | ",
                                     bar(x)["2"] ," = ",Dii," | ",
                                     sigma["1"] ," = ",si, " | ",
                                     sigma["2"] ," = ",sii, "")),list(Dii=round(Dii,dig),Di=round(Di,dig),si=round(si,dig),sii=round(sii,dig)))),
    4,line=1, cex=0.9)
  
  # add a nice legend to explain curves and associated probability densities
  legend("topleft", legend=c(expression(paste("p(H[",x["1"],"] | S" ["1"] ,"/ N" ["1"],", I)")),
                             expression(paste("p(H[",x["2"],"] | S" ["2"] ,"/ N" ["2"],", I)"))),
         text.col=c("red","blue"), bty="n", cex=1, y.intersp=1.4,
         title="probability densities", title.col="black", title.adj=1.1)
  if(pdfout) dev.off()
  # alternative: add another nice legend with information
  #legend("topright", legend=c(eval(substitute(expression(paste(bar(delta)["1"] ," = ",Di," (mean)")),list(Di=round(Di,dig)))),
  #                            eval(substitute(expression(paste(bar(delta)["2"] ," = ",Dii," (mean)")),list(Dii=round(Dii,dig)))),
  #                            eval(substitute(expression(paste(sigma["1"] ," = ",si, " (SD)")),list(si=round(si,dig)))),
  #                            eval(substitute(expression(paste(sigma["2"] ," = ",sii, " (SD)")),list(sii=round(sii,dig))))),
  #                            text.col=c("black"),
  #                            bty="n", cex=1,
  #                            y.intersp=1.4,
  #                            title="values", title.col="black")
  
}
################################################################################
#call
#UMSplot(inputv=res.SIB)
#
#call with pdf output
#UMSplot(inputv=res.SIB,pdfout=TRUE,fname="UMSplot_1998_p48.pdf")

# --------------------------------------------------
# End: ON THE DIFFERENCE IN MEANS
# --------------------------------------------------

