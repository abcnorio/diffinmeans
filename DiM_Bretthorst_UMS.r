### (C) 2005-2023 by Leo Guertler 
### R-code supplement
### to the book
###
### "Subjektive Ansichten und objektive Betrachtungen"
###
### written by Gürtler & Huber (2023)
###
### All R-code is published under the GPL v3 license:
###
### https://www.gnu.org/licenses/gpl-3.0.en.html
###
### except for 'borrowed' code - see links and references.
### For this R-code the original license of the respective
### authors is valid.
###
### R-code published on
###
### https://osdn.net/projects/mixedmethod-rcode
### https://github.com/abcnorio/mixedmethod-rcode


################################################################################
# --------------------------------------------------
# ON THE DIFFERENCE IN MEANS
# --------------------------------------------------


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
# call:
# res.SIB <- SucRatesIntBounds(Si=11, Ni=15, Sii=10, Nii=16, smin=0, snames=c("voluntary","non-voluntary"))
################################################################################ 


#########################################################################
# DIFFERENCE IN MEANS INTEGRALS
#
### PMV - NORMAL
#
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
# call:
# PMV <- PMV.hypo(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH)


### PMV - LOG
PMV.hypo.log <- function(NN,dd,upnum,lownum,sL,sH, Nsteps=100, ...)
{
  half <- 1/2
  integpmv.log <- function(s, ...)
  { 
    #      1 / (s^NN) * exp(-dd / (2 * s^2)) *
    #        ( pgamma(upnum/(2*s^2),1/2)*gamma(1/2) + pgamma(lownum/(2*s^2),1/2)*gamma(1/2) )
    sqrtpi.log <- lgamma(1/2)
    #      logs <- c(pgamma(upnum/(2*s^2),1/2,log=TRUE) + sqrtpi.log,
    #                pgamma(lownum/(2*s^2),1/2,log=TRUE) + sqrtpi.log
    #                )
    #pmv <- integrate(integpmv, lower=sL, upper=sH)$value
    #PMV <- pmv / sqrt(2*NN)
    #res <- log(1) - log(s)*NN + (-dd / (2*s^2)) + logSumExp(logs)
    s2.2 <- 2*s^2
#       v <- c( pgamma(upnum/(2*s^2),1/2,log=TRUE,lower=TRUE) + sqrtpi.log,
#              pgamma(lownum/(2*s^2),1/2,log=TRUE,lower=TRUE) + sqrtpi.log
     
    v <- c( pgamma(upnum/s2.2,half,log=TRUE,lower=TRUE) + sqrtpi.log,
            pgamma(lownum/s2.2,half,log=TRUE,lower=TRUE) + sqrtpi.log
          )
    #v
    -log(s)*NN + (-dd / (2*s^2)) + .llog.2add.short(v)
    #.llog.2add.short(pgamma(upnum/(2*s^2),1/2,log=TRUE,lower=TRUE) + sqrtpi.log,
    #               pgamma(lownum/(2*s^2),1/2,log=TRUE,lower=TRUE) + sqrtpi.log)
    #logSumExp(logs)
    }
    pmv.log <- simpsonrule.nlb(integpmv.log, sL, sH,
                               type="log",
                               Nsteps=Nsteps,
                               half=half)
    PMV.log <- pmv.log - log( sqrt(2*NN) )
return(PMV.log)
}
# call:
# PMV.log <- PMV.hypo.log(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH)
  

### PMV - BROB
PMV.hypo.brob <- function(NN,dd,upnum,lownum,sL,sH, Nsteps=100, ...)
{
  integpmv.brob <- function(s, ...)
  {
    1 / (as.brob(s)^NN) * exp(as.brob(-dd) / (2 * s^2)) *
      ( gamma_inc_gen(a=1/2,z0=0,z1=upnum/(2*s^2),log=TRUE)*sqrt(pi) +
        gamma_inc_gen(a=1/2,z0=0,z1=lownum/(2*s^2),log=TRUE)*sqrt(pi) )
  }
  pmv.brob <- simpsonrule.nlb(integpmv.brob, sL, sH, type="brob",Nsteps=Nsteps)
  PMV.brob <- pmv.brob / sqrt(2*NN)
return(PMV.brob)
}
# call:
# PMV.brob <- PMV.hypo.brob(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH, Nsteps=100)

# debug calls:
# PMV <- PMV.hypo(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH)
# PMV.log <- PMV.hypo.log(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH)
# PMV.brob <- PMV.hypo.brob(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH)
#
# log(PMV)
# PMV.log
# PMV.brob@x
    
  
### PMbarV - NORMAL
#  
# Bretthorst, p.6f.
# sign in inner brackets between the two pgamma*gamma's
# sign + in case if U1H, H1L are of different sign
# sign - in case if U1H, U1L are of same sign
# same for U2H, U2L
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
# call:
# PMbarV <- PMbarV.hypo(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)

  
### PMbarV - LOG  
PMbarV.hypo.log <- function(NN,zz,dd,upinum,lowinum,sL,sH,H,L,Ni,Nii, Nsteps=100, ...)
{
  sqrtpi.log <- lgamma(1/2)
  half <- 1/2
  integpmbarv.log <- function(s, ...)
  {
    #(* + if (L-DD) < 0 *) 
    #1/(s^(NN-1)) * exp(-zz/(2*s^2)) *
    #  ( pgamma(upinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowinum/(2*s^2),1/2)*gamma(1/2) ) *
    #  ( pgamma(upiinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowiinum/(2*s^2),1/2)* gamma(1/2) )
    s2.2 <- 2*s^2

#      v1 <- c( pgamma(upinum/(2*s^2),1/2, log=TRUE) + sqrtpi.log, 
#               pgamma(lowinum/(2*s^2),1/2, log=TRUE) + sqrtpi.log )
#      v2 <- c( pgamma(upiinum/(2*s^2),1/2, log=TRUE) + sqrtpi.log, 
#               pgamma(lowiinum/(2*s^2),1/2, log=TRUE) + sqrtpi.log )
 
    v1 <- c( pgamma(upinum/s2.2,half, log=TRUE) + sqrtpi.log, 
             pgamma(lowinum/s2.2,half, log=TRUE) + sqrtpi.log )
    v2 <- c( pgamma(upiinum/s2.2,half, log=TRUE) + sqrtpi.log, 
             pgamma(lowiinum/s2.2,half, log=TRUE) + sqrtpi.log )

    # log(1) == 0
    log(1) - ( (NN-1)*log(s) )  +
             (-zz/(2*s^2)) +
             .llog.2add.short( v1 ) +
             .llog.2add.short( v2 )
  }
  #pmbarv <- integrate(integpmbarv, lower=sL, upper=sH)$value
  pmbarv.log <- simpsonrule.nlb(integpmbarv.log, sL, sH,
                                sqrtpi.log=sqrtpi.log,
                                half=half,
                                type="log", Nsteps=Nsteps)
  #PMbarV <- pmbarv / (2*(H-L) * sqrt(Ni*Nii))
  PMbarV.log <- pmbarv.log - ( log(2*(H-L)) +
                (1/2)*( log(Ni) + log(Nii)) ) #log(sqrt(Ni*Nii)) )
return(PMbarV.log)
}
# call:
# PMbarV <- PMbarV.hypo.loog(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)

  
### PMbarV - BROB
PMbarV.hypo.brob <- function(NN,zz,dd,upinum,lowinum,sL,sH,H,L,Ni,Nii, Nsteps=100, ...)
{
  integpmbarv.brob <- function(s)
  {
    1/(as.brob(s)^(NN-1)) * exp(as.brob(-zz)/(2*s^2)) *
      ( pgamma(upinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowinum/(2*s^2),1/2)*gamma(1/2) ) *
      ( pgamma(upiinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowiinum/(2*s^2),1/2)*gamma(1/2) )
  }
  #    pmbarv.brob <- simpsonrule.brob(fx=integpmbarv.brob, sL=sL, sH=sH, Nsteps=Nsteps)
  pmbarv.brob <- simpsonrule.nlb(fx=integpmbarv.brob, sL, sH, type="brob",Nsteps=Nsteps)
  PMbarV.brob <- pmbarv.brob / (2*(H-L) * sqrt(Ni*Nii))
return(PMbarV.brob)  
}  
# call:
# PMbarV.brob <- PMbarV.hypo.brob(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=100)
  
# debug calls:  
# PMbarV <- PMbarV.hypo(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
# PMbarV.log <- PMbarV.hypo.log(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
# PMbarV.brob <- PMbarV.hypo.brob(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
#
# log(PMbarV)
# PMbarV.log
# PMbarV.brob
  
  
### PMVbar - NORMAL
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
# call:
# PMVbar <- PMVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
  

### PMVbar - LOG
#UiA <- function(A) Ni*(Dsi-2*Di*A+A^2)/2
#UiiA <- function(A) Nii*(Dsii-2*Dii*A+A^2)/2
UiA.log <- function(A) log(Ni) + log(Dsi-2*Di*A+A^2) -log(2) #ok
UiiA.log <- function(A) log(Nii) + log(Dsii-2*Dii*A+A^2) -log(2) #ok
  
PMVbar.hypo.log <- function(Dsi,Dsii,sL,sH,H,L,Ni,Nii, Nsteps=100, ...)
{
  #const.
  lg.Ni.2 <- lgamma(Ni/2)
  lg.Nii.2 <- lgamma(Nii/2)
  Ni.2 <- Ni/2
  Nii.2 <- Nii/2
  sH.2 <- sH^2
  sL.2 <- sL^2
  g.Ni.2 <- gamma(Ni/2)
  g.Nii.2 <- gamma(Nii/2)
  
  integpmvbar.log <- function(A, ...)
  {
    #1/UiA(A)^(Ni/2) *
    #  1/UiiA(A)^(Nii/2) *
    #  ( pgamma(UiA(A)/(sH^2),Ni/2)*gamma(Ni/2) - pgamma(UiA(A)/(sL^2),Ni/2)*gamma(Ni/2) ) *
    #  ( pgamma(UiiA(A)/(sH^2),Nii/2)*gamma(Nii/2) - pgamma(UiiA(A)/(sL^2),Nii/2)*gamma(Nii/2) )
    
    UiA.A <- UiA(A)
    UiiA.A <- UiiA(A)
    #Ni.2 <- Ni/2
    #Nii.2 <- Nii/2
    v1 <- c( pgamma( UiA.A/sH.2, Ni.2, log=TRUE,lower=FALSE ) + lg.Ni.2, #a1
             pgamma( UiA.A/sL.2, Ni.2, log=TRUE,lower=FALSE ) + lg.Ni.2
           )
    v2 <- c( pgamma( UiiA.A/sH.2, Nii.2, log=TRUE,lower=FALSE ) + lg.Nii.2, #a2
             pgamma( UiiA.A/sL.2, Nii.2, log=TRUE,lower=FALSE ) + lg.Nii.2
           )
    log(1) - ( Ni.2 * log(UiA.A) ) + #ok
    log(1) - ( Nii.2 * log(UiiA.A) ) + #ok
    .llog.2sub.short(v1) + .llog.2sub.short(v2)
  }
    
  pmvbar.log <- simpsonrule.nlb(integpmvbar.log,
                                L, H,
                                sL.2=sL.2, sH.2=sH.2,
                                lg.Ni.2=lg.Ni.2, lg.Nii.2=lg.Nii.2,
                                Ni.2=Ni.2, Nii.2=Nii.2,
                                g.Ni.2=g.Ni.2,
                                g.Ni.2=g.Nii.2,
                                type="log", Nsteps=Nsteps)
  PMVbar.log <- pmvbar.log - ( log(4) + log( log(sH/sL) ) )
return(PMVbar.log)  
}  
# call:
# PMVbar.log <- PMVbar.hypo.log(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
  

### PMVbar - BROB
# not possible for pgamma due to its definition without brob numbers  
UiA.brob <- function(A) as.brob( Ni*(Dsi-2*Di*A+A^2)/2 )
UiiA.brob <- function(A) as.brob( Nii*(Dsii-2*Dii*A+A^2)/2 )
  
PMVbar.hypo.brob <- function(Dsi,Dsii,sL,sH,H,L,Ni,Nii, Nsteps=100, ...)
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
  #pmvbar.brob <- simpsonrule.brob(fx=integpmvbar.brob, sL=L, sH=H, Nsteps=Nsteps)
  pmvbar.brob <- simpsonrule.nlb(fx=integpmvbar.brob, L, H, type="brob", Nsteps=Nsteps)
  PMVbar.brob <- pmvbar.brob / (4*log(sH/sL))    
return(PMVbar.brob)
}
# call:
# PMVbar.brob <- PMVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=100)
  
# debug calls:
# PMVbar <- PMVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
# PMVbar.log <- PMVbar.hypo.log(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
# PMVbar.brob <- PMVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
#
# log(PMVbar)
# PMVbar.log
# PMVbar.brob

    
### PMbarVbar - NORMAL
# Bretthorst, p.11
# sign: minus if U1H, U1L have same sign
# sign: plus if U1H, U1L have different sign
# same true for U2H, U2L
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
# call:
# PMbarVbar <- PMbarVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
  

### PMbarVbar - LOG
# Bretthorst, p.11
# sign: minus if U1H, U1L have same sign
# sign: plus if U1H, U1L have different sign
# same true for U2H, U2L
  
#UiiB <- function(B) Nii*(Dsii-2*Dii*B+B^2)/2
#  UiiB.log <- function(B, ...)
#  {
#    #a <- log(Nii) - log(2)
#    #b <- log(Dsii)
#    #ce <- log(2) + log(Dii) + log(B)
#    #d <- 2*log(B)
#    v1 <- c(log(Dsii), 2*log(B))
#    v2 <- log(2) + log(Dii) + log(B)
#    #return( a + .llog.2sub.short( c( .llog.2add.short(c(b,d)),ce) ) )
#    return( log(Nii) - log(2) + .llog.2sub.short( c( .llog.2add.short(v1),v2) ) )
#  }
UiiB.log <- function(B, ...)
{
  v1 <- c(log(Dsii), 2*log(B))
  v2 <- log(2) + log(Dii) + log(B)
return( log(Nii) - log(2) + .llog.2sub.short( c( .llog.2add.short(v1),v2) ) )
}
# call:
# UiiB.log(B)
# log(UiiB(B))
  
PMbarVbar.hypo.log <- function(Dsi,Dsii,sL,sH,H,L,Ni,Nii, Nsteps=100, ...)
{
  lg.Ni.2 <- lgamma(Ni/2)
  lg.Nii.2 <- lgamma(Nii/2)
  Ni.2 <- Ni/2
  Nii.2 <- Nii/2
  sH.2 <- sH^2
  sL.2 <- sL^2
    
  integpmbarvbar1.log <- function(A, ...)
  {
    #1/(UiA(A)^(Ni/2)) *
    #lower=FALSE for integral z^inf
    #the absolute value is ok, only the sign would differ
    #the multiplication of both integrals if handled identical
    #will lead to positive values
    #and with prob distributions we do not have negative probs!
    #  ( pgamma(UiA(A)/(sH^2),Ni/2,lower=FALSE)*gamma(Ni/2) -
    #      pgamma(UiA(A)/(sL^2),Ni/2,lower=FALSE)*gamma(Ni/2) )
      
    #v1 <- c( pgamma(UiA(A)/sH.2, Ni.2, lower=FALSE, log=TRUE) + lg.Ni.2,
    #         pgamma(UiA(A)/sL.2, Ni.2, lower=FALSE, log=TRUE) + lg.Ni.2 )
    UiA.A <- UiA(A)
    v1 <- c( pgamma(UiA.A/sH.2, Ni.2, lower=FALSE, log=TRUE) + lg.Ni.2,
             pgamma(UiA.A/sL.2, Ni.2, lower=FALSE, log=TRUE) + lg.Ni.2 )
    #log(1) - ( Ni.2 * log(UiA(A)) ) + .llog.2sub.short(v1) #todo check which one is greater! -> put into function!!!!!!!!!!!!!
    log(1) - ( Ni.2 * log(UiA.A) ) + .llog.2sub.short(v1) #todo check which one is greater! -> put into function!!!!!!!!!!!!!
  }
    
  integpmbarvbar2.log <- function(B, ...)
  {
    #ok#1/(UiiB(B)^(Nii/2)) *
    #  #lower=FALSE for integral z^inf
    #ok#  ( pgamma(UiiB(B)/(sH^2),Nii/2,lower=FALSE)*gamma(Nii/2) -
    #ok#      pgamma(UiiB(B)/(sL^2),Nii/2,lower=FALSE)*gamma(Nii/2) )
      
    UiiB.B <- UiiB(B)
 #     v2 <- c( pgamma(UiiB(B)/sH.2, Nii.2, lower=FALSE, log=TRUE) + lg.Nii.2,
 #              pgamma(UiiB(B)/sL.2, Nii.2, lower=FALSE, log=TRUE) + lg.Nii.2 )
       
    v2 <- c( pgamma(UiiB.B/sH.2, Nii.2, lower=FALSE, log=TRUE) + lg.Nii.2,
             pgamma(UiiB.B/sL.2, Nii.2, lower=FALSE, log=TRUE) + lg.Nii.2 )
    log(1) - ( Nii.2 * log(UiiB(B)) ) + .llog.2sub.short(v2)
  }
    
  pmbarvbar1.log <- simpsonrule.nlb(integpmbarvbar1.log,
                                    lower=L, upper=H,
                                    lg.Ni.2=lg.Ni.2,
                                    Ni.2=Ni.2,
                                    sL.2=sL.2, sH.2=sH.2,
                                    sL=sL, sH=sH,
                                    type="log", Nsteps=Nsteps)
  pmbarvbar2.log <- simpsonrule.nlb(integpmbarvbar2.log,
                                    lower=L, upper=H,
                                    lg.Ni.2=lg.Ni.2,
                                    Ni.2=Ni.2,
                                    sL.2=sL.2, sH.2=sH.2,
                                    sL=sL, sH=sH,
                                    type="log", Nsteps=Nsteps)
    
  #pmbarvbar <- pmbarvbar1 * pmbarvbar2
  #PMbarVbar <- pmbarvbar / (4*(H-L) * log(sH/sL))
  v1 <- c(log(H), log(L))
  v2 <- log(log(sH)-log(sL))
  pmbarvbar.log <- pmbarvbar1.log + pmbarvbar2.log
  PMbarVbar.log <- pmbarvbar.log - 
                   ( log(4) + .llog.2sub.short(v1) + v2 )
    
return(PMbarVbar.log)
}  
# call:
# PMbarVbar.log <- PMbarVbar.hypo.log(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
  

### PMbarVbar - BROB
UiiB.brob <- function(B) as.brob( Nii*(Dsii-2*Dii*B+B^2)/2 )
  
PMbarVbar.hypo.brob <- function(Dsi,Dsii,sL,sH,H,L,Ni,Nii, Nsteps=100, ...)
{
  #first integral
  #lower=FALSE for new simsponrule.nlb!!!
  integpmbarvbar1.brob <- function(A)
  {
    1/( UiA(as.brob(A))^as.brob(Ni/2) )  *
      ( brob( pgamma(UiA(A)/(sH^2),Ni/2,log.p=TRUE,lower=FALSE)+lgamma(Ni/2) ) -
        brob( pgamma(UiA(A)/(sL^2),Ni/2,log.p=TRUE,lower=FALSE)+lgamma(Ni/2) ) )  
  }  
  #second integral
  #lower=FALSE for new simsponrule.nlb!!!
  integpmbarvbar2.brob <- function(B)
  {
    1/( UiiB.brob(B)^(Nii/2)) *
      ( brob(pgamma(UiiB(B)/(sH^2),Nii/2,log.p=TRUE,lower=FALSE)+lgamma(Nii/2)) -
        brob(pgamma(UiiB(B)/(sL^2),Nii/2,log.p=TRUE,lower=FALSE)+lgamma(Nii/2)) )
  }
  #pmbarvbar1.brob <- simpsonrule.brob(fx=integpmbarvbar1.brob, sL=L, sH=H, Nsteps=Nsteps)
  #pmbarvbar2.brob <- simpsonrule.brob(fx=integpmbarvbar2.brob, sL=L, sH=H, Nsteps=Nsteps)
    
  pmbarvbar1.brob <- simpsonrule.nlb(fx=integpmbarvbar1.brob, L, H, type="brob", Nsteps=Nsteps)
  pmbarvbar2.brob <- simpsonrule.nlb(fx=integpmbarvbar2.brob, L, H, type="brob", Nsteps=Nsteps)
    
  pmbarvbar.brob <- pmbarvbar1.brob * pmbarvbar2.brob
  PMbarVbar.brob <- pmbarvbar.brob / (4*(H-L) * log(sH/sL))
return(PMbarVbar.brob)
}
# call:
# PMbarVbar.brob <- PMbarVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=100)
  
# debug calls:
# PMbarVbar <- PMbarVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
# PMbarVbar.log <- PMbarVbar.hypo.log(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
# PMbarVbar.brob <- PMbarVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
#
# log(PMbarVbar)
# PMbarVbar.log
# PMbarVbar.brob
  
################################################################################

  
################################################################################
# create constants from initial values
.DiM_ccs <- function(inval=NULL, smin=0)
{
  # create constsns
  # extract values for calculations from inputv object (class "UMS")
  
  # general
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
  
  # PMV
  dd <- NN * (DsD - DD^2)
  lownum <- NN * (L-DD)^2
  upnum  <- NN * (H-DD)^2
  
  # PMbarV
  zz <- Ni * (Dsi - Di^2) + Nii * (Dsii - Dii^2)
  lowinum <- Ni * (L - Di)^2
  upinum <- Ni * (H - Di)^2
  lowiinum <- Nii * (L - Dii)^2
  upiinum <- Nii * (H - Dii)^2

  # compile input results independent from method
  constants <- data.frame(sample1=snames[1],sample2=snames[2],
                          Di,Dii,Ni,si,Nii,sii,
                          NN,DD,ss,L,H,sL,sH,Dsi,Dsii,
                          DsD,dd,
                          lownum,upnum,
                          lowinum,upinum,
                          lowiinum,upiinum,
                          smin,
                          zz)
  # resulting table/ dataframe with input values
  input.df <- data.frame("No." = c(Ni,Nii,NN),
                         "Standard Deviation" = c(si,sii,ss),
                         "Mean" = c(Di,Dii,DD),
                         "Data set" = c(snames,"combined"),
                         check.names = FALSE)
  # table/ dataframe with prior values/ information
  prior.df <- data.frame("Numerical Example" = c("Prior Mean lower bound","Prior Mean upper bound",
                                                 "Prior Standard Deviation lower bound","Prior Standard Deviation upper bound"),
                         "Value" = c(L,H,sL,sH), check.names = FALSE) 
  
  res <- list(constants=constants,
              input.df=input.df,
              prior.df=prior.df
              )
return(res)
}
# call:
# 
################################################################################


################################################################################
# calculate the four integrals
# PMV, PMbarV, PMVbar, PMbarVbar
.DiM_callintegrals <- function(DiM.ccs,
                               TYPE="normal",
                               Nsteps=100,
                               parallel=FALSE
                              )
{
  
  # general note
  # necessary variables in the function
  # NN, DD, Dsi, Dsii, DsD, ss,
  # dd, lownum, upnum, low, up, pmv, PMV,
  # zz, lowinum, upinum, lowiinum, upiinum, pmbarv, PMbarV,
  # pmvbar, PMVbar, pmbarvbar, PMbarVbar, cc,
  # mv, mbarv, mvbar, mbarvbar,
  # samemeans, diffmeans, samevars, diffvars, diffsets

  constants <- DiM.ccs[["constants"]]
  attach(constants)
    
  # notes    
  # Bretthorst, p.6f.
  # sign in inner brackets between the two pgamma*gamma's
  # sign + in case if U1H, H1L are of different sign
  # sign - in case if U1H, U1L are of same sign
  # same for U2H, U2L
  
  # PMVbar
  # no constants
  
  # PMbarVbar
  # no constants

  # Bretthorst, p.11
  # sign: minus if U1H, U1L have same sign
  # sign: plus if U1H, U1L have different sign
  # same true for U2H, U2L
  
  # integpmbarvbar1
  # lower=FALSE for integral z^inf
  # the absolute value is ok, only the sign would differ
  # the multiplication of both integrals if handled identical
  # will lead to positive values
  # and with prob distributions we do not have negative probs!
  #
  # brob version:
  # first and second integral with
  # lower=FALSE for new simpsonrule.nlb!!!  

  # calculate integrals  
  if(TYPE=="normal")
  {
    cat("\nWork based on normal values\n")
    cat("\nCalculate PMV\n")
    PMV <- PMV.hypo(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH)
    cat("Calculate PMbarV\n")
    PMbarV <- PMbarV.hypo(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
    cat("Calculate PMVbar\n")
    PMVbar <- PMVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)
    cat("Calculate PMbarVbar\n")
    PMbarVbar <- PMbarVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii)

  } else if(TYPE=="log")
  {
    cat("\nWork based on log values\n")
    if(parallel) require(parallel)
    
    cat("\nCalculate PMV\n")
    PMV <- PMV.hypo.log(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH, Nsteps=Nsteps)
    cat("Calculate PMbarV\n")
    PMbarV <- PMbarV.hypo.log(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)
    cat("Calculate PMVbar\n")
    PMVbar <- PMVbar.hypo.log(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)
    cat("Calculate PMbarVbar\n")
    PMbarVbar <- PMbarVbar.hypo.log(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)

  } else if(TYPE=="brob")
  {
    cat("\nWork based on Brobdingnag values\n")
    require(Brobdingnag)

    cat("\nCalculate PMV\n")
    PMV <- PMV.hypo.brob(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH, Nsteps=Nsteps)
    cat("Calculate PMbarV\n")
    PMbarV <- PMbarV.hypo.brob(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)
    cat("Calculate PMVbar\n")
    PMVbar <- PMVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)
    cat("Calculate PMbarVbar\n")
    PMbarVbar <- PMbarVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)

  }
  detach(constants)
  
  DiM.integrals <- c(PMV, PMbarV, PMVbar, PMbarVbar)
  names(DiM.integrals) <- c("PMV","PMbarV","PMVbar","PMbarVbar")

return(DiM.integrals)
}
# call:
# 
################################################################################

  
################################################################################
# compile results bsed on the four integrals
.DiM_prepres <- function(DiM.integrals,
                         TYPE=NA,
                         convback=FALSE, # convert probs back from log-scale
                         percfac=1, # not 100% but sum(ps) = 1
                         dig=5
                        )
{
  
  stopifnot(!is.na(TYPE))
  
  if(TYPE=="normal" && convback==TRUE)
  {
    cat("\nTYPE=normal & convback=TRUE -> no conversion\nwill be done for values on the normal-scale.\n\n")
    convback <- FALSE
  }
  
  if(TYPE=="normal")
  {
    denom <- sum(DiM.integrals) # PMV, PMbarV, PMVbar, PMbarVbar
    cc <- 1 / denom
    
    # output not as p*100%, just p*1
    mv <- percfac * cc * DiM.integrals["PMV"]
    mbarv <- percfac * cc * DiM.integrals["PMbarV"]
    mvbar <- percfac * cc * DiM.integrals["PMVbar"]
    mbarvbar <- percfac * cc * DiM.integrals["PMbarVbar"]
    
    samemeans <- mv + mvbar
    diffmeans <- mbarv + mbarvbar
    samevars <- mv + mbarv
    diffvars <- mvbar + mbarvbar
    diffsets <- mvbar + mbarv + mbarvbar
    samesets <- 1 - diffsets
    
    OR.diffmeans <- diffmeans/samemeans
    OR.samemeans <- samemeans/diffmeans
    OR.diffvars <- diffvars/samevars
    OR.samevars <- samevars/diffvars
    OR.diffsets <- diffsets/mv
    OR.samesets <- mv/diffsets

  } else if(TYPE=="log")
  {
    denom <- .llog.2add.short( DiM.integrals )
    #.llog.2add.short(c(PMV, PMbarV, PMVbar, PMbarVbar))
    cc <- - denom #log(1)
    
    # no *100%, just *1
    percfac.l <- log(percfac)
    mv <- percfac.l + cc + DiM.integrals["PMV"]
    mbarv <- percfac.l + cc + DiM.integrals["PMbarV"]
    mvbar <- percfac.l + cc + DiM.integrals["PMVbar"]
    mbarvbar <- percfac.l + cc + DiM.integrals["PMbarVbar"]
    
    samemeans <- .llog.2add.short( c(mv, mvbar) )
    diffmeans <- .llog.2add.short( c(mbarv, mbarvbar) )
    samevars <- .llog.2add.short( c(mv, mbarv) )
    diffvars <- .llog.2add.short( c(mvbar, mbarvbar) )
    diffsets <- .llog.2add.short( c(mvbar, mbarv, mbarvbar) )
    
    samesets <- .llog.2sub.short( c(log(1),diffsets) )
    
    OR.diffmeans <- diffmeans - samemeans
    OR.samemeans <- samemeans - diffmeans
    OR.diffvars <- diffvars - samevars
    OR.samevars <- samevars - diffvars
    OR.diffsets <- diffsets - mv
    OR.samesets <- mv - diffsets
  
  } else if(TYPE=="brob")
  {
    # calculate total probability (denominator of Bayes' Theorem)
    # sum of all four hypotheses (probabilities)
    denom <- sum(list2vec.brob(DiM.integrals))
    cc <- 1 / denom
    
    # output not as p*100%, just p*1
    mv <- percfac * cc * DiM.integrals[["PMV"]]
    mbarv <- percfac * cc * DiM.integrals[["PMbarV"]]
    mvbar <- percfac * cc * DiM.integrals[["PMVbar"]]
    mbarvbar <- percfac * cc * DiM.integrals[["PMbarVbar"]]
    
    samemeans <- mv + mvbar
    diffmeans <- mbarv + mbarvbar
    samevars <- mv + mbarv
    diffvars <- mvbar + mbarvbar
    diffsets <- mvbar + mbarv + mbarvbar
    samesets <- 1 - diffsets
    
    OR.diffmeans <- diffmeans/samemeans
    OR.samemeans <- samemeans/diffmeans
    OR.diffvars <- diffvars/samevars
    OR.samevars <- samevars/diffvars
    OR.diffsets <- diffsets/mv
    OR.samesets <- mv/diffsets
  
  } 
  
  # unnormalized + total prob
  DiM.integrals.ext <- c(DiM.integrals, cc, denom)
  names(DiM.integrals.ext)[c(5,6)] <- c("cc","denom")
  
  # normalized + Odds Ratios (OR)
  DiM.probsOR.res <- c(mv, mbarv, mvbar, mbarvbar,
                       samemeans, diffmeans, samevars, diffvars,
                       samesets, diffsets,
                       OR.diffmeans, OR.samemeans,
                       OR.diffvars, OR.samevars,
                       OR.diffsets, OR.samesets
                      )
  names(DiM.probsOR.res) <- c("mv", "mbarv", "mvbar", "mbarvbar",
                              "samemeans", "diffmeans", "samevars", "diffvars",
                              "samesets", "diffsets",
                              "OR.diffmeans", "OR.samemeans",
                              "OR.diffvars", "OR.samevars",
                              "OR.diffsets", "OR.samesets"
                             )
  # probs = 1:10
  # ORs = 11:16

  # diffmeans vs samemeans
  # diffvars vs samevars
  # diffsets vs samesets
  diff.v <- c("diffmeans","diffvars","diffsets")
  same.v <- c("samemeans","samevars","mv")
  diff.terms <- c(rep("different",2),"different means and/ or")
  same.terms <- c(rep("the same",2), "the same means and")
  diffvssame.v <- c("diffvssamemeans","diffvssamevars","diffvssamesets")
  samevsdiff.v <- c("samevsdiffmeans","samevsdiffvars","samevsdiffsets")
  diffsame.terms <- diff.terms
  names(diffsame.terms) <- diff.v

  # diff greater/ same for means/vars/sets  
  if(TYPE=="normal")
  {
    diffgreaterthansame <- DiM.probsOR.res[ diff.v ] > DiM.probsOR.res[ same.v ]
    diffbysame <- DiM.probsOR.res[ diff.v ] / DiM.probsOR.res[ same.v ]
    diffbysame.ID <- which(diffgreaterthansame == FALSE)
    if(length(diffbysame.ID) > 0) diffbysame[ diffbysame.ID ] <- 1 / diffbysame[ diffbysame.ID ]

  } else if(TYPE=="log")
  {
    diffgreaterthansame <- DiM.probsOR.res[ diff.v ] > DiM.probsOR.res[ same.v ]
    diffbysame <- DiM.probsOR.res[ diff.v ] - DiM.probsOR.res[ same.v ]
    diffbysame.ID <- which(diffgreaterthansame == FALSE)
    if(length(diffbysame.ID) > 0) diffbysame[ diffbysame.ID ] <- log(1) - diffbysame[ diffbysame.ID ]
    
  } else if(TYPE=="brob")
  {

    diffbysame <- sapply(seq_along(1:length(diff.v)), function(i)
          {
            DiM.probsOR.res[diff.v[i]][[1]] / DiM.probsOR.res[same.v[i]][[1]]
          })
    names(diffbysame) <- diff.v
    
    diffgreaterthansame <- sapply(seq_along(1:length(diff.v)), function(i)
          {
            DiM.probsOR.res[diff.v[i]][[1]] > DiM.probsOR.res[same.v[i]][[1]]
          })
    names(diffgreaterthansame) <- diff.v
    
    diffbysame.ID <- which(diffgreaterthansame == FALSE)
    if(length(diffbysame.ID) > 0)
    {
      diffbysame[ diffbysame.ID ][[1]] <- 1 / diffbysame[ diffbysame.ID ][[1]]
    }

  }
  
  diff.t <- c("Different vs. same means","Different vs. same variances","Different vs. same sets")
  same.t <- c("Same vs. different means","Same vs. different variances","Same vs. different sets")

  if(length(diffbysame.ID) > 0) 
  {
    diffsame.terms[ diffbysame.ID ] <- same.terms[ diffbysame.ID ]
    diffvssame.v[ diffbysame.ID ] <- samevsdiff.v [ diffbysame.ID ]
    diff.t[ diffbysame.ID ] <- same.t[ diffbysame.ID ]
  }
  
  DiM.probsOR.orig <- c(DiM.integrals.ext, DiM.probsOR.res, diffbysame)
  
  # convert back to normal scale (log, brob) - but not case normal
  if(convback==FALSE || TYPE=="normal")
  {
    cat("\nMethod = ",TYPE,"\nconvback = ",convback,"\n\nProbabilities and Odds Ratios (= OR) will not be\nconverted from log-scale back to non-log-scale.\n\n",sep="")
    DiM.probsOR.exp <- NA
  } else
  {
    if(TYPE=="brob")
    {
      DiM.probsOR.exp <- sapply(DiM.probsOR.orig, as.numeric)
    } else if(TYPE=="log")
    {
      DiM.probsOR.exp <- exp(DiM.probsOR.orig)
    }
  }

  # correct wrong names
  names(diffgreaterthansame) <- names(DiM.probsOR.orig)[23:25] <- diffvssame.v
  if(!is.na(DiM.probsOR.exp[1])) names(DiM.probsOR.exp) <- names(DiM.probsOR.orig)

  # table/ dataframe with resulting probabilities
  if(TYPE=="brob")
  {
    p.H.D1D2I <- brobgiveback(DiM.probsOR.res[1:10])
    ORs <- sapply(DiM.probsOR.res[11:16], Brobdingnag:::getX)
    ORs <- brobgiveback(DiM.probsOR.res[11:16])
  } else if(TYPE %in% c("normal","log"))
  {
    p.H.D1D2I <- DiM.probsOR.res[1:10]
    ORs <- DiM.probsOR.res[11:16]
  } else
  {
    stop("Unknown TYPE=",TYPE)
  }

  # TRUE for log and brob
  cat("\nMethod = ",TYPE,"\nconvback = ",convback,"\n\n",sep="")
  tempnams <- c(rep("integrals (unnormalized)",4),
                c("cc","denom"),
                rep("integrals (normalized)",4),
                rep("hypothesis",6),
                rep("Odds Ratio",6),
                rep("Odds Ratio (diff vs same)",3)
                )
    
  # give out table to show conversion log/ brob back to normal
  if(TYPE=="brob") 
  {
    proutvalues <-brobgiveback(DiM.probsOR.orig)
  } else
  {
    proutvalues <- DiM.probsOR.orig
  }
    
  DiM.probs.res.df <- data.frame("Hypothesis (Probabilities)" = c(
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
                                 "p(H|D1,D2,I)"=p.H.D1D2I,
                                 check.names=FALSE
                                )
      
  # table/ dataframe with Odds Ratio results
  DiM.OR.res.df <- data.frame("Odds Ratio - Hypothesis in favor of ..." = c(
                              "A difference in Means",
                              "The same Means",
                              "A difference in Standard Deviations",
                              "The same Standard Deviations",
                              "A difference in the Sets (different Means and/ or Standard Deviations)",
                              "The same Sets (same Means and/ or Standard Deviations)"),
                              "Odds Ratio (OR)"=ORs,
                              check.names = FALSE
                             )

  if(TYPE == "brob") diffbysame <- brobgiveback(diffbysame)
  DiM.diffvssame.df <- data.frame("Odds Ratio - Hypothesis in favor of ..." =diff.t,
                                  "Odds Ratio (OR)"=diffbysame,
                                  check.names=FALSE
                                 )

  # TODO
  # infinity check

  # Inf test
  if(TYPE=="brob")
  {
    DiM.probs.inf.tf <- sapply(DiM.probsOR.exp, is.infinite) 
  } else # log or normal
  {
    DiM.probs.inf.tf <- DiM.probsOR.exp == abs(Inf)
  }
  
  if(is.logical(DiM.probs.inf.tf))
  {
    DiM.probs.inf.ID <- which(DiM.probs.inf.tf)
    if(length(DiM.probs.inf.ID) == 0) DiM.probs.inf.ID <- 0
  } else
  {
    DiM.probs.inf.ID <- NA
  }

  # results into a list
  res.short  <- list(DiM.probs.res.df=DiM.probs.res.df,
                     DiM.OR.res.df=DiM.OR.res.df,
                     DiM.diffvssame.df=DiM.diffvssame.df,
                     DiM.probsOR.exp=DiM.probsOR.exp
                    )
        
  res <- list(
              DiM.probsOR.orig=DiM.probsOR.orig,
              DiM.probsOR.exp=DiM.probsOR.exp,
              DiM.probs.inf.tf=DiM.probs.inf.tf,
              DiM.probs.inf.ID=DiM.probs.inf.ID,
              diffsame.terms=diffsame.terms,
              percfac=percfac,
              convback=convback,
              TYPE=TYPE,
              res.short=res.short #,
              #diffbysame=diffbysame,
              #diffbysame.exp=diffbysame.exp,
              #DiM.integrals.ext=DiM.integrals.ext,
              #DiM.probsOR.res=DiM.probsOR.res
             )

return(res)
}
# call:
# 
################################################################################


################################################################################
# meta function to calculate DiM
DiM <- function(initials=NA,
                cinput=NA,
                smin=0,
                cmethod="normal", # log, brob
                convback=FALSE,
                percfac=1,
                prout=TRUE,
                plotti=TRUE,
                legendplace="topleft",
                dig=3)
{
  
  if(cmethod == "normal" && convback == TRUE)
  {
    cat("\nmethod = 'normal' and convback == TRUE -> convback will be ignored.\n")
    convback <- FALSE
  }
  
  if(cmethod %in% c("log","brob") && is.na(convback))
  {
    stop("Choose convback = TRUE or FALSE, variable not set.")
  }
  
  if( !is.list(initials) || !(cinput %in% c("qual", "quan")) )
  {
    stop("Input either not list or unknown input method.")
  }

  if( any(is.na(initials)) || !any(sapply(initials, is.numeric)) )
  {
    stop("Either NAs or non-numeric values in 'initials'.")
  }
  
  initials.n <- names(initials)
  
  if("snames" %in% names(initials) &&
     length(initials$snames) > 0 &&
     is.character(initials$snames))
  { 
    snames <- initials$snames
  } else snames <- c("group1","group2")
  
  if(cinput == "qual")
  {
    # calculate sd + mean from success/ totals for each group
    cat("\nCalculate means and standard deviations from successes/ totals.\n")
    
    listofvars_qual <- c("Si","Ni","Sii","Nii")
    if( !all(listofvars_qual %in% initials.n) )
    {
      stop("Values are missing in 'initials' - we need:\n\n",paste(listofvars_qual,collapse="\n"),"\n",sep="")
    }

    inval <- with(initials,
      SucRatesIntBounds(Si=Si, Ni=Ni, Sii=Sii, Nii=Nii, smin=smin, snames=snames)
    )
  } else if(cinput == "quan")
  {
    # calculate based on sd + mean as well sa lower/ upper limits
    cat("\nCalculate based on means, standard deviations, and lower/ upper limits.\n")

    listofvars_quan <- c("Ni","Di","si","Nii","Dii","sii","L","H","sH","sL")
    if( !all(listofvars_quan %in% initials.n) )
    {
      stop("Values are missing in 'initials' - we need:\n\n",paste(listofvars_quan,collapse="\n"),"\n",sep="")
    }

    inval <- with(initials,
        inval <- list(Ni=Ni, Di=Di, si=si, Nii=Nii, Dii=Dii, sii=sii, L=L, H=H, sL=sL, sH=sH, snames=snames, smin=smin)
    )
    attr(inval, "typ") <- "QUAN"
  } else stop("Wrong oinput method specified.")
    
  # prepare constants
  cat("\nCalculate constants.\n")
  DiM.ccs <- .DiM_ccs(inval=inval)
  
  # calculate main integrals
  cat("\nCalculate main integrals:\n")
  DiM.integrals <- .DiM_callintegrals(DiM.ccs, TYPE=cmethod)
  
  # prepare output values and handle conversion back to normal
  # in case of different method etc
  cat("\nPrepare probabilities of the hypotheses.\n")
  DiM.out <- .DiM_prepres(DiM.integrals, TYPE=cmethod, convback=convback, percfac=percfac, dig=dig)
  
  # print nice output
  if(prout == "long")
  {
    cat("\nPreparing nice long output.\n")
    UMSprint(DiM.out, inval=inval, dig=dig)
  } else if(prout == "short")
  {
    cat("\nPreparing short and efficient output.\n")
    DiM.prout.short <- list(DiM.probs.res.df=DiM.out$res.short$DiM.probs.res.df,
                            DiM.OR.res.df=DiM.out$res.short$DiM.OR.res.df,
                            DiM.diffvssame.df=DiM.out$res.short$DiM.diffvssame.df,
                            DiM.probsOR.exp=DiM.out$res.short$DiM.probsOR.exp,
                            convback=convback)
    .DiM_prout.short(DiM.prout.short, DiM.ccs, dig=dig)
  } else cat("\nNo output.\n")
  
  # plot qual
  if(cinput == "qual" && plotti == TRUE)
  {
    cat("\nPlot curves for success/ failures.\n")
    UMSplot(inval=inval, legendplace=legendplace)
  }
  
  # prepare results
  DiM.res <-list(initials=initials,
                 inval=inval,
                 DiM.ccs=DiM.ccs,
                 DiM.integrals=DiM.integrals,
                 DiM.out=DiM.out
                 #DiM.prout.short=DiM.prout.short
                 )
  
return(DiM.res)
}
# call:
#
###############################################################################

  
###############################################################################
# plot graphical comparison
# only qualitative Bretthorst
UMSplot <- function(inval, dig=4, pdfout=FALSE, legendplace="topleft",
                    fname="UMSplot.pdf", loga=TRUE, fac=1.1, Nsteps=100)
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
  mtext(eval(substitute(expression(paste(bar(x)["1"] ," = ",Di," | ",
                                         bar(x)["2"] ," = ",Dii," | ",
                                         sigma["1"] ," = ",si, " | ",
                                         sigma["2"] ," = ",sii, "")),list(Dii=round(Dii,dig),Di=round(Di,dig),si=round(si,dig),sii=round(sii,dig)))),
                        4,line=1, cex=0.9)
  
  # add a nice legend to explain curves and associated probability densities
  legend(legendplace, legend=c(expression(paste("p(H[",x["1"],"] | S" ["1"] ,"/ N" ["1"],", I)")),
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
# call:
# UMSplot(inval=res.SIB)
#
# call with pdf output:
# UMSplot(inval=res.SIB,pdfout=TRUE,fname="UMSplot_1998_p48.pdf")
################################################################################


################################################################################
# print out DiM results - long version
UMSprint <- function(results, inval, dig=3)
{
  TYPE <- results[["TYPE"]]
  consts <- DiM.ccs[["constants"]]
  consts.out <- signif( consts[,-c(1,2)], digits=dig )
  convback <- results[["convback"]]
    
  if(attr(inval,"typ") == "SRIB")
  {
    srib.base <- TRUE
  } else
  {
    srib.base <- FALSE
  }

  if(TYPE %in% c("log","brob") && convback == TRUE)
  {
    res.out <- signif( results$DiM.probsOR.exp, digits=dig)
  } else if( TYPE == "log" && convback == FALSE)
  {
    res.out <- paste("exp(", signif( results$DiM.probsOR.orig, digits=dig), ")", sep="")
    names(res.out) <- names(results$DiM.probsOR.orig)
  } else if(TYPE == "brob" && convback == FALSE)
  {
    # as vector and character, we won't calculate anything anymore
    # just output
    if(dig < 5) dig <- 5
    res.out <- brobgiveback(results$DiM.probsOR.orig)
    #diffbysame <- brobgiveback(results$diffbysame)
  } else if(TYPE == "normal")
  {
    res.out <- signif( results$DiM.probsOR.orig, digits=dig)
  }
  
  cat("\n\n#########################################################################################\n")
  cat("###\n")
  cat("### ON THE DIFFERENCE IN MEANS\n")
  cat("### G.L. Bretthorst (1993)\n")
  cat("###\n### original Mathematica code by U.M. Studer (90s, Switzerland)\n")
  cat("###\n### methods = normal, log, brob\n### log and brob are for large numbers and infinity errors\n")
  cat("\nNote:\nIf any probability is printed as '1' (= one) or '0' (= zero), it means that the\nprobability is practically that value by giving respect to limited computer precision.\n")

  cat("\n...based on METHOD = ",TYPE," and CONVBACK = ",convback,"\n\n",sep="")

  # data (input) = descriptive statistics
  cat("\n-------------------------------- Data (Input) ------------------------------------------\n\n")
  if(!srib.base)
  {
    cat(paste("N_1 = ",consts.out$Ni ," :\t\t\tMean_1 ± SD_1\t\t= ", consts.out$Di," ± ",consts.out$si, "\n", sep=""))
    cat(paste("N_2 = ",consts.out$Nii," :\t\t\tMean_2 ± SD_2\t\t= ", consts.out$Dii," ± ",consts.out$sii, "\n", sep=""))
  } else
  {
    cat(paste("s_1/N_1 = ",inval$Si,"/",consts.out$Ni," :\t\t\tMean_1 ± SD_1\t\t= ", consts.out$Di," ± ",consts.out$si, "\n", sep=""))
    cat(paste("s_2/N_2 = ",inval$Sii,"/",consts.out$Nii," :\t\t\tMean_2 ± SD_2\t\t= ", consts.out$Dii," ± ",consts.out$sii, "\n", sep=""))
  }
  cat(paste("N_total = N_1 + N_2  = ", consts.out$NN ," :\t\tMean_comb ± SD_comb\t= ", consts.out$DD," ± ",consts.out$ss, "\n", sep=""))
  cat("\n")
  cat(paste("Bounds on the Mean (s_min = ",consts.out$smin,"):\t\tMean_L = ",consts.out$L, ",\tMean_H = ",consts.out$H,"\n",sep=""))
  cat(paste("Bounds on the Standard Deviation:\t  SD_L = ",consts.out$sL,",\t  SD_H = ",consts.out$sH,"\n",sep=""))
  if(consts.out$L < consts.out$DD)
  {
    cat(paste("\nMean_L - Mean_comb < 0 = ", (consts.out$L < consts.out$DD), "\t(-> '+'-sign between Gamma-fcts o.k.)","\n", sep=""))
  } else
  {
    cat(paste("\nMean_L - Mean_comb < 0 = ", (consts.out$L < consts.out$DD), "\t(-> '+'-sign between Gamma-fcts false!)","\n", sep=""))
  }    
  
  # results
  cat("\n-------------------------------- Results ------------------------------------------------\n\n")
  
  cat(paste("p(mv | D_1, D_2, I)\t\t= const. ", res.out["PMV"], "\n", sep=""))
  cat(paste("p(mbarv | D_1, D_2, I)\t\t= const. ", res.out["PMbarV"], "\n", sep=""))
  cat(paste("p(mvbar | D_1, D_2, I)\t\t= const. ", res.out["PMVbar"], "\n", sep=""))
  cat(paste("p(mbarvbar | D_1, D_2, I)\t= const. ", res.out["PMbarVbar"],"\n\n", sep=""))
  cc.x <- signif( 1/ (4 * (consts.out$H - consts.out$L) * log(consts.out$sH / consts.out$sL) * (2*pi)^(consts.out$NN/2) ), digits=dig)
  cat(paste("where\t\tconst.\t= ", cc.x," / p(D_1,D_2|I)\n\t\t\t= ",res.out["cc"], "\n",sep=""))
  
  # model probabilities
  cat(paste("\n--------------- Model --------------------------------- Probability ---------------------\n\n", sep=""))
  cat(paste("mv:\t\tSame Mean,      Same Variance:\t\t", res.out["mv"], "\n", sep=""))
  cat(paste("mbarv:\t\tDifferent Mean, Same Variance:\t\t", res.out["mbarv"], "\n", sep=""))
  cat(paste("mvbar:\t\tSame Mean,      Different Variance:\t", res.out["mvbar"], "\n", sep=""))
  cat(paste("mbarvbar:\tDifferent Mean, Different Variance:\t", res.out["mbarvbar"], "\n", sep=""))
  
  # odds ratios
  cat("\n------------------------------ Odds Ratios ----------------------------------------------\n\n")
  cat(paste("The probability the means are the same is:  ",res.out["samemeans"], "\n", sep=""))
  cat(paste("The probability the means are different is: ",res.out["diffmeans"], "\n", sep=""))
  cat(paste("The odds ratio is ",res.out[23]," to 1 in favor of ",results$diffsame.terms["diffmeans"]," means.\n", sep=""))
  cat("\n")
  cat(paste("The probability the variances are the same is:  ",res.out["samevars"],
            "\nThe probability the variances are different is: ",res.out["diffvars"], "\n", sep=""))
  cat(paste("The odds ratio is ",res.out[24]," to 1 in favor of ",results$diffsame.terms["diffvars"]," variances\n", sep=""))
  cat("\n")
  cat(paste("The probability the data sets are the same is:  ",res.out["mv"], "\n", sep=""))
  cat(paste("The probability the data sets are different is: ",res.out["diffsets"], "\n", sep=""))
  cat(paste("The odds ratio is ",res.out[25]," to 1 in favor of ",results$diffsame.terms["diffsets"]," variances.\n", sep=""))
  cat("\n-------------------------------- End ----------------------------------------------------\n\n")
  cat("\n#########################################################################################\n\n")
  
}
# call:
# UMSprint(results=DiM.results)
# UMSprint(results=DiM.results, SRIB=SRIB.results)
################################################################################

  
####################################################  
# print out DiM results - short version
.DiM_prout.short <- function(DiM.prout.short, DiM.ccs, dig=3)
{

  cat("\nShort output of 'the difference in means' based on method = ",DiM.prout.short$TYPE,":\n",sep="")
  cat("\nValues were converted back to normal = ",DiM.prout.short$convback,"\n\n",sep="")
  print(DiM.ccs$input.df, right=FALSE, row.names=FALSE, digits=dig)
  cat(paste("\n"))
  print(DiM.ccs$prior.df, right=FALSE, row.names=FALSE, digits=dig)
  cat(paste("\n"))

  print(DiM.prout.short$DiM.probs.res.df, right=FALSE, row.names=FALSE, digits=dig)
  cat(paste("\n"))
  print(DiM.prout.short$DiM.OR.res.df, right=FALSE, row.names=FALSE, digits=dig)
  cat(paste("\n"))
  print(DiM.prout.short$DiM.diffvssame.df, right=FALSE, row.names=FALSE, digits=dig)
  cat(paste("\n"))
  
  if(DiM.prout.short$TYPE != "normal" && DiM.prout.short$convback==TRUE)
  {
    cat("\nconvback == TRUE, values based on exp(x):\n")
    print( data.frame("exp(x)"=DiM.prout.short$DiM.probsOR.exp[1:6], check.names=FALSE), digits=dig)
    print( data.frame("exp(x)"=DiM.prout.short$DiM.probsOR.exp[7:16], check.names=FALSE), digits=dig)
    print( data.frame("exp(x)"=DiM.prout.short$DiM.probsOR.exp[17:22], check.names=FALSE), digits=dig)
    print( data.frame("exp(x)"=DiM.prout.short$DiM.probsOR.exp[23:25], check.names=FALSE), digits=dig)
  }
}
# call:
#.DiM_prout.short(DiM.out, DiM.ccs)
#.DiM_prout.short(DiM.out.convbackT.normal$res.short, DiM.ccs)
#.DiM_prout.short(DiM.out.convbackF, DiM.ccs, convback=FALSE)
# DiM.out <- DiM.out.convbackT
# DiM.out <- DiM.out.convbackF
################################################################################


################################################################################
# small speed test to compare methods (normal, log, brob)
.DiM_speedtest <- function(DiM.ccs, Nsteps=100)
{

  constants <- DiM.ccs[["constants"]]
  attach(constants)
    
  mat.t <- matrix(data=NA, nrow=4, ncol=3)
  rownames(mat.t) <- c("PMV","PMbarV","PMVbar","PMbarVbar")
  colnames(mat.t) <- c("normal","log","brob")
  mat.v <- mat.t

  cat("\n...processing normal\n")
  mat.t[1,"normal"] <- system.time( mat.v[1,"normal"] <- PMV.hypo(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH) )["elapsed"]
  mat.t[2,"normal"] <- system.time( mat.v[2,"normal"] <- PMbarV.hypo(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii) )["elapsed"]
  mat.t[3,"normal"] <- system.time( mat.v[3,"normal"] <- PMVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii) )["elapsed"]
  mat.t[4,"normal"] <- system.time( mat.v[4,"normal"] <- PMbarVbar.hypo(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii) )["elapsed"]

  cat("\n...processing log\n")
  mat.t[1,"log"] <- system.time( mat.v[1,"log"] <- PMV.hypo.log(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH, Nsteps=Nsteps) )["elapsed"]
  mat.t[2,"log"] <- system.time( mat.v[2,"log"] <- PMbarV.hypo.log(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps) )["elapsed"]
  mat.t[3,"log"] <- system.time( mat.v[3,"log"] <- PMVbar.hypo.log(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps) )["elapsed"]
  mat.t[4,"log"] <- system.time( mat.v[4,"log"] <- PMbarVbar.hypo.log(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps) )["elapsed"]

  cat("\n...processing brob\n")
  mat.t[1,"brob"] <- system.time( mat.v[1,"brob"] <- PMV.hypo.brob(NN=NN, dd=dd, upnum=upnum, lownum=lownum, sL=sL, sH=sH, Nsteps=Nsteps)@x )["elapsed"]
  mat.t[2,"brob"] <- system.time( mat.v[2,"brob"] <- PMbarV.hypo.brob(NN=NN, zz=zz, dd=dd, upinum=upinum, lowinum=lowinum, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)@x )["elapsed"]
  mat.t[3,"brob"] <- system.time( mat.v[3,"brob"] <- PMVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)@x )["elapsed"]
  mat.t[4,"brob"] <- system.time( mat.v[4,"brob"] <- PMbarVbar.hypo.brob(Dsi=Dsi, Dsii=Dsii, sL=sL, sH=sH, H=H, L=L, Ni=Ni, Nii=Nii, Nsteps=Nsteps)@x )["elapsed"]
  
  detach(constants)
  
  mat.comp <- rbind( mat.t[,"brob"]/mat.t[,"log"],
                     mat.t[,"brob"]/mat.t[,"normal"],
                     mat.t[,"log"]/mat.t[,"normal"]
                   )
  rownames(mat.comp) <- c("brob-vs-log","brob-vs-normal","log-vs-normal")
return( list(mat.t=mat.t,mat.v=mat.v,mat.comp=mat.comp) )
}
# call:
# 
################################################################################


###############################################################################
# plot graphical comparison
# only qualitative Bretthorst
DiM.quan.plot <- function(stype=NA,
                          v1=NA, v2=NA, # vector with raw values
                          xbar1=NA, sd1=NA, n1=NA, xbar2=NA, sd2=NA, n2=NA, # summary statistics
                          groupnames=paste("group_",c(1,2),sep=""),
                          fname="DiM_plot.quan.pdf", 
                          # limits
                          # default = +/- 3*sd                          
                          nsd=3,
                          # length of sequence / points to draw
                          l.out=100,
                          # scale before plotting?
                          center=FALSE, scale=FALSE,
                          simulate=FALSE, # makes only sense if no raw values are present
                          seed=9999,
                          fac.ext=0.1,
                          dig=2, pdfout=FALSE, legendplace="left"
                          )
{
  
  
  
  if(stype == "raw" && simulate == TRUE) stop("No simulation with raw values.")
  
  if(stype == "sumstat") # summary statistics
  {
  
    if(simulate == TRUE)
    {
      set.seed(seed)
      v1 <- rnorm(n=n1, mean=xbar1, sd=sd1)
      v2 <- rnorm(n=n2, mean=xbar2, sd=sd2)
    } else if(simulate == FALSE)
    {  
      if( any(is.na(c(xbar1,sd1,n1,xbar2,sd2,n2))) ) stop("Missing values in summary statistics present.")
      lv.stats <- data.frame(mean=c(xbar1,xbar2),
                             sd=c(sd1,sd2),
                             N=c(n1,n2)
                            )
      rownames(lv.stats) <- groupnames
      
      # limits and ranges x-axis for both groups
      lim1.r <- xbar1 + c(-1,1) * nsd * sd1
      lim2.r <- xbar2 + c(-1,1) * nsd * sd2
      xlim.r <- range(c(lim1.r,lim2.r)) * c(1-fac.ext,1+fac.ext)
      # x values
      lx1 <- seq(lim1.r[1], lim1.r[2], length.out=l.out)
      lx2 <- seq(lim2.r[1], lim2.r[2], length.out=l.out)
      # density values y axis
      ly1 <- dnorm(x=lx1, mean=lv.stats[1,"mean"], sd=lv.stats[1,"sd"])
      ly2 <- dnorm(x=lx2, mean=lv.stats[2,"mean"], sd=lv.stats[2,"sd"])
      # plot range y axis
      ylim.max <- max(ly1,ly2) * (1+fac.ext)
    }
  }
  
  if(stype == "raw" || stype == "sumstat" && simulate == TRUE) # raw values or summary stats + simulation
  {
    if( any(is.na(c(v1,v2))) ) stop("Missing data in raw values.")
    n1 <- length(v1)
    n2 <- length(v2)
    lv <- list(v1,v2)
    names(lv) <- groupnames
    
    if(center == TRUE)
    {
      lv <- lapply(lv, function(x) scale(x, center=TRUE, scale=FALSE))
    }
    if(scale == TRUE)
    {
      lv <- lapply(lv, function(x) scale(x, center=FALSE, scale=TRUE))
    } 
    
    # summary statistics 
    lv.stats <- t(sapply(lv, function(x) c(mean(x),sd(x),length(x))))
    colnames(lv.stats) <- c("mean","sd","N")
    # densities x and y axis
    lv.dens <- lapply(lv, density)
    lx1 <- lv.dens[[1]]$x
    ly1 <- lv.dens[[1]]$y
    lx2 <- lv.dens[[2]]$x
    ly2 <- lv.dens[[2]]$y
    # range x axis
    xlim.r <- range(sapply(lv.dens, function(x) x$x)) * c(1-fac.ext,1+fac.ext)
    # plot range y axis
    ylim.max <- max(sapply(lv.dens, function(x) x$y)) * (1+fac.ext)
  }
  
  cat("\nPlot quantitative values - summary statistics -\n\n")
  cat("stype\t\t= ",stype,"\n",sep="")
  cat("simulate\t= ",simulate,"\n\n",sep="")
  # print descriptive characteristics
  print(lv.stats)
  cat("\n")
  
  # Cohens d for different sds and sample sizes
  cd <- function(xbar1, sd1, n1, xbar2, sd2, n2)
  {
    spooled <- sqrt( ( (n1-1)*(sd1^2) + (n2-1)*(sd2^2) ) / (n1+n2-2) )
    return( (xbar1-xbar2) / spooled )
  }
  cd.v <- cd(lv.stats[1,"mean"],lv.stats[1,"sd"],lv.stats[1,"N"],
             lv.stats[2,"mean"],lv.stats[2,"sd"],lv.stats[2,"N"])
  
  ### plot graphical comparison p(H[delta]|S_1/N_1,I) versus p(H[delta]|S_2/N_2,I)
  if(pdfout) pdf(fname,width=9,height=6,paper="A4r")
  par(mar=c(5,6,5,5))
  plot(0,
       xlab="", ylab="probabilities",
       main="",
       type="l", lty=1, lwd=1.75, col="red", bty="l",
       xlim=xlim.r,
       ylim=c(0,ylim.max)
       )
  lines(lx1,ly1, lty=1, lwd=1.75, col="red")
  lines(lx2,ly2, lty="dashed", lwd=1.75, col="blue")
  
  mtext(expression(paste(delta," Mean Difference",sep="")), 3, line=2, cex=1.5)
  mtext(eval(substitute( expression(paste(bar(delta)["1vs2"]," = ",cd.v)),
             list(cd.v=round(cd.v,dig)))),
        3, line=0.5, cex=1)
  mtext(expression(paste(delta," = ", mu["1"]," - ",mu["2"])), 1, line=3, cex=1.5)

  if(sumstats == "v") # vertical
  {    
    mtext(eval(substitute(expression(paste(bar(x)["1"] ," = ",Di," | ",
                                           sigma["1"] ," = ",si, " | ",
                                           "n1 = ",n1, "")),
                          list(Di=round(lv.stats[1,"mean"],dig),
                               si=round(lv.stats[1,"sd"],dig),
                               n1=n1))),
                          4,line=-1, cex=0.9)
    mtext(eval(substitute(expression(paste(bar(x)["2"] ," = ",Dii," | ",
                                           sigma["2"] ," = ",sii, " | ",
                                           "n2 = ",n2, "")),
                          list(Dii=round(lv.stats[2,"mean"],dig),
                               sii=round(lv.stats[2,"sd"],dig),
                               n2=n2))),
                          4,line=0, cex=0.9)
  } else if(sumstats == "h") # horizontal
  {
    # alternative: add another nice legend with information
    legend("right", legend=c(eval(substitute(expression(paste(bar(x)["1"] ," = ",Di)),list(Di=round(lv.stats[1,"mean"],dig)))),
                             eval(substitute(expression(paste(bar(x)["2"] ," = ",Dii)),list(Dii=round(lv.stats[2,"mean"],dig)))),
                             "",
                             eval(substitute(expression(paste(sigma["1"] ," = ",si)),list(si=round(lv.stats[1,"sd"],dig)))),
                             eval(substitute(expression(paste(sigma["2"] ," = ",sii)),list(sii=round(lv.stats[2,"sd"],dig)))),
                             "",
                             eval(substitute(expression(paste("n1 = ",ni)),list(ni=n1))),
                             eval(substitute(expression(paste("n2 = ",nii)),list(nii=n2)))
                             #"",
                             #eval(substitute(expression(paste(bar(delta)["1vs2"] ," = ",cd.v)),list(cd.v=round(cd.v,dig))))
                             ),
                    text.col=c("black"),
                    bty="n", cex=1,
                    y.intersp=1.4,
                    title="statistics", title.col="darkred"
                    )
  }  
  # add a nice legend to explain curves and associated probability densities
  legend(legendplace, legend=groupnames,
         text.col=c("red","blue"), bty="n", cex=1, y.intersp=1.4,
         lty=c(1,2), col=c("red","blue"), lwd=2,
         #title="Groups", title.col="black", title.adj=1.1
         )
  if(pdfout) dev.off()

  }
# call:
# DiM.quan.plot()
#
# call with pdf output:
# DiM.quan.plot(,pdfout=TRUE,fname="DiM.quan.plot.pdf")
################################################################################


# --------------------------------------------------
# End: ON THE DIFFERENCE IN MEANS
# --------------------------------------------------
