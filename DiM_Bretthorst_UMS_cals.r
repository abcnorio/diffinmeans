### (C) 2005-2023 by LG
### R-code supplement
### to the book
###
### "Subjektive Ansichten und objektive Betrachtungen"
###
### written by LG & Huber (2023)
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
# On the difference in means
# paper: G.L. Bretthorst "On the difference of means" (1993)

#
# R code based on Mathematica code by UMStuder (90's, Zürich/ CH)
# R code by LG
# first: 12-06-05, 21-06-06, 20-04-17
# brob: 07-09-20
# use log-only: 2023-06-29
# latest: 2023-07-26

# load helpe functions and all what is necessary
source("DiM_Bretthorst_UMS.r")
source("DiM_Bretthorst_helperfuncs-general.r")


################################################################################
# example values DiM

# works
res.SIB.NRFmf <- SucRatesIntBounds(Si=20, Ni=47, Sii=13, Nii=28, smin=0, snames=c("male","female"))
inval <- res.SIB.NRFmf
inval

# example when normal version fails due to infinity errors - requires large number fix (log, Brobdingnag)
res.SIB.NRFtotal <- SucRatesIntBounds(Si=(20+13), Ni=(47+28), Sii=(338 %/% 4), Nii=338, smin=0, snames=c("male","female"))
# normal - does not work
DiM.res <- DiM(initials=res.SIB.NRFtotal, cinput="qual", cmethod="normal", prout="long", convback=FALSE, plotti=FALSE) 
# brob version
DiM.res <- DiM(initials=res.SIB.NRFtotal, cinput="qual", cmethod="brob", prout="long", convback=FALSE, plotti=FALSE) 
# log version
DiM.res <- DiM(initials=res.SIB.NRFtotal, cinput="qual", cmethod="log", prout="long", convback=FALSE, plotti=FALSE) 
# plot to see why it fails
UMSplot(inval=res.SIB.NRFtotal, pdfout=FALSE)


inval <- res.SIB.NRFtotal
inval


# UMS works
res.SIB.UMS <- SucRatesIntBounds(Si=11, Ni=15, Sii=10, Nii=16, smin=0, snames=c("male","female"))
inval <- res.SIB.UMS
inval



### step by step procedure

# constants
DiM.ccs <- .DiM_ccs(inval=inval)
DiM.ccs

# speedtest
DiM.comp <- .DiM_speedtest(DiM.ccs=DiM.ccs)
DiM.comp

# calculate integrals
DiM.integrals.normal <- .DiM_callintegrals(DiM.ccs, method="normal")
DiM.integrals.normal
#DiM.integrals <- DiM.integrals.normal

DiM.integrals.log <- .DiM_callintegrals(DiM.ccs, method="log")
exp(DiM.integrals.log)
#DiM.integrals <- DiM.integrals.log

DiM.integrals.brob <- .DiM_callintegrals(DiM.ccs, method="brob")
sapply(DiM.integrals.brob, as.numeric)
#DiM.integrals <- DiM.integrals.brob


# prepare results before printing them based on integrals
DiM.out.convbackT.normal <- .DiM_prepres(DiM.integrals.normal, method="normal", percfac=1)
DiM.out.convbackT.normal <- .DiM_prepres(DiM.integrals.normal, method="normal", convback=TRUE, percfac=1)
DiM.out.convbackT.normal
results <- DiM.out.convbackT.normal
UMSprint(results=results, ccs=DiM.ccs, inval=inval, dig=3)
 
DiM.out.convbackT.log <- .DiM_prepres(DiM.integrals.log, method="log", convback=TRUE, percfac=1)
UMSprint(results=DiM.out.convbackT.log, ccs=DiM.ccs, inval=inval, dig=3)

DiM.out.convbackF.log <- .DiM_prepres(DiM.integrals.log, method="log", convback=FALSE, percfac=1)
UMSprint(results=DiM.out.convbackF.log, ccs=DiM.ccs, inval=inval, dig=3)

DiM.out.convbackT.brob <- .DiM_prepres(DiM.integrals.brob, method="brob", convback=TRUE, percfac=1)
UMSprint(results=DiM.out.convbackT.brob, ccs=DiM.ccs, inval=inval, dig=3)

DiM.out.convbackF.brob <- .DiM_prepres(DiM.integrals.brob, method="brob", convback=FALSE, percfac=1)
UMSprint(results=DiM.out.convbackF.brob, ccs=DiM.ccs, inval=inval, dig=3)

# plot
UMSplot(inval=inval)
UMSplot(inval=inval, legendplace="topright")
################################################################################


################################################################################
# use DiM call

### EXAMPLE 1 - UMS 1998, p.47
# QUAL
initials <- list(Si=11, Ni=15, Sii=10, Nii=16,
                 snames=c("voluntary","non-voluntary"))
cinput <- "qual"
SucRatesIntBounds(Si=11, Ni=15, Sii=10, Nii=16, smin=0, snames=c("male","female"))
# QUAN
initials <- list(Ni=15, Di=0.7058824, si=0.1073966,
                 Nii=16, Dii=0.6111111, sii=0.1118397,
                 L=0.05, H=0.95, sL=0.052, sH=0.118,
                 snames=c("voluntary","non-voluntary"))
cinput <- "quan"

### EXAMPLE 2 - UMS 1998, p.48
initials <- list(Si=20, Ni=31, Sii=17, Nii=27,
                 snames=c("voluntary","non-voluntary"))
cinput <- "qual"
SucRatesIntBounds(Si=20, Ni=31, Sii=17, Nii=27, smin=0, snames=c("male","female"))
# QUAN
initials <- list(Ni=31, Di=0.6363636, si=0.08249866,
                 Nii=27, Dii=0.6206897, sii=0.08858781,
                 L=0.03, H=0.97, sL=0.029, sH=0.092,
                 snames=c("voluntary","non-voluntary"))
cinput <- "quan"



# choose method
cmethod <- "normal"
cmethod <- "log"
cmethod <- "brob"

# calls
DiM.res <- DiM(initials=initials, cinput=cinput, cmethod="normal", prout="short", convback=FALSE, plotti=TRUE) 
DiM.res <- DiM(initials=initials, cinput=cinput, cmethod="normal", prout="short", convback=FALSE, plotti=FALSE)
DiM.res <- DiM(initials=initials, cinput=cinput, cmethod="normal", prout="long", convback=FALSE, plotti=FALSE)

DiM.res <- DiM(initials=initials, cinput=cinput, cmethod="log", prout="long", convback=TRUE, plotti=FALSE)
DiM.res <- DiM(initials=initials, cinput=cinput, cmethod="log", prout="long", convback=FALSE, plotti=FALSE)

DiM.res <- DiM(initials=initials, cinput=cinput, cmethod="brob", prout="long", convback=TRUE, plotti=FALSE)
DiM.res <- DiM(initials=initials, cinput=cinput, cmethod="brob", prout="long", convback=FALSE, plotti=FALSE)


################################################################################
# quantitative example + plot

# example 1
# UMS, s.a.
initials <- list(Ni=15, Di=0.7058824, si=0.1073966,
                 Nii=16, Dii=0.6111111, sii=0.1118397,
                 L=0.05, H=0.95, sL=0.052, sH=0.118,
                 snames=c("voluntary","non-voluntary"))

attach(initials)
xbar1 <- Di
sd1 <- si
n1 <- Ni
xbar2 <- Dii
sd2 <- sii
n2 <- Nii
groupnames <- snames
detach(initials)

xbar1
sd1
n1
xbar2
sd2
n2
groupnames

DiM.quan.plot(stype="sumstat", xbar1=xbar1, sd1=sd1, n1=n1, xbar2=xbar2, sd2=sd2, n2=n2, simulate=FALSE, groupnames=groupnames)
DiM.quan.plot(stype="sumstat", xbar1=xbar1, sd1=sd1, n1=n1, xbar2=xbar2, sd2=sd2, n2=n2, simulate=TRUE, groupnames=groupnames)
nfac <- 100
DiM.quan.plot(stype="sumstat", xbar1=xbar1, sd1=sd1, n1=n1*nfac, xbar2=xbar2, sd2=sd2, n2=n2*nfac, simulate=TRUE, groupnames=groupnames)


# example 2
# group1 
xbar1 <-4
sd1 <- 2
n1 <- 24

# group 2
xbar2 <- 3.8 
sd2 <- 2.2
n2 <- 28

set.seed(1423)
v1 <- rnorm(n=n1, mean=xbar1, sd=sd1)
v2 <- rnorm(n=n2, mean=xbar2, sd=sd2)
groupnames <- c("test","control")

DiM.quan.plot(stype="sumstat", xbar1=xbar1, sd1=sd1, n1=n1, xbar2=xbar2, sd2=sd2, n2=n2, simulate=FALSE, groupnames=groupnames)
DiM.quan.plot(stype="sumstat", xbar1=xbar1, sd1=sd1, n1=n1, xbar2=xbar2, sd2=sd2, n2=n2, simulate=TRUE, groupnames=groupnames)
nfac <- 100
DiM.quan.plot(stype="sumstat", xbar1=xbar1, sd1=sd1, n1=n1*nfac, xbar2=xbar2, sd2=sd2, n2=n2*nfac, simulate=TRUE, groupnames=groupnames)
DiM.quan.plot(stype="raw", v1=v1, v2=v2, groupnames=groupnames)


################################################################################


# NOT RUN BELOW THIS POINT
################################################################################
# perform several speedtests
SAVE <- FALSE
reps <- 20
mat <- matrix(data=NA, ncol=4, nrow=reps)
colnames(mat) <- c("PMV","PMbarV","PMVbar","PMbarVbar")
mat[1,] <- .DiM_speedtest(DiM.ccs=DiM.ccs)$mat.comp[1,]
for(i in 2:reps)
{
  print(i)
  mat[i,] <- .DiM_speedtest(DiM.ccs=DiM.ccs)$mat.comp[1,]
}
mat.fn <- rbind( apply(mat,2,fivenum), apply(mat,2,mean), apply(mat,2,sd) )
rownames(mat.fn) <- c("min", "lH","med","uH","max","mean","sd")
t(mat.fn)
if(SAVE) save(mat.fn, file=paste("DiM_reps-",reps,sep=""))
################################################################################


################################################################################
# TODO run simulation grid approach

# check success_1 of total_1 vs success_2 of total_2
size <- 10
tab <- expand.grid(si=1:size, Ni=1:size, sii=1:size, Nii=1:size)
dim.tab <- dim(tab)
dim.tab[1]/3600
head(tab)

################################################################################
