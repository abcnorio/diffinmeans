################################################################################
# ON THE DIFFERENCE OF MEANS
################################################################################

source("DiM_Bretthorst_PG.r")
source("DiM_Bretthorst_UMS.r")
source("DiM_Bretthorst_helperfuncs-general.r")

# call PG scheme
inputvalues <- list(snames = c("riverB.1","riverB.2"),
                    # sample 1
                    d1 = c(13.2,13.8,8.7,9,8.6,9.9,14.2,9.7,10.7,8.3,8.5,9.2),
                    # sample 2
                    d2 = c(8.9,9.1,8.3,6,7.7,9.9,9.9,8.9),
                    
                    # Input priors and no. of steps in evaluation of p(r|D_1,D_2,I) & p(delta|D_1,D_2,I)
                    # ndelta = number of steps in delta parameter (mean difference)
                    ndelta = 1000, #100
                    # nr = number of steps in r parameter (ratio of the standard deviations)
                    nr = 1000, # 100
                    
                    # Set prior limits (assumed the same for each data set) on mean (low,high),
                    # and prior limits (assumed the same for each data set) on the
                    # standard deviation (sigmalow, sigmahigh).
                    # upper mean
                    high = 12,
                    # lower mean
                    low = 7,
                    # upper sd
                    sigma.high = 4,
                    # lower sd
                    sigma.low = 1)
inputvalues
# normal
dim.res <- DiM.pg(invtyp="pg", inputvalues, print.res=TRUE)
DiM.print.pg(dim.res)
dim.res.newlimits <- DiM.extract.limits(dim.res, scaleL=30, scaleH=4, change=TRUE)
DiM.extract.limits(dim.res.newlimits, change=FALSE)
dim.res.calc <- DiM.plot.calc.pg(dim.res.newlimits, type="normal")
DiM.plot.pg(dim.res.calc, filling=FALSE, type="normal")

# brob
dim.res <- DiM.pg(invtyp="pg", inputvalues, print.res=TRUE, type="brob")
DiM.print.pg(dim.res)

# call according to UMS scheme
inputvalues.ums <- list(snames=c("Jaynes.1","Jaynes.2"), si=6.48, Ni=4, sii=7.48, Nii=9, Di=50, Dii=42, L=34, H=58, sL=3, sH=10, ndelta=1000, nr=1000)
dim.res <- DiM.pg(invtyp="ums", inputvalues=inputvalues.ums, print.res=TRUE)

# works
res.SIB <- SucRatesIntBounds(Si=11, Ni=15, Sii=10, Nii=16, smin=0, snames=c("voluntary","non-voluntary"))
res.SIB
#
DIM.pg.res <- DiM.pg(invtyp="ums", inputvalues=res.SIB, print.res=TRUE, type="normal")
DIM.pg.res
DiM.res <- DIM.pg.res
DiM.print.pg(DIM.pg.res)
# ratio of SD requires input values scaleL (low) and scaleH (high), otherwise the script breaks
DiM.plotvalues.res.nonbrob <- DiM.plot.calc.pg(DIM.pg.res, scaleL=2, scaleH=8, type="normal")
DiM.plot.pg(DiM.plotvalues.res.nonbrob, filling=TRUE, type="normal")


# as brob
DIM.pg.res <- DiM.pg(invtyp="ums", inputvalues=res.SIB, print.res=TRUE, type="brob")
DIM.pg.res
DiM.print.pg(DIM.pg.res)
DiM.plotvalues.res.nonbrob.brob <- DiM.plot.calc.pg(DIM.pg.res, scaleL=2, scaleH=8, type="normal")
DiM.plot.pg(DiM.plotvalues.res.nonbrob.brob, filling=TRUE, type="normal")


# normally does not work
res.SIB.NRFtotal <- SucRatesIntBounds(Si=(20+13), Ni=(47+28), Sii=(338 %/% 4), Nii=338, smin=0, snames=c("male","female"))
res.SIB.NRFtotal
# results in INF values in integrals - script breaks without BROB
DIM.pg.res.brob <- DiM.pg(invtyp="ums", inputvalues=res.SIB.NRFtotal, print.res=TRUE, type="normal")
# works with BROB
DIM.pg.res.brob <- DiM.pg(invtyp="ums", inputvalues=res.SIB.NRFtotal, print.res=TRUE, type="brob")
# DIM.pg.res.brob <- DiM.pg(invtyp="ums", inputvalues=res.SIB.NRFtotal, print.res=TRUE, BROB=TRUE)
DiM.print.pg(DIM.pg.res.brob)
# try several limits...
DiM.extract.limits(DIM.pg.res.brob, scaleL=10, scaleH=1, change=FALSE)
DiM.extract.limits(DIM.pg.res.brob, scaleL=10, scaleH=2, change=FALSE)
DiM.extract.limits(DIM.pg.res.brob, scaleL=5, scaleH=2, change=FALSE)
DiM.extract.limits(DIM.pg.res.brob, scaleL=100, scaleH=4, change=FALSE)
# apply new limits
DiM.newlimits <- DiM.extract.limits(DIM.pg.res.brob, scaleL=10, scaleH=1, change=TRUE)
DiM.newlimits <- DiM.extract.limits(DIM.pg.res.brob, scaleL=100, scaleH=4, change=TRUE)
# check
DiM.extract.limits(DiM.newlimits, change=FALSE)
# calc plot values
DiM.newlimits.calc.plot <- DiM.plot.calc.pg(DiM.newlimits, type="brob")
# plot
DiM.plot.pg(DiM.newlimits.calc.plot, filling=FALSE, by1=TRUE, type="brob")

