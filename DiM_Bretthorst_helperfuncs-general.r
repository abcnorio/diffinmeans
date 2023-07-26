################################################################################
# --------------------------------------------------
# ON THE DIFFERENCE IN MEANS
# --------------------------------------------------
# based on G.L. Bretthorst (1993) "On the difference of means"
#
# helper functions


################################################################################ 
# add two values/ vector on the log scale
.llog.2add.short <- function(v, log1p.crit=1e-09) #1e-15
{
  max.v <- max(v)
  sum.v <- sum(exp(v - max.v))
  ifelse(sum.v < log1p.crit,
         return( max.v + log1p(sum.v) ),
         return( max.v + log(sum.v) )
        )
}
# call:
# 
################################################################################ 


################################################################################
# subtract log values  
.llog.2sub.short <- function(v)
{
  return( max(v) + log1p(-exp(-abs(v[1]-v[2]))) )
}
# call:
# 
################################################################################


################################################################################
# Gamma funcs
# (generalized) incomplete gamma function
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
    require(Brobdingnag)
    if(is.na(z1)) return( pgamma(z0,a, lower=FALSE, log=TRUE))
    else return( brob(pgamma(z0,a, lower=FALSE,log=TRUE)) - brob(pgamma(z1,a, lower=FALSE,log=TRUE)) )
  }  
}
# call:
# 
################################################################################ 


################################################################################ 
gamma_inc_gen.alt <- function(a,z0,z1=NA, log=FALSE, lower=FALSE)
{
  
  if(is.na(z1))
  {
    res <- pgamma(z0,a,lower=lower,log=log)
    ifelse(log, return(brob), return(res) )
  } else
  {
    p1 <- pgamma(z0,a,lower=lower,log=log)
    p2 <- pgamma(z1,a,lower=lower,log=log)
  }
  ifelse(log, return(brob(p1)-brob(p2)), return(p1-p2))
}
# call:
# 
################################################################################ 


################################################################################ 
# simple Simpson rule for
# normal, log, brob functions ie. output
#
#simpsonrule.nlb
#
# call:
#
#f <- function(x) sin(x)
#f.log <- function(x) log(sin(x))
#f.brob <- function(x) brob(log(sin(x)))
##
#lower <- 0
#upper <- pi
#Nsteps <- 100
##
#sprintf("%.16e", simpsonrule.nlb(fx=f, lower=lower, upper=upper, method="normal",Nsteps=Nsteps) )
#sprintf("%.16e", exp(simpsonrule.nlb(fx=f.log, lower=lower, upper=upper, method="log",Nsteps=Nsteps)) )
#sprintf("%.16e", as.numeric(simpsonrule.nlb(fx=f.brob, lower=lower, upper=upper, method="brob",Nsteps=Nsteps)) )
#sprintf("%.16e", integrate(f,lower,upper)$v )
################################################################################ 


#########################################################################
# helper function to replace '%*%' scalar product that does not work for brob objects
#
scalarprod.brob <- function(c1,c2)
{
  return( sum(c1*c2) )
}
# call:
# 
#########################################################################


#########################################################################
# helper function to convert a brob list to a brob vector
list2vec.brob <- function(list.brob, check.list=FALSE)
{
  if(!is.list(list.brob)) stop("Input is not a list.")
  if(check.list)
  {
    if(any(unlist(lapply(list.brob, function(x) attr(x, "class"))) != "brob"))
    {
      stop("There are non-brob elements in the list. Stopping.")
    }  
  }
return( brob( sapply(list.brob,getX), sapply(list.brob,getP) ) )
}
# call:
# 
#########################################################################


################################################################################
# extract and format brob numbers from brob list as character
brobgiveback <- function(br,names=NA,dig=3)
{
  res <- paste("[",ifelse(sapply(br, getP),"+","-"),"] exp(", signif(sapply(br, getX),digits=dig),")",sep="")
  if(sum(is.na(names))==0 && length(res)==length(names))
  {
    names(res) <- names
  } else names(res) <- names(br)
return(res)
}
# call:
# brobgiveback(DiM.integrals)
################################################################################


################################################################################
# Simpson rule for normal, log, and brob
simpsonrule.nlb <-  function(fx,
                             lower, upper,
                             method="normal", # log, brob
                             eps=.Machine$double.xmin,
                             #eps=.Machine$double.eps,
                             log1p.crit=1e-09, # 1e-15
                             Nsteps=100,
                             parallel=FALSE,
                             fixINF=TRUE, ...)
{
  
  if(lower == upper)
  {
    cat("\nLower and upper limits of the integral are identical.\nNo calculations done.\n")
    return(0)
  } else if(lower > upper)
  {
    cat("\nBe aware of the limits chosen:\tlower > upper\n.Negative value of the integral will result.\n")
    return( -1 * simpsonrule.nlb(fx,
                                 upper, lower,
                                 method=method,
                                 eps=eps,
                                 log1p.crit=log1p.crit,
                                 Nsteps=Nsteps,
                                 parallel=parallel,
                                 fixINF=fixINF) )
  }    
  
  sek.l <- 2*Nsteps+1
  hml.s <- sign(upper-lower)
  
  # avoid infinite values (zero with log(0), etc.)
  if(fixINF == TRUE)
  {
    if(method %in% c("normal","log"))
    {
      if(is.infinite(fx(lower))) lower <- lower + eps * hml.s
      if(is.infinite(fx(upper))) upper <- upper - eps * hml.s      
    } else if(method == "brob")
    {
      if(is.infinite(fx(lower)@x)) lower <- lower + eps * hml.s
      if(is.infinite(fx(upper)@x)) upper <- upper - eps * hml.s
    }
  }
  sek <- seq(lower,upper,length=sek.l)
  
  if(method == "normal")
  {
    xfxap <- fx(sek)
    h <- xfxap[2] - xfxap[1]
    res <- sum(h*( xfxap[2 * (1:Nsteps) - 1] +
                     4*xfxap[2 * (1:Nsteps)] +
                     xfxap[2 * (1:Nsteps) + 1] )
               /3)
  } else if(method == "log")
  {
    if(parallel)
    {
      require(parallel)
      mccores <- detectCores(all.tests=FALSE, logical=TRUE)
      fx.log <- simplify2array(mclapply(seq_along(1:sek.l),
                                        function(i) fx(sek[i]),
                                        mc.cores=mccores
      ) )
    } else
    {
      fx.log <- Vectorize(fx)(sek)
    }
    h.log <- log(sek[2] - sek[1])
    s.log <- c( fx.log[2 * (1:Nsteps) - 1],
                log(4) + fx.log[2 * (1:Nsteps)],
                fx.log[2 * (1:Nsteps) + 1]
    )
    max.s.log <- max(s.log)
    sum.log <- sum(exp(s.log - max.s.log))
    #res <- h.log - log(3) + max.s.log + log1p(sum(exp(s.log[-which.max(s.log)] - max.s.log)))
    #
    # when do log1p(x) vs. log(x)
    # https://scicomp.stackexchange.com/questions/20629/when-should-log1p-and-expm1-be-used
    # IF 1+x=1 in floating point accuracy
    # in case of x = 1e-15 which is -34.53878 on log scale
    # ie. usual double precision arithmetic
    # if is already useful for
    # x = 1e−09 which is -20.72327 on log scale
    ifelse(sum.log < log1p.crit,
           res <- h.log - log(3) +  max.s.log + log1p(sum.log),
           res <- h.log - log(3) +  max.s.log + log(sum.log)
    )
    
  } else if(method=="brob")
  {
    fx.brob <- list2vec.brob(lapply(seq_along(1:sek.l), function(i) fx(sek[i])))
    h <- as.brob(sek[2] - sek[1]) # diff
    res <- sum(h/as.brob(3) * ( fx.brob[2 * (1:Nsteps) - 1] +
                                  4 * fx.brob[2 * (1:Nsteps)] +
                                  fx.brob[2 * (1:Nsteps) + 1]
    ) )
  }
  return(res)
}
################################################################################
