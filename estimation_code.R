
# Estmation code

system("R CMD SHLIB cgj_new.c biv-nt.c")
dyn.load(paste("cgj_new",.Platform$dynlib.ext,sep=""))

######################################
# Montserrat Fuentes' program  (http://www4.stat.ncsu.edu/~fuentes/)
# to compute geodetic distance
rdistearth<- function(loc1, loc2, miles = FALSE ){
  if (miles) 
    R <- 3963.34
  else R <- 6378.388	
  if(missing(loc2))
    loc2 <- loc1
  R <- 6371
  lat <- loc1[, 2]
  lon <- loc1[, 1]
  coslat1 <- cos((lat * pi)/180)
  sinlat1 <- sin((lat * pi)/180)
  coslon1 <- cos((lon * pi)/180)
  sinlon1 <- sin((lon * pi)/180)
  lat <- loc2[, 2]
  lon <- loc2[, 1]
  coslat2 <- cos((lat * pi)/180)
  sinlat2 <- sin((lat * pi)/180)
  coslon2 <- cos((lon * pi)/180)
  sinlon2 <- sin((lon * pi)/180)
  PP1 <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1)
  PP2 <- cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2)
  pp <- (PP1 %*% t(PP2))
  R * acos(ifelse(pp > 1, 1, pp))
}
# transform the coordinates
# Montserrat Fuentes' projection on the plane Centered around the center of gravity
### lon.lat represent the matrix of the logitude and latitude of the data
lonlat.to.planar<-function(lon.lat, miles =FALSE){
  x <- lon.lat[, 1]
  y <- lon.lat[, 2]
  mx <- mean(x)
  my <- mean(y)
  temp <- cbind(rep(mx, 2), range(y))
  sy <- rdistearth(temp)[2, 1]
  temp <- cbind(range(x), rep(my, 2))
  sx <- rdistearth(temp ,miles = miles)[2, 1]
  temp <- list(x = sx/(max(x) - min(x)), y = sy/(max(y) - min(y)))
  COORD=cbind((x - mx) * temp$x, (y - my) * temp$y)
  return(COORD)
}
############ empirical transformation
data2frechet<-function(data)
{
 n.site<-dim(data)[2]
 n.observation<-dim(data)[1]
 newdata<-matrix(NA,n.observation,n.site)

 for (i in 1: n.site) {
   newdata[,i]<-gev2frech(data[,i],emp=TRUE)	
 }
 return(newdata)
}	

### changed to inverted data
tx<- function(x){
  return(-1/log(1-exp(-1/x))) 
}

overlarea<-function(d,radius)
  #Area of Overlapping Circles of the same radius 
{
  n<-length(d)
  val<-numeric(n)
  tmp<-.C("overlareawrap", d = as.double(d),  radius=as.double(radius),val=as.double(val), n =as.integer(n),NAOK = TRUE, DUP=TRUE)
  return(tmp$val)
  
}

distance.point.to.set<-function(x,y){
  #  Calculating all distances between one point  and a group of points
  # x point (vector)
  # y group of points (matrix each row represent a point)
  
  n<-dim(y)[1]
  x1<-x[1]
  x2<-x[2]
  y1<-y[,1]
  y2<-y[,2]
  d<-numeric(n)
  tmp<-.C("dispointset", x1 = as.double(x1), x2 = as.double(x2),y1 = as.double(y1),y2 = as.double(y2),  
          n = as.integer(n), d=as.double(d), NAOK = TRUE, DUP=TRUE)
  return(tmp$d)
  
}



check.corrmodel<-function(corrmodel)
{
  if (corrmodel == 'exponential') cr <-1
  if (corrmodel == 'wave') cr <-2
  if (corrmodel == 'spherical') cr <-3
  if (corrmodel == 'cauchy') cr <-4
  return(cr)
}

###############################################################################################
####################### For smith model
# pairwise censored likelihood for smith
smith.pwl.censored<-function(param, y, x1, x2,  nsites, ntimes, delta=Inf,threshold) 
{
  val <- 0
  sigma<-param[1] 
  
  tmp<-.C("smithpwl_censored",y = as.double(y), x1 = as.double(x1), x2 = as.double(x2),  
          n = as.integer(nsites), ntimes = as.integer(ntimes),sigma = as.double(sigma),
          val = as.double(val), delta = as.double(delta), threshold = as.double(threshold), NAOK = TRUE, DUP=TRUE)
  
  return(-tmp$val)	
}


# Fit smith using pairwise censored likelihood 
smith.fit.censored<-function(data, coords, param = NULL, delta=Inf,
                             prob =NULL,method ="nlminb",lower=NULL,upper=NULL,control = list(kkt = TRUE),std=TRUE, m=1) {
  
  # data matrix of data each row contains
  # coords matrix of coordinates
  # m number of observations per block
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  nsites <- dim(data)[2]
  ntimes <- dim(data)[1]
  
  y <- as.numeric(t(data))
  
  threshold<--1/log(prob)	
  res<-list()
  res$nsites<-nsites
  res$ntimes<-ntimes		
  res$start<-param
  res$prob<-prob
  res$threshold<-threshold	
  parnames<-c("sigma")	
  if (std){ 
    print("maximum pairwise likelihood with std: first step")
  }
  else {
    print("maximum pairwise likelihood without std")	
  }
  
  a <-optimx(par=param, smith.pwl.censored, y = y, x1 = x1, x2 = x2, nsites = nsites,
             ntimes =ntimes,delta = delta,threshold=threshold,method =method,lower=lower,upper=upper,control =control)
  res$negplik <- a$value	
  res$hessian<-attr(a,"details")[1,]$nhatend
  res$grad<-attr(a,"details")[1,]$ngatend
  res$param<-as.numeric(coef(a))
  h.hat<-res$hessian
  res$delta<-delta
  
  if (std) 
  {
    
    print("maximum pairwise likelihood with std:  standard errors")	
    
    inv.h.hat<-ginv(h.hat)
    res$h.hat<-h.hat
    j.hat<-j.hat.smith.censored.func(res$param,data=data,coords=coords, delta=delta, m=m, threshold=threshold)
    res$j.hat<-j.hat
    pen<-sum(diag(j.hat%*%inv.h.hat))
    res$pen<-pen
    var.hat<-inv.h.hat%*%j.hat%*%inv.h.hat
    res$var.hat<-var.hat	
    std.err<-sqrt(diag(var.hat))
    res$std.err<-std.err
    res$clic<-2*(res$negplik+pen)		
  }
  
  return(res)	
}

j.hat.smith.censored.func<-function(param,data,coords, delta, m=1, threshold) {
  
  # m number of observations per block
  #
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  ntimes <- dim(data)[1]
  nsites <- dim(data)[2]
  j.hat<-matrix(0,length(param),length(param))
  nblocks<-seq(1,ntimes-m+1)
  mean.dlpwl<-rep(0,length(param))
  for (i in nblocks)
  {
    y <- as.numeric(t(data[i:(i+m-1),]))
    
    
    dlpwl<-grad(smith.pwl.censored,param,y = y, x1 = x1, x2 = x2, nsites = nsites, ntimes = m,
                delta = delta,threshold=threshold)
    
    mean.dlpwl<-mean.dlpwl+dlpwl/length(nblocks)
    j.hat<-j.hat+dlpwl%*%t(dlpwl)/length(nblocks)
  }
  
  j.hat<-j.hat-mean.dlpwl%*%t(mean.dlpwl)
  j.hat<-j.hat/m
  
  
  return(j.hat)
}



####################### For inverted smith model

# pairwise censored likelihood for smith model
invsmith.pwl.censored<-function(param, y, x1, x2,  nsites, ntimes, delta=Inf,threshold) 
{
  val <- 0
  sigma<-param[1] 
  
  tmp<-.C("invsmithpwl_censored",y = as.double(y), x1 = as.double(x1), x2 = as.double(x2),  
          n = as.integer(nsites), ntimes = as.integer(ntimes),sigma = as.double(sigma),
          val = as.double(val), delta = as.double(delta), threshold = as.double(threshold), NAOK = TRUE, DUP=TRUE)
  
  return(-tmp$val)	
}


# Fit smith model using pairwise censored likelihood 
invsmith.fit.censored<-function(data, coords, param = NULL, delta=Inf,
                                prob =NULL,method ="nlminb",lower=NULL,upper=NULL,control=list(kkt=TRUE),std=TRUE, m=1) {
  
  # data matrix of data each row contains
  # coords matrix of coordinates
  # m number of observations per block
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  nsites <- dim(data)[2]
  ntimes <- dim(data)[1]
  
  y <- as.numeric(t(data))
  
  threshold<--1/log(prob)	
  res<-list()
  res$nsites<-nsites
  res$ntimes<-ntimes		
  res$start<-param
  res$prob<-prob
  res$threshold<-threshold	
  parnames<-c("sigma")		
  if (std){ 
    print("maximum pairwise likelihood with std: first step")
  }
  else {
    print("maximum pairwise likelihood without std")	
  }
  
  a <-optimx(par=param, invsmith.pwl.censored, y = y, x1 = x1, x2 = x2, nsites = nsites,
             ntimes =ntimes,delta = delta,threshold=threshold,method =method,lower=lower,upper=upper,control=control)
  res$negplik <- a$value	
  res$hessian<-attr(a,"details")[1,]$nhatend
  res$grad<-attr(a,"details")[1,]$ngatend
  res$param<-as.numeric(coef(a))
  h.hat<-res$hessian
  res$delta<-delta
  
  if (std) 
  {
    
    print("maximum pairwise likelihood with std:  standard errors")	
    
    inv.h.hat<-ginv(h.hat)
    res$h.hat<-h.hat
    j.hat<-j.hat.invsmith.censored.func(res$param,data=data,coords=coords, delta=delta, m=m, threshold=threshold)
    res$j.hat<-j.hat
    pen<-sum(diag(j.hat%*%inv.h.hat))
    res$pen<-pen
    var.hat<-inv.h.hat%*%j.hat%*%inv.h.hat
    res$var.hat<-var.hat	
    std.err<-sqrt(diag(var.hat))
    res$std.err<-std.err
    res$clic<-2*(res$negplik+pen)		
  }
  
  return(res)	
}



j.hat.invsmith.censored.func<-function(param,data,coords, delta, m=1, threshold) {
  # m number of observations per block
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  ntimes <- dim(data)[1]
  nsites <- dim(data)[2]
  j.hat<-matrix(0,length(param),length(param))
  nblocks<-seq(1,ntimes-m+1)
  mean.dlpwl<-rep(0,length(param))
  for (i in nblocks)
  {
    y <- as.numeric(t(data[i:(i+m-1),]))
    
    
    dlpwl<-grad(invsmith.pwl.censored,param,y = y, x1 = x1, x2 = x2, nsites = nsites, ntimes = m,
                delta = delta,threshold=threshold)
    
    mean.dlpwl<-mean.dlpwl+dlpwl/length(nblocks)
    j.hat<-j.hat+dlpwl%*%t(dlpwl)/length(nblocks)
  }
  
  j.hat<-j.hat-mean.dlpwl%*%t(mean.dlpwl)
  j.hat<-j.hat/m
  
  
  return(j.hat)
}



####################### For TEG model
davghol.pwl.censored<-function(param, y, x1, x2,  nsites, ntimes, delta=Inf, corrmodel = "exponential",threshold) 
{
  # pairwise censored likelihood for  Truncated Extremal Gaussian process (TEG)
  
  
  val <- 0
  theta<-param[1]
  radius<-param[2]	
  cr<-check.corrmodel(corrmodel)
  
  
  tmp<-.C("davgholpwl_censored",y = as.double(y), x1 = as.double(x1), x2 = as.double(x2),  
          n = as.integer(nsites), ntimes = as.integer(ntimes),theta = as.double(theta), radius = as.double(radius),
          val = as.double(val), delta = as.double(delta), corrmodel = as.integer(cr), threshold = as.double(threshold), NAOK = TRUE, DUP=TRUE)
  
  return(-tmp$val)	
}


# Fit Truncated Extremal Gaussian process (TEG) using pairwise censored likelihood 
davghol.fit.censored<-function(data, coords, param = NULL, delta=Inf, corrmodel='exponential' ,
                               prob =NULL,method ="nlminb",lower=NULL,upper=NULL,control=list(kkt=TRUE),std=TRUE, m=1) {
  
  # data matrix of data each row contains
  # coords matrix of coordinates
  # m number of observations per block
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  nsites <- dim(data)[2]
  ntimes <- dim(data)[1]
  
  y <- as.numeric(t(data))
  
  threshold<--1/log(prob)	
  res<-list()
  res$nsites<-nsites
  res$ntimes<-ntimes		
  res$start<-param
  res$prob<-prob
  res$threshold<-threshold	
  parnames<-c("theta","radius")	
  if (std){ 
    print("maximum pairwise likelihood with std: first step")
  }
  else {
    print("maximum pairwise likelihood without std")	
  }
  
  a <-optimx(par=param, davghol.pwl.censored, y = y, x1 = x1, x2 = x2, nsites = nsites,
             ntimes =ntimes,delta = delta, corrmodel = corrmodel,threshold=threshold,method =method,lower=lower,upper=upper,control=control)
  res$negplik <- a$value	
  res$hessian<-attr(a,"details")[1,]$nhatend
  res$grad<-attr(a,"details")[1,]$ngatend
  res$param<-as.numeric(coef(a))
  h.hat<-res$hessian
  res$corrmodel<-corrmodel
  res$delta<-delta
  
  if (std) 
  {
    
    print("maximum pairwise likelihood with std:  standard errors")	
    
    inv.h.hat<-ginv(h.hat)
    res$h.hat<-h.hat
    j.hat<-j.hat.davghol.censored.func(res$param,data=data,coords=coords, delta=delta, m=m,corrmodel=corrmodel, threshold=threshold)
    res$j.hat<-j.hat
    pen<-sum(diag(j.hat%*%inv.h.hat))
    res$pen<-pen
    var.hat<-inv.h.hat%*%j.hat%*%inv.h.hat
    res$var.hat<-var.hat	
    std.err<-sqrt(diag(var.hat))
    res$std.err<-std.err
    res$clic<-2*(res$negplik+pen)		
  }
  
  return(res)	
}


j.hat.davghol.censored.func<-function(param,data,coords, delta, m=1, corrmodel, threshold) {
  
  # m number of observations per block
  #
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  ntimes <- dim(data)[1]
  nsites <- dim(data)[2]
  j.hat<-matrix(0,length(param),length(param))
  nblocks<-seq(1,ntimes-m+1)
  mean.dlpwl<-rep(0,length(param))
  for (i in nblocks)
  {
    y <- as.numeric(t(data[i:(i+m-1),]))
    
    
    dlpwl<-grad(davghol.pwl.censored,param,y = y, x1 = x1, x2 = x2, nsites = nsites, ntimes = m,
                delta = delta, corrmodel = corrmodel,threshold=threshold)
    
    mean.dlpwl<-mean.dlpwl+dlpwl/length(nblocks)
    j.hat<-j.hat+dlpwl%*%t(dlpwl)/length(nblocks)
  }
  
  j.hat<-j.hat-mean.dlpwl%*%t(mean.dlpwl)
  j.hat<-j.hat/m
  
  return(j.hat)
}






####################### For brown-resnick model
# pairwise censored likelihood for brown-resniek
brown.pwl.censored<-function(param, y, x1, x2,  nsites, ntimes, delta=Inf,threshold) 
{
  val <- 0
  sigma<-param[1] 
  theta<-param[2]	 
  
  tmp<-.C("brownpwl_censored",y = as.double(y), x1 = as.double(x1), x2 = as.double(x2),  
          n = as.integer(nsites), ntimes = as.integer(ntimes),sigma = as.double(sigma), theta = as.double(theta),
          val = as.double(val), delta = as.double(delta), threshold = as.double(threshold), NAOK = TRUE, DUP=TRUE)
  
  return(-tmp$val)	
}


# Fit brown-resniek using pairwise censored likelihood 
brown.fit.censored<-function(data, coords, param = NULL, delta=Inf,
                             prob =NULL,method ="nlminb",lower=NULL,upper=NULL,control = list(kkt = TRUE),std=TRUE, m=1) {
  
  # data matrix of data each row contains
  # coords matrix of coordinates
  # m number of observations per block
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  nsites <- dim(data)[2]
  ntimes <- dim(data)[1]
  
  y <- as.numeric(t(data))
  
  threshold<--1/log(prob)	
  res<-list()
  res$nsites<-nsites
  res$ntimes<-ntimes		
  res$start<-param
  res$prob<-prob
  res$threshold<-threshold	
  parnames<-c("sigma","theta")	
  if (std){ 
    print("maximum pairwise likelihood with std: first step")
  }
  else {
    print("maximum pairwise likelihood without std")	
  }
  
  a <-optimx(par=param, brown.pwl.censored, y = y, x1 = x1, x2 = x2, nsites = nsites,
             ntimes =ntimes,delta = delta,threshold=threshold,method =method,lower=lower,upper=upper,control =control)
  res$negplik <- a$value	
  res$hessian<-attr(a,"details")[1,]$nhatend
  res$grad<-attr(a,"details")[1,]$ngatend
  res$param<-as.numeric(coef(a))
  h.hat<-res$hessian
  res$delta<-delta
  
  if (std) 
  {
    
    print("maximum pairwise likelihood with std:  standard errors")	
    
    inv.h.hat<-ginv(h.hat)
    res$h.hat<-h.hat
    j.hat<-j.hat.brown.censored.func(res$param,data=data,coords=coords, delta=delta, m=m, threshold=threshold)
    res$j.hat<-j.hat
    pen<-sum(diag(j.hat%*%inv.h.hat))
    res$pen<-pen
    var.hat<-inv.h.hat%*%j.hat%*%inv.h.hat
    res$var.hat<-var.hat	
    std.err<-sqrt(diag(var.hat))
    res$std.err<-std.err
    res$clic<-2*(res$negplik+pen)		
  }
  
  return(res)	
}

j.hat.brown.censored.func<-function(param,data,coords, delta, m=1, threshold) {
  
  # m number of observations per block
  #
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  ntimes <- dim(data)[1]
  nsites <- dim(data)[2]
  j.hat<-matrix(0,length(param),length(param))
  nblocks<-seq(1,ntimes-m+1)
  mean.dlpwl<-rep(0,length(param))
  for (i in nblocks)
  {
    y <- as.numeric(t(data[i:(i+m-1),]))
    
    
    dlpwl<-grad(brown.pwl.censored,param,y = y, x1 = x1, x2 = x2, nsites = nsites, ntimes = m,
                delta = delta,threshold=threshold)
    
    mean.dlpwl<-mean.dlpwl+dlpwl/length(nblocks)
    j.hat<-j.hat+dlpwl%*%t(dlpwl)/length(nblocks)
  }
  
  j.hat<-j.hat-mean.dlpwl%*%t(mean.dlpwl)
  j.hat<-j.hat/m
  
  
  return(j.hat)
}
  


####################### For inverted brown-resnick model

# pairwise censored likelihood for brown-resniek
invbrown.pwl.censored<-function(param, y, x1, x2,  nsites, ntimes, delta=Inf,threshold) 
{
  val <- 0
  sigma<-param[1] 
 theta<-param[2]	
  
  tmp<-.C("invbrownpwl_censored",y = as.double(y), x1 = as.double(x1), x2 = as.double(x2),  
          n = as.integer(nsites), ntimes = as.integer(ntimes),sigma = as.double(sigma), theta = as.double(theta),
          val = as.double(val), delta = as.double(delta), threshold = as.double(threshold), NAOK = TRUE, DUP=TRUE)
  
  return(-tmp$val)	
}


# Fit brown-resniek using pairwise censored likelihood 
invbrown.fit.censored<-function(data, coords, param = NULL, delta=Inf,
                                prob =NULL,method ="nlminb",lower=NULL,upper=NULL,control=list(kkt=TRUE),std=TRUE, m=1) {
  
  # data matrix of data each row contains
  # coords matrix of coordinates
  # m number of observations per block
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  nsites <- dim(data)[2]
  ntimes <- dim(data)[1]
  
  y <- as.numeric(t(data))
  
  threshold<--1/log(prob)	
  res<-list()
  res$nsites<-nsites
  res$ntimes<-ntimes		
  res$start<-param
  res$prob<-prob
  res$threshold<-threshold	
  parnames<-c("sigma","theta")		
  if (std){ 
    print("maximum pairwise likelihood with std: first step")
  }
  else {
    print("maximum pairwise likelihood without std")	
  }
  
  a <-optimx(par=param, invbrown.pwl.censored, y = y, x1 = x1, x2 = x2, nsites = nsites,
             ntimes =ntimes,delta = delta,threshold=threshold,method =method,lower=lower,upper=upper,control=control)
  res$negplik <- a$value	
  res$hessian<-attr(a,"details")[1,]$nhatend
  res$grad<-attr(a,"details")[1,]$ngatend
  res$param<-as.numeric(coef(a))
  h.hat<-res$hessian
  res$delta<-delta
  
  if (std) 
  {
    
    print("maximum pairwise likelihood with std:  standard errors")	
    
    inv.h.hat<-ginv(h.hat)
    res$h.hat<-h.hat
    j.hat<-j.hat.invbrown.censored.func(res$param,data=data,coords=coords, delta=delta, m=m, threshold=threshold)
    res$j.hat<-j.hat
    pen<-sum(diag(j.hat%*%inv.h.hat))
    res$pen<-pen
    var.hat<-inv.h.hat%*%j.hat%*%inv.h.hat
    res$var.hat<-var.hat	
    std.err<-sqrt(diag(var.hat))
    res$std.err<-std.err
    res$clic<-2*(res$negplik+pen)		
  }
  
  return(res)	
}



j.hat.invbrown.censored.func<-function(param,data,coords, delta, m=1, threshold) {
  # m number of observations per block
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  ntimes <- dim(data)[1]
  nsites <- dim(data)[2]
  j.hat<-matrix(0,length(param),length(param))
  nblocks<-seq(1,ntimes-m+1)
  mean.dlpwl<-rep(0,length(param))
  for (i in nblocks)
  {
    y <- as.numeric(t(data[i:(i+m-1),]))
    
    
    dlpwl<-grad(invbrown.pwl.censored,param,y = y, x1 = x1, x2 = x2, nsites = nsites, ntimes = m,
                delta = delta,threshold=threshold)
    
    mean.dlpwl<-mean.dlpwl+dlpwl/length(nblocks)
    j.hat<-j.hat+dlpwl%*%t(dlpwl)/length(nblocks)
  }
  
  j.hat<-j.hat-mean.dlpwl%*%t(mean.dlpwl)
  j.hat<-j.hat/m
  
  
  return(j.hat)
}
  
  
  

####################### For mixture models with constant a
###################### For mixture between TEG and inverted brown-resniek 

TEGbrown.pwl.censored<-function(param, y, x1, x2,  nsites, ntimes, delta=Inf, corrmodel = "exponential",threshold=NULL) 
{
  # pairwise cencord likelihood for the mixture between TEG and inverted brown-resniek
  
  val <- 0
  mix<-param[1]
  theta1<-param[2]
  radius<-param[3]
  sigma<-param[4]
  theta2<-param[5]
  cr<-check.corrmodel(corrmodel)
  
  tmp<-.C("TEGbrownpwl_censored",y = as.double(y), x1 = as.double(x1), x2 = as.double(x2),  
          n = as.integer(nsites), ntimes = as.integer(ntimes),theta1 = as.double(theta1), radius = as.double(radius), sigma = as.double(sigma),
          theta2 = as.double(theta2), val = as.double(val), delta = as.double(delta), corrmodel = as.integer(cr), mix = as.double(mix),
          threshold=as.double(threshold),NAOK = TRUE, DUP=TRUE)
  
  return(-tmp$val)
  
}


TEGbrown.fit.censored<-function(data, coords, param = NULL, delta=Inf, corrmodel='exponential' ,prob =NULL,method ="nlminb",
                                lower=NULL,upper=NULL,control = list(kkt = TRUE),std=TRUE,m=1 ) {
  # data matrix of data each row contains
  # coords matrix of coordinates
  # corrmodel correlation model
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  nsites <- dim(data)[2]
  ntimes <- dim(data)[1]
  
  y <- as.numeric(t(data))
  
  threshold<- -1/log(prob)
  res<-list()
  res$nsites<-nsites
  res$ntimes<-ntimes
  res$prob<-prob
  res$start<-param
  res$threshold<-threshold
  
  parnames<-c("mix","theta1","radius","sigma","theta2")	
  
  if (std) {
    print("maximum pairwise likelihood with std: first step")	
  }
  else {
    print("maximum pairwise likelihood without std")
    
  }
  
  a <-optimx(par=param, TEGbrown.pwl.censored, y = y, x1 = x1, x2 = x2, nsites = nsites, ntimes = ntimes, delta = delta, 
             corrmodel = corrmodel,threshold=threshold,method = method,lower=lower,upper=upper,control = control)
  
  
  res$negplik <- a$value	
  res$hessian<-attr(a,"details")[1,]$nhatend
  res$grad<-attr(a,"details")[1,]$ngatend
  res$param<-as.numeric(coef(a))
  h.hat<-res$hessian
  res$corrmodel<-corrmodel
  res$delta<-delta
  
  if (std) 
  {
    
    print("maximum pairwise likelihood with std: standard errors")	
    
    inv.h.hat<-ginv(h.hat)
    res$h.hat<-h.hat
    param<-res$param
    j.hat<-j.hat.TEGbrown.censored.func(param,data=data,coords=coords, delta=delta, m=m,corrmodel=corrmodel,threshold=threshold)
    res$j.hat<-j.hat
    pen<-sum(diag(j.hat%*%inv.h.hat))
    res$pen<-pen
    var.hat<-inv.h.hat%*%j.hat%*%inv.h.hat
    res$var.hat<-var.hat	
    std.err<-sqrt(diag(var.hat))
    res$std.err<-std.err
    res$clic<-2*(res$negplik+pen)		
  }
  
  return(res)	 
}


j.hat.TEGbrown.censored.func<-function(param,data,coords, delta, m=1, corrmodel,threshold) {
  # m number of observations per block
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  ntimes <- dim(data)[1]
  nsites <- dim(data)[2]
  j.hat<-matrix(0,length(param),length(param))
  nblocks<-seq(1,ntimes-m+1)
  mean.dlpwl<-rep(0,length(param))
  for (i in nblocks)
  {
    y <- as.numeric(t(data[i:(i+m-1),]))
    
    dlpwl<-grad(TEGbrown.pwl.censored,param,y = y, x1 = x1, x2 = x2, nsites = nsites, ntimes = m,
                delta = delta, corrmodel = corrmodel,threshold=threshold)
    
    mean.dlpwl<-mean.dlpwl+dlpwl/length(nblocks)
    j.hat<-j.hat+dlpwl%*%t(dlpwl)/length(nblocks)
  }
  
  j.hat<-j.hat-mean.dlpwl%*%t(mean.dlpwl)
  
  j.hat<-j.hat/m
  
  
  return(j.hat)
}


###################### For mixture between TRUNCATED WITH inverted smith
sgj.pwl.censored<-function(param, y, x1, x2,  nsites, ntimes, delta=Inf, corrmodel = "exponential",threshold=NULL) 
{
  # pairwise likelihood for TRUNCATED WITH inverted smith
  
  val <- 0
  mix<-param[1]
  theta<-param[2]
  radius<-param[3]
  sigma<- param[4]
  cr1<-check.corrmodel(corrmodel)
  
  
  tmp<-.C("sgjpwl_censored",y = as.double(y), x1 = as.double(x1), x2 = as.double(x2),  
          n = as.integer(nsites), ntimes = as.integer(ntimes),theta = as.double(theta), radius = as.double(radius), sigma = as.double(sigma),
          val = as.double(val), delta = as.double(delta), corrmodel = as.integer(cr1), mix = as.double(mix),
          threshold=as.double(threshold),NAOK = TRUE, DUP=TRUE)
  
  
  return(-tmp$val)
  
}
##################################################################################################
sgj.fit.censored<-function(data, coords, param = NULL, delta=Inf, corrmodel='exponential' ,prob =NULL,method ="nlminb",
                               lower=NULL,upper=NULL,control = list(kkt = TRUE),std=TRUE,m=1) {

  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  nsites <- dim(data)[2]
  ntimes <- dim(data)[1]
 
  y <- as.numeric(t(data))
  
  threshold<- -1/log(prob)
  res<-list()
  res$nsites<-nsites
  res$ntimes<-ntimes
  res$prob<-prob
  res$start<-param
  res$threshold<-threshold
	
  parnames<-c("theta","radius","sigma")	
  if (std) {
    print("maximum pairwise likelihood with std: first step")	
  }
  else {
    print("maximum pairwise likelihood without std")
    
  }
 
  a <-optimx(par=param,sgj.pwl.censored, y = y, x1 = x1, x2 = x2, nsites = nsites, ntimes = ntimes, delta = delta, 
             corrmodel = corrmodel,threshold=threshold,method = method,lower=lower,upper=upper,control = control)
  
  
  res$negplik <- a$value	
  res$hessian<-attr(a,"details")[1,]$nhatend
  res$grad<-attr(a,"details")[1,]$ngatend
  res$param<-as.numeric(coef(a))
  h.hat<-res$hessian
  res$corrmodel<-corrmodel
  res$delta<-delta
  
  if (std) 
  {
    
    print("maximum pairwise likelihood with std: standard errors")	
    
    inv.h.hat<-ginv(h.hat)
    res$h.hat<-h.hat
    param<-res$param
    j.hat<-j.hat.sgj.censored.func(param,data=data,coords=coords, delta=delta, m=m,corrmodel=corrmodel,threshold=threshold)
    res$j.hat<-j.hat
    pen<-sum(diag(j.hat%*%inv.h.hat))
    res$pen<-pen
    var.hat<-inv.h.hat%*%j.hat%*%inv.h.hat
    res$var.hat<-var.hat	
    std.err<-sqrt(diag(var.hat))
    res$std.err<-std.err
    res$clic<-2*(res$negplik+pen)		
  }
  
  return(res)	 
}		
  
##########################################
j.hat.sgj.censored.func<-function(param,data,coords, delta, m=1, corrmodel,threshold) {
  # m number of observations per block
  x1 <- as.numeric(coords[,1])
  x2 <- as.numeric(coords[,2])
  ntimes <- dim(data)[1]
  nsites <- dim(data)[2]
  j.hat<-matrix(0,length(param),length(param))
  nblocks<-seq(1,ntimes-m+1)
  mean.dlpwl<-rep(0,length(param))
  for (i in nblocks)
  {
    y <- as.numeric(t(data[i:(i+m-1),]))
    
    dlpwl<-grad(sgj.pwl.censored,param,y = y, x1 = x1, x2 = x2, nsites = nsites, ntimes = m,
                delta = delta, corrmodel = corrmodel,threshold=threshold)
    
    mean.dlpwl<-mean.dlpwl+dlpwl/length(nblocks)
    j.hat<-j.hat+dlpwl%*%t(dlpwl)/length(nblocks)
  }
  
  j.hat<-j.hat-mean.dlpwl%*%t(mean.dlpwl)
  
  j.hat<-j.hat/m
  
  
  return(j.hat)
}


# estimation the parameters for three models 
# Fitting the proposed Madogramme by Non-linear least square
#########################################################
mado.fit.models <- function(par,data,coord,Model='mixture',structure='SC',miss=TRUE)
{
  emp.mado  <- function(data, coord) {
    Fr <- pfrechet(data, shape = 1)
    loc.all <- coord
    pair.w <- combn(ncol(Fr), 2)
    
    n_pairs <- ncol(pair.w)
    n_obs <- nrow(Fr)
    abs.F <- array(NA, dim = c(n_obs, n_pairs, 1))
    
    dist1 <- as.matrix(dist(loc.all, method = "euclidean"))
    dist_vec <- sapply(1:n_pairs, function(j) dist1[pair.w[1,j], pair.w[2,j]])
    
    for (i in 1:n_pairs) {
      v <- 0.5 * abs(Fr[, pair.w[1, i]] - Fr[, pair.w[2, i]])
      abs.F[, i, 1] <- v  # optionally: pmin(v, 0.1666667)
    }
    Y<-  abs.F[, , 1]
    return(list(dist=dist_vec,Y=Y))
  }
  dep <- emp.mado(data,coord)
  dist <- dep$dist
  big.M  <- dep$Y
  in.obsr <- dim(big.M)[1]
  n.pair <- dim(big.M)[2]
  
  if(Model=='mixture'){
    if      (structure=='SM')
    {
      lower<-c(0,min(dist),min(dist),min(dist))
      upper<-c(1,Inf,Inf,Inf)
      
    }
    else if (structure=='SC' |structure=='BR' ) 
    {
      lower<-c(0,min(dist),min(dist),min(dist),min(dist))
      upper<-c(1,Inf,Inf,Inf,Inf)
    }
  }
  
  if(Model=='max'| Model=='imax'){
    if (structure=='SM') 
    {
      lower<-min(dist)
      upper<-max(dist)
    }
    else if (structure=='SC'|structure=='BR')
    {
      
      lower<-c(min(dist),min(dist))
      upper<-c(max(dist),max(dist))
      
    }
  }
  
  
  myfun <- function(par)
  {
    #par <-par.tr(par) 
    if (Model=='mixture'){
      x1<-par[1]
      x2<-par[2]
      x3<-par[3]
      
      rho.x <- exp(-(dist/x2))
      alpha.x <- overlarea(dist,x3)#(1-(0.5*dist/x3))*(dist<(2*x3))#
      ext.x <- (2-(alpha.x*(1-sqrt(0.5*(1-rho.x)))))
      
      if(structure=="SM")
      {
        x4<-par[4]
        ext.y <-(2*pnorm((dist)/(2*sqrt(x4^2))))
      }
      else if(structure=="SC")
      {
        x4<-par[4]    
        x5<-par[5]
        rho.y <- exp(-(dist/x4))#1/(1+(dist/x4)^2)#  
        alpha.y <- overlarea(dist,x5)#(1-(0.5*dist/x4))*(dist<(2*x4))#
        ext.y <- (2-(alpha.y*(1-sqrt(0.5*(1-rho.y)))))      
        
      }
      else if(structure=="BR")
      {
        x4<-par[4]
        x5<-par[5]
        vario<- (x5)*sqrt(2*(1-exp(-(dist/x4))))
        ext.y <- 2*pnorm(vario/2)      
        
      }
      
      mado <-sum((as.matrix((1/in.obsr^2)*(apply(big.M, 1,function(y) (y-(((x1*(ext.x-1))/((x1*(ext.x-1))+2))-((x1*ext.x-1)/(2*x1*ext.x+2))-((ext.y/(x1*ext.x + (1-x1)*ext.y + 1))*beta(((x1*ext.x+1)/(1-x1)),ext.y))))^2)))))
    }
    if(Model=='max')
    {
      if(structure=="SM")
      {
        x<-par[1]
        vario<- dist/sqrt(x)
        ext <- 2*pnorm(vario/2)  
      }
      if(structure=="BR")
      {
        x1<-par[1]
        x2<-par[2]
        vario<- (x2)*sqrt(2*(1-exp(-(dist/x1)^2)))
        ext <- 2*pnorm(vario/2)  
      }
      if(structure=="SC")
      {
        x<-par[1]
        y<-par[2]
        
        rho <- exp(-(dist/x)^2)
        alpha <- overlarea(dist,y)#(1-(0.5*dist/x4))*(dist<(2*x4))#
        ext <- (2-(alpha*(1-sqrt(0.5*(1-rho)))))
      }
      mado <-sum((as.matrix((1/in.obsr)*(apply(big.M, 1,function(y) (y-((ext-1)/(2*(ext+1))))^2)))))
    }
    if(Model=='imax')
    {
      if(structure=="SM")
      {
        x<-par[1]
        vario<- dist/sqrt(x^2)
        eta <- 1/(2*pnorm(vario/2))  
      }
      if(structure=="BR")
      {
        x1<-par[1]
        x2<-par[2]
        vario<- (x2)*sqrt(2*(1-exp(-(dist/x1))))
        eta <- 1/(2*pnorm(vario/2))  
      }
      if(structure=="SC")
      {
        x<-par[1]
        y<-par[2]
        
        rho <- exp(-(dist/x))
        alpha <- overlarea(dist,y)#(1-(0.5*dist/x4))*(dist<(2*x4))#
        eta <- 1/(2-(alpha*(1-sqrt(0.5*(1-rho)))))
      }
      mado <-sum((as.matrix((1/in.obsr^2)*(apply(big.M, 1,function(y) (y-((1-eta)/(2*(1+eta))))^2)))))
      
    }
    return(na.omit(mado))
  }
  one<- optimx::optimx(par, myfun,method ="nlminb",lower =lower ,upper =upper,control = list(kkt = TRUE))
  return(one)
}
########################################################################
# Begin the computation the simulation F-madogram #
#======================================================================================

# Compute Empirical F-monogram
def.mad <- function(data) {
  Fr <- pfrechet(data$data, shape=1)
  loc.all <- data$loc.all
  pair.w <- combn(ncol(Fr), 2)
  
  n_pairs <- ncol(pair.w)
  n_obs <- nrow(Fr)
  
  dist  <- array(NA, dim = c(1, n_pairs, 1))
  abs.F <- array(NA, dim = c(n_obs, n_pairs, 1))
  
  dist1 <- as.matrix(dist(loc.all, method = "euclidean"))
  dist[1,,1] <- sapply(1:n_pairs, function(j) dist1[pair.w[1,j], pair.w[2,j]])
  
  for (i in 1:n_pairs) {
    v <- 0.5 * abs(Fr[, pair.w[1, i]] - Fr[, pair.w[2, i]])
    abs.F[, i, 1] <- v  # optionally: pmin(v, 0.1666667)
  }
  
  return(list(dist = dist, abs.fr.dif = abs.F))
}

mado.fit.Simu <- function(par,data,AI.model="SM"){
  AI.model=AI.model
  if(AI.model=="SM")
  {
    n.par<- 4
    one <- matrix(NA,(dim(def.mad(data)$abs.fr.dif)[3]),n.par)
  }
  if(AI.model=="SC")
  {
    n.par<- 5
    one <- matrix(NA,(dim(def.mad(data)$abs.fr.dif)[3]),n.par)
  }
  #========================================================
  
  itr<-0
  for(i in 1:(dim(def.mad(data)$abs.fr.dif)[3]))
  {
    dist <- (def.mad(data)$dist)[,,i]
    big.M  <- (def.mad(data)$abs.fr.dif)[,,i]
    
    in.obsr <-dim(big.M)[1]
    
    
    myfun <- function(par)
    {
      x1<-par[1]
      x2<-par[2]
      x3<-par[3]
      
      rho.x <- exp(-dist/x3)
      alpha.x <- overlarea(dist,x2)#(1-(0.5*dist/x2))*(dist<(2*x2))#
      ext.x <- (2-(alpha.x*(1-sqrt(0.5*(1-rho.x)))))
      
      if(AI.model=="SM")
      {
        x4<-par[4]
        ext.y <-(2*pnorm((dist)/(2*sqrt(x4))))
      }
      else if(AI.model=="SC")
      {
        x4<-par[4]
        x5<-par[5]
        rho.y <- exp(-dist/x5)
        alpha.y <- overlarea(dist,x4)#(1-(0.5*dist/x4))*(dist<(2*x4))#
        ext.y <- (2-(alpha.y*(1-sqrt(0.5*(1-rho.y)))))      
        
      }
      mado <-sum((as.matrix((1/in.obsr^2)*(apply(big.M, 1,function(y) (y-(((x1*(ext.x-1))/((x1*(ext.x-1))+2))-((x1*ext.x-1)/(2*x1*ext.x+2))-((ext.y/(x1*ext.x + (1-x1)*ext.y + 1))*beta(((x1*ext.x+1)/(1-x1)),ext.y))))^2)))))
      return(mado)
    }
    
    
    one[i,] <- optim(par, myfun,method="L-BFGS-B",lower =c(0,rep(0.01,n.par)) ,upper =c(1,rep(2,n.par)))$par
    
    itr<-itr+1
    print(itr)
    print(one[i,])
  }
  return(one)
}


  
