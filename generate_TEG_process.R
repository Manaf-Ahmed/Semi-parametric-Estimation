
system("R CMD SHLIB rdavghol.c")
dyn.load(paste("rdavghol",.Platform$dynlib.ext,sep=""))
###########################################################################
rdavghol<-function(coords, n.obs=1, corrmodel="exponential", 
                   range.param=stop("missing range parameter"), bound =3, radius=NULL)
{  
  
  ans<-rdavghol.C(coords, n.obs, corrmodel,range.param, bound, radius)   
  return(ans)
}

rdavghol.C<-function(coords, n.obs=1, corrmodel="exponential", 
                     range.param=stop("missing range parameter"), bound =3, radius=NULL)
{  
  # This function generates a random field 
  # on W^+\times (xlims[1],xlims[2])\times(ylims[1],ylims[2]) for the Davison Gholam Razee model
  #    coords coordinates
  #    cov.model the covariance model
  #    range.param  parameter for the range
  #    radius  radius of the ball
  #    bound for the extremal Gaussian process we choose (bound) C=3.
  #    xlims   Limits of the area in the X direction. Defaults to [0,1].
  #    ylims  Limits of the area in the Y direction. Defaults to [0,1].
  #    n.obs    Number of simulations. Defaults to 1.
  #    dataZ simulation
  
  cov.mod <- switch(corrmodel, "whitmat" = 1, "cauchy" = 2,
                    "powexp" = 3, "bessel" = 4, "exponential" = 5,"spherical" = 6)
  
  const<-sqrt(2)/(sqrt(pi)*radius^2) 
  uBound<-bound*const
  # we specify the window W i.e. we compute the convex hull of the points
  w<-convexhull.xy(coords)
  n.site<-nrow(coords)
  
  q<-ncol(coords)
  
  limitsW<-c(min(coords[,1])-radius,max(coords[,1])+radius,min(coords[,2])-radius, max(coords[,2])+radius)
  scalefact<-sqrt(2)/(sqrt(pi)*radius^2) 
  sigma2<-1
  nugget<-0
  smooth<-0.5
  ans <- rep(-1e10, n.obs * n.site)
  ans<-.C("rdavghol",coord = as.double(coords), nObs=as.integer(n.obs), 
          nSite=as.integer(n.site),  
          dim = as.integer(q), covmod=as.integer(cov.mod), limitsW =as.double(limitsW), radius= as.double(radius), range = as.double(range.param),
          scalefact = as.double(scalefact), smooth = as.double(smooth), uBound =
            as.double(uBound), ans = as.double(ans), NAOK = TRUE, DUP=TRUE)$ans
  ans <- matrix(ans, nrow = n.obs, ncol = n.site)
  return(list(data=ans,coords=coords))  
}  

# Generate Max-mtixture process
Simu.data<- function(par,n.site,n.obs,R.area,simu.model,corrmodel.x,corrmodel.y)
{
  
  m.p<-par[1]
  range.param.x<-par[2]
  radius.x<-par[3]
  
  
  loc.all <-cbind(runif(n.site,0,R.area),runif(n.site,0,R.area))
  X  <-rdavghol(loc.all,n.obs=n.obs,corrmodel=corrmodel.x ,range.param=range.param.x, radius=radius.x)$data
  
  if(simu.model=="simu.MM1")
  {
    range.param.y<-par[4]
    Y1 <- rmaxstab(n=n.obs, coord=loc.all, cov.mod ="gauss", 
                   cov11 = range.param.y, cov12 = 0, cov22 = range.param.y,grid = FALSE)
    
  }
  if (simu.model=="simu.MM2")
  {
    radius.y<-par[4]
    range.param.y<-par[5]
    
    Y1  <-rdavghol(loc.all, n.obs=n.obs,corrmodel=corrmodel.y,
                   range.param=range.param.y, 
                   radius=radius.y)$data 
    
  }
  Y <- -1/log(1-exp(-1/Y1))
  Z <- pmax(m.p*X,(1-m.p)*Y)
  
  return(list(loc.all=loc.all,data=Z))
}
