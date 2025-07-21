# estimation the parameters for three models 
#==================================================================================
# Load required packages
load_required_packages <- function() {
  pkgs <- c(
    "SpatialExtremes", "geoR", "mgcv", "polyclip", "spatstat", "numDeriv", "MASS",
    "mnormt", "ucminf", "optimx", "mev", "POT", "scales", "cubature",
    "copula", "VGAM", "automap", "gstat", "R.matlab", "extRemes", "matrixStats",
    "reshape","fossil", "ggplot2","dplyr","tidyr","cowplot"
  )
  for (pkg in unique(pkgs)) {
    if (!require(pkg, character.only = TRUE)) {
      warning(sprintf("Package '%s' is not installed.", pkg))
    }
  }
}
load_required_packages()

# Load required sources 
source("generate_TEG_process.R") # Estmation function MLE and NLS of Madograme
source("estimation_code.R")
source("box_bar_density_plots.R")


#====================================================================================
# Begin the computation the simulation F-madogram #
#======================================================================================
n.itr<-100
n.obs<-1000
n.site <- 50
prob<- 0.9

res.ls<-res.lik<- matrix(NA,n.itr,4)

for (i in 1:n.itr){
  dat <- Simu.data(c(0.25,0.20,0.25,.6),n.site,n.obs,1,simu.model="simu.MM1",
                   corrmodel.x="exponential")
  data <- dat$data
  loc<- dat$loc.all
  #   Likelihood estimation of  the max-mixture of the TEG model 
      #with the inverse Smith model
 est <- sgj.fit.censored(data, loc, param = c(0.2,0.2,0.2,0.4), delta=Inf, 
                         corrmodel='exponential',
                         prob =prob, method ="nlminb",
                         lower=c(0,0.01,0.01,0.01),
                          upper=c(1,Inf,Inf,Inf),
                         control = list(kkt =TRUE,maxit=1000),
                         std=F, m=1)
 ######## Fit models using a simplified version of the F-madogram
 
 # Structure options (dependence structure for the spatial model):
 #   AI.model = 'SM' – refers to the max-mixture of the TEG model with the inverse Smith model
 #   AI.model = 'SC' – refers to the max-mixture of the TEG model with the inverse TEG model
  mst<- mado.fit.Simu(par=c(.1,.1,.15,0.4),dat,AI.model="SM")
  res.ls[i,] <- as.numeric(mst)
  res.lik[i,] <- est$param
  print(i)
}


true_params <- c(a = 0.25, rx = 0.25, theta_X = 0.2, sigma_Y = 0.6)
# Subtract true parameters from each column of estimates
res.diff.ls <- sweep(res.ls, 2, true_params, FUN = "-")
res.diff.likl <- sweep(res.lik, 2, true_params, FUN = "-")
boxplot_mado<-make_boxplot (res.diff.ls, c("a", "rx", "theta_X", "sigma_Y"), 
                            expression(hat(psi)[M] - psi))
boxplot_likl<- make_boxplot(rres.diff.likl, c("a", "rx", "theta_X", "sigma_Y"), 
                            expression(hat(psi)[M] - psi))

rmse_values_ls <- matrix(sqrt(colMeans((res.ls - matrix(rep(true_params,
                                              each = nrow(res.ls)), ncol = 4))^2)), 
                         nrow = 1)
                         
rmse_values_likl <- matrix(sqrt(colMeans((res.lik - matrix(rep(true_params,
                                                each = nrow(res.lik)), ncol = 4))^2)),
                           nrow = 1)
colnames(rmse_values_ls) <- c("a", "rx", "theta_X", "sigma_Y")
colnames(rmse_values_likl) <- c("a", "rx", "theta_X", "sigma_Y")

barplot_mado <- make_barplot(rmse_values_ls,
                             colnames_vec = c("a", "rx", "theta_X", "sigma_Y"),  
                             ylab_expr = expression(hat(psi)[M] - psi), 
                             ymax = 0.1)

barplot_likl <- make_barplot(rmse_values_likl, 
                             colnames_vec = c("a", "rx", "theta_X", "sigma_Y"),  
                             ylab_expr = expression(hat(psi)[L] - psi), 
                             ymax = 0.1)

densityplot_mado<- make_densityplot(res.diff.ls, c("a", "rx", "theta_X", "sigma_Y"), expression(hat(psi)[M] - psi))
#====================================================================================
# Begin fitting the real dataset
#======================================================================================

mydata <- load_data()
coord <- load_coord()
COORD<- lonlat.to.planar(coord, miles =FALSE)
COORD<- cbind(rescale(COORD[,1]),rescale(COORD[,2]))
X=data2frechet(mydata)

######## fit  modeles bt F-Madogramm 
# Model options:
#   'mixture'  – refers to the Max-Mixture model
#   'max'      – refers to Max-Stable models
#   'imax'     – refers to Inverted Max-Stable models
#
# Structure options (dependence structure for the spatial model):
#   'SM' – refers to the Smith model
#   'SC' – refers to the Schlather model
#   'BR' – refers to the Brown–Resnick model

fit_mado<- mado.fit.models(par=c(.2,50,50,50,100),mydata,COORD,
                            Model='mixture',structure='SM',miss=TRUE)

######## fit smith model

f_smith=smith.fit.censored(X, COORD, param = 300, delta=Inf,
                           prob =prob,method ="nlminb",lower=0.01,
                           upper=Inf,control = list(kkt =TRUE,maxit=1000),std=TRUE, m=1)


######## fit inverted smith model

f_invsmith=invsmith.fit.censored(X, COORD, param = 300, delta=Inf,
                                 prob =prob,method ="nlminb",lower=0.01,
                                 upper=Inf,control = list(kkt =TRUE,maxit=1000),std=TRUE, m=1) 



######### fit TEG model
f_TEG=davghol.fit.censored(X, COORD, param = c(300,500), delta=Inf, corrmodel='exponential' ,
                           prob =prob,method ="nlminb",lower=c(0.01,0.01),
                           upper=c(Inf,Inf),control = list(kkt =TRUE,maxit=1000),std=TRUE, m=1)


######## fit Brown-resiek model

f_BR=brown.fit.censored(X, COORD, param = c(300,100), delta=Inf,
                        prob =prob,method ="nlminb",lower=c(1,1),
                        upper=c(Inf,Inf),control = list(kkt =TRUE,maxit=1000),std=TRUE, m=1)


######## fit inverted Brown-resiek model

f_invBR=invbrown.fit.censored(X, COORD, param = c(10,100), delta=Inf,
                              prob =prob,method ="nlminb",lower=c(0.01,0.01),
                              upper=c(Inf,Inf),control = list(kkt =TRUE,maxit=1000),std=TRUE, m=1) 


##################### fit max-mixture between TEG and inverted BR

f_TEG_invBR= TEGbrown.fit.censored(X, COORD, param = c(0.25,100,300,400,0.5), delta=Inf, corrmodel='exponential',
                                   prob =prob, method ="nlminb",
                                   lower=c(0.01,0.01,0.01,0.01,0.01),
                                   upper=c(0.99,Inf,Inf,Inf,Inf),
                                   control = list(kkt =TRUE,maxit=1000),
                                   std=TRUE, m=1)

##################### fit max-mixture between TEG and inverted smith

f_TEG_invSM= sgj.fit.censored(X, COORD, param = c(0.25,100,300,400), delta=Inf, corrmodel='exponential',
                              prob =prob, method ="nlminb",
                              lower=c(0.01,0.01,0.01,0.01),
                              upper=c(0.99,Inf,Inf,Inf),
                              control = list(kkt =TRUE,maxit=1000),
                              std=TRUE, m=1)

