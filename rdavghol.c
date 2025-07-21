#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
//#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <time.h>

#define MINF -1.0e15
#define EPS DBL_EPSILON


void distance(double *coord, int *nDim, int *nSite,
	      int *vec, double *dist){

  //This function computes either the euclidean distance or the
  //distance vector between each pair of locations
  const int nPair = *nSite * (*nSite - 1) / 2;
  int i, j, k, currentPair = 0;

  if (*vec){
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	for (k=0;k<*nDim;k++)
	  dist[k * nPair + currentPair] = coord[k * *nSite + j] -
	    coord[k * *nSite + i];

	currentPair++;
      }
    }
  }

  else{
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	dist[currentPair] = 0;
	for (k=0;k<*nDim;k++)
	  dist[currentPair] += (coord[i + k * *nSite] -	coord[j + k * *nSite]) *
	    (coord[i + k * *nSite] - coord[j + k * *nSite]);

	dist[currentPair] = sqrt(dist[currentPair]);
	currentPair++;
      }
    }
  }
}





int getCurrentPair(int site1, int site2, int nSite){
  //site1 has to be smaller than site2
  int ans = site1 * nSite - site1 * (site1 + 1) / 2 + site2 -
    site1 - 1;
  return ans;
}

void getSiteIndex(int currentPair, int nSite, int *site1, int *site2){
  int nFree = nSite - 2,
      cum = nSite - 2;

  *site1 = 0;

  while (currentPair > cum){
    *site1 = *site1 + 1;
    cum += nFree;
    nFree--;
  }

  *site2 = *site1 + currentPair - cum + nFree + 1;

  return;
}

double whittleMatern(double *dist, int n, double nugget, double sill, double range,
		     double smooth, double *rho){

  //This function computes the whittle-matern covariance function
  //between each pair of locations.
  //When ans != 0.0, the whittle-matern parameters are ill-defined.

  const double cst = sill * R_pow(2, 1 - smooth) / gammafn(smooth),
    irange = 1 / range;

  //Some preliminary steps: Valid points?
  if (smooth < EPS)
    return (1 - smooth + EPS) * (1 - smooth + EPS) * MINF;

  else if (smooth > 100)
    /* Not really required but larger smooth parameters are unlikely
       to occur */
    return (smooth - 99) * (smooth - 99) * MINF;

  if (range <= 0)
    return (1 - range) * (1 - range) * MINF;

  if (sill <= 0)
    return (1 - sill) * (1 - sill) * MINF;

  if (nugget < 0)
    return (1 - nugget) * (1 - nugget) * MINF;

  for (int i=0;i<n;i++){
    double cst2 = dist[i] * irange;

    if (cst2 == 0)
      rho[i] = sill + nugget;

    else
      rho[i] = cst * R_pow(cst2, smooth) * bessel_k(cst2, smooth, 1);
  }

  return 0.0;
}

double covExp(double *dist, int n, double nugget, double sill, double range, double *rho){

  //This function computes the exponential covariance function
  //between each pair of locations.
 

  const double cst = sill, irange = 1 / range;
  double cst2;

  for (int i=0;i<n;i++){
     cst2 = -dist[i] * irange;

    if (cst2 == 0)
      rho[i] = sill + nugget;

    else
      rho[i] = cst * exp(cst2);
     
  }

  return 0.0;
}


double cauchy(double *dist, int n, double nugget, double sill, double range,
	      double smooth, double *rho){

  //This function computes the cauchy covariance function between each
  //pair of locations.
  //When ans != 0.0, the cauchy parameters are ill-defined.

  const double irange2 = 1 / (range * range);

  //Some preliminary steps: Valid points?
  if (smooth < 0)
    return (1 - smooth) * (1 - smooth) * MINF;

  else if (smooth > 100)
    return (smooth - 99) * (smooth - 99) * MINF;

  if (range <= 0.0)
    return (1 - range) * (1 - range)* MINF;

  if (sill <= 0.0)
    return (1 - sill) * (1 - sill) * MINF;

  if (nugget < 0)
    return (1 - nugget) * (1 - nugget) * MINF;

  for (int i=0;i<n;i++){

    if (dist[i] == 0)
      rho[i] = nugget + sill;

    else
      rho[i] = sill * R_pow(1 + dist[i] * dist[i] * irange2, -smooth);
  }

  return 0.0;
}

double caugen(double *dist, int n, double nugget, double sill, double range,
	      double smooth, double smooth2, double *rho){

  /*This function computes the generalized cauchy covariance function
    between each pair of locations.  When ans != 0.0, the parameters
    are ill-defined. */

  const double irange = 1 / range, ratioSmooth = -smooth / smooth2;

  //Some preliminary steps: Valid points?
  if (smooth < 0)
    return (1 - smooth) * (1 - smooth) * MINF;

  /*else if (smooth1 > 500)
    return (smooth1 - 499) * (smooth1 - 499) * MINF; */

  if ((smooth2 > 2) || (smooth2 <= 0))
    return (1 - smooth2) * (1 - smooth2) * MINF;

  if (range <= 0)
    return (1 - range) * (1 - range)* MINF;

  if (sill <= 0)
    return (1 - sill) * (1 - sill) * MINF;

  if (nugget < 0)
    return (1 - nugget) * (1 - nugget) * MINF;

  for (int i=0;i<n;i++){
    if (dist[i] == 0)
      rho[i] = nugget + sill;

    else
      rho[i] = sill * R_pow(1 + R_pow(dist[i] * irange, smooth2),
			    ratioSmooth);
  }

  return 0.0;
}

double powerExp(double *dist, int n, double nugget, double sill, double range,
		double smooth, double *rho){

  //This function computes the powered exponential covariance function
  //between each pair of locations.
  //When ans != 0.0, the powered exponential parameters are ill-defined.

  const double irange = 1 / range;

  //Some preliminary steps: Valid points?
  if ((smooth < 0) || (smooth > 2))
    return (1 - smooth) * (1 - smooth) * MINF;

  if (range <= 0)
    return (1 - range) * (1 - range) * MINF;

  if (sill <= 0)
    return (1 - sill) * (1 - sill) * MINF;

  if (nugget < 0)
    return (1 - nugget) * (1 - nugget) * MINF;

  for (int i=0;i<n;i++){
    if (dist[i] == 0)
      rho[i] = nugget + sill;

    else
      rho[i] = sill * exp(-R_pow(dist[i] * irange, smooth));
  }

  return 0.0;
}

double bessel(double *dist, int n, int dim, double nugget, double sill,
	      double range, double smooth, double *rho){
  //This function computes the bessel covariance function
  //between each pair of locations.
  //When ans != 0.0, the powered exponential parameters are ill-defined.

  const double irange = 1 / range, cst = sill * R_pow(2, smooth) * gammafn(smooth + 1);

  //Some preliminary steps: Valid points?
  if (smooth < (0.5 * (dim - 2)))
    return (1 + 0.5 * (dim - 2) - smooth) * (1 + 0.5 * (dim - 2) - smooth) * MINF;

  /* else if (smooth > 100)
    //Require as bessel_j will be numerically undefined
    return (smooth - 99) * (smooth - 99) * MINF; */

  if (range <= 0)
    return (1 - range) * (1 - range) * MINF;

  if (sill <= 0)
    return (1 - sill) * (1 - sill) * MINF;

  if (nugget < 0)
    return (1 - nugget) * (1 - nugget) * MINF;

  for (int i=0;i<n;i++){
    double cst2 = dist[i] * irange;

    if (cst2 == 0)
      rho[i] = nugget + sill;

    else if (cst2 <= 1e5)
      rho[i] = cst * R_pow(cst2, -smooth) * bessel_j(cst2, smooth);

    else
      // approximation of the besselJ function for large x
      rho[i] = cst * R_pow(cst2, -smooth) * M_SQRT_2dPI *
	cos(cst2 - smooth * M_PI_2 - M_PI_4);

    /*if (!R_FINITE(rho[i]))
      return MINF;*/
  }

  return 0.0;
}

double mahalDistFct(double *distVec, int n, double *cov11,
		    double *cov12, double *cov22, double *mahal){
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 2D
  //When ans != 0.0, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.

  const double det = *cov11 * *cov22 - *cov12 * *cov12,
    idet = 1 / det;

  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if (*cov11 <= 0)
    return (1 - *cov11) * (1 - *cov11) * MINF;

  if (*cov22 <= 0)
    return (1 - *cov22) * (1 - *cov22) * MINF;

  if (det <= 0)
    return (1 - det) * (1 - det) * MINF;

  for (int i=0;i<n;i++)
    mahal[i] = sqrt((*cov11 * distVec[n + i] * distVec[n + i] -
		     2 * *cov12 * distVec[i] * distVec[n + i] +
		     *cov22 * distVec[i] * distVec[i]) * idet);

  return 0.0;
}

double covSpherical(double *dist, int n, double nugget, double sill, double range, double *rho){

  //This function computes the spherical covariance function
  //between each pair of locations.
 

  const double cst = sill, irange = 1 / range;
  double cst2;

  for (int i=0;i<n;i++){
	if (dist[i] < range) {
      cst2 = dist[i] * irange;

    if (cst2 == 0)
      rho[i] = sill + nugget;

    else
      rho[i] = cst * (1 - (1.5 * cst2) + 0.5 * R_pow(cst2,3));
    }
    else 
    {
		rho[i] = 0.0;
	} 
  }

  return 0.0;
}




void buildcovmat(int *nSite, int *covmod, double *coord, int *dim,
		 double *nugget, double *sill, double *range,
		 double *smooth, double *covMat){

  int currentPair, nPairs, effnSite = *nSite, zero = 0;
  const double one = 1.0, dzero = 0.0;
  double flag = 0;

  nPairs = effnSite * (effnSite - 1) / 2;

  double *dist = malloc(nPairs * sizeof(double)),
    *rho = malloc(nPairs * sizeof(double));

  distance(coord, dim, nSite, &zero, dist);

  switch (*covmod){
  case 1:
    flag = whittleMatern(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 2:
    flag = cauchy(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 3:
    flag = powerExp(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 4:
    flag = bessel(dist, nPairs, *dim, dzero, one, *range, *smooth, rho);
    break;
  case 5:
    flag = covExp(dist, nPairs, dzero, one, *range, rho);
    break;  
  case 6:
    flag = covSpherical(dist, nPairs, dzero, one, *range, rho);
    break;  

  }
  
  if (flag != 0.0)
    error("The covariance parameters seem to be ill-defined. Please check\n");
  
  //Fill the non-diagonal elements of the covariance matrix
  currentPair = -1;
  for (int i = 0; i < (effnSite-1); i++){
    for (int j = i + 1; j < effnSite; j++){
      currentPair++;
      covMat[effnSite * i + j] = covMat[effnSite * j + i] =*sill * rho[currentPair];
    }
  }
  
  //Fill the diagonal elements of the covariance matrix
    
  for (int i = 0; i < effnSite; i++)
      covMat[i * (effnSite + 1)] = *sill + *nugget;


  free(dist); free(rho);
  return;
}





void rdavghol(double *coord, int *nObs, int *nSite, int *dim,
		 int *covmod, double *limitsW, double *radius,  
		 double *range,double *scalefact, double *smooth, double *uBound,
		 double *ans){
  
  int i, j, k, neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double area,  d, nugget = 0, runifx, runify, sigma2 = 1, one = 1, zero = 0;

  
    neffSite = *nSite;
    lagj = *nObs;
  
  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double));

  buildcovmat(nSite, covmod, coord, dim, &nugget, &sigma2, range,
	      smooth, covmat);
  
  /* Compute the Cholesky decomposition of the covariance matrix */
  int info = 0;
  //F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info);
  F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info, 1);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();
  
  for (i=*nObs;i--;){
    double poisson = 0;
    int nKO = neffSite;
    
    while (nKO) {
      /* The stopping rule is reached when nKO = 0 i.e. when each site
	 satisfies the condition in Eq. (8) of Schlather (2002) */
      poisson += exp_rand();
      
	  /* We simulate uniformly one realisation of  on  the dilated window W */
 	  runifx = limitsW[0]+(limitsW[1]-limitsW[0])*unif_rand();
 	  runify = limitsW[2]+(limitsW[3]-limitsW[2])*unif_rand();
      /* We simulate one realisation of a Gaussian random field with
	 the required covariance function */
	 
      for (j=neffSite;j--;)
			gp[j] = norm_rand();
      
      //F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt);
      F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt, 1, 1, 1);
     
     for (j=neffSite;j--;) {
		 d =0.0;
		 
	     d = (coord[j]-runifx) *(coord[j] - runifx)
	          + (coord[j + neffSite] -	runify)*(coord[j+neffSite] - runify);
	      d = sqrt(d) ;     
		  gp[j]= (*scalefact)*gp[j] * (gp[j] > 0)*(d <= *radius);
		  //Rprintf("%f \n",gp[j]);
	      }
	
      nKO = neffSite;
      
      for (j=neffSite;j--;){
			ans[j * lagj + i * lagi] = fmax2(gp[j]/poisson,
				      ans[j * lagj + i * lagi]);
				      
			
			nKO -= (*uBound/poisson <= ans[j * lagj + i * lagi]);
      }
    }
  }
  
  PutRNGstate();
   
  area = (limitsW[1]-limitsW[0])*(limitsW[3]-limitsW[2]);
  for (i=*nObs * neffSite;i--;)
    ans[i] = area*ans[i];

  
  free(covmat); free(gp);

  return;
}

