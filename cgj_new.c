#include <stdio.h>
#include <stdlib.h> /* for malloc & free */
#include "f2c.h"
#include <math.h> 
#include <R.h>
#include <Rinternals.h>

#include <R_ext/Arith.h>	/* NA handling */
#include <Rmath.h>	
#include <R_ext/Random.h>	/* ..RNGstate */
#include <R_ext/Applic.h>	/* NA handling */
#define SMALL 1e-320
#define MINF -1.0e15

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

double F77_NAME(mvbvn)(double *lower, 
			double *upper, 
			int* infin, 
			double *correl);

double dist(double a1, double b1, double a2, double b2 )
/*
evaluate the euclidean distance
*/
{

	double val, d1 , d2;
	
	d1 = a1 - a2;
	d2 = b1 - b2;
	val = R_pow_di(d1, 2) + R_pow_di(d2, 2);
	val = sqrt(val);   

	return(val);
}

void  dispointset(double *x1, double *x2, double *y1,  double *y2, int *n, double *d)
{
/*

*/

	int  i;
	for (i = 0; i < *n; i++) {
		d[i]= dist(*x1, *x2, y1[i], y2[i]);
		}
}

double corr( double d, double *theta, int *corrmodel) {
        double val;
        

        switch(*corrmodel) {
		case 1: val  = exp(-d/theta[0]);
                break;

		case 2: if (d == 0) val = 1.0;
	                else val = (theta[0]/d)*sin(d/(theta[0]));
                break;

 		case 3: if (d <= theta[0]) val = 1.0-1.5*d/theta[0]+0.5*R_pow_di(d/theta[0], 3);
	                else val = 0.0;
                break;

 		case 4: val = 1.0/(1+R_pow(d/theta[0],2));
                break;

		default: if (d == 0) val = 1.0;
	                else val  = 0.0; /*  only for safety*/
	}
	if (val >= 1.0 ) {
	   	Rprintf(" rho func theta  %f \n", theta[0]);
		}


	
        return(val);
}    


double h12(double y1, double y2, double rho)
{
	double val;
	
	val = 1.0-rho*rho;
	val = (y2-rho*y1)/sqrt(val);
	return(val);
	

}


double e2f(double a)
/*
Marginal transformation
from exponential distribution to Frechet distribution
*/
{

	double val;
	
	val = -1.0 /a;
	val = -1.0/log(1-exp(val));   

	return(val);
}


double pfr(double y)
/*
evaluate the CDF of a unit Frechet r.v. at point y.

*/
{
	double val;
	

        val = exp(-1.0/y);
	
	return(val);
}



double overlarea(double d, double *radius)
  /*Area of Overlapping Circles of the same radius */
{
  double val, radius2;
  
  if (2 *(*radius) > d) {
    radius2 =(*radius)*(*radius);
    
    val =0.5*d/ (*radius);
    val = 2.0 * radius2*acos(val)-0.5*d*sqrt(4*radius2-d*d);
    val = val/(radius2*M_PI);
    
  }	
  else {
    val = 0.0;	
  }
  
  return (val);
  
  
}

void overlareawrap(double *d, double *radius, double *val, int *n)
  /*Area of Overlapping Circles of the same radius wrapper */
{
  int i;
  
  for (i = 0; i < *n; i++) {
    val[i] = overlarea(d[i], radius);
  }
}



/* TEG model*/
double V(double y1, double y2, double alpha, double rho)
/*
evaluate    V where exp(-V) is the bivariate CDF 
of TEG at point y1,y2.
with 
alpha 
rho correlation

*/

{
	double  sumy1y2, val;
	
	sumy1y2 = y1+y2;
	val = 1.0- 2.0*(1.0+rho)*y1*y2/(sumy1y2*sumy1y2);
	val = sqrt(val);
	val = (1.0/y1+1.0/y2)*(1.0-0.5*alpha*(1.0-val));

	return(val);
}



double V1(double y1, double y2, double alpha, double rho)

{
	double  sumy1y2, val;
	

	sumy1y2 = y1+y2;
	val = 1.0- 2.0*(1.0+rho)*y1*y2/R_pow_di(sumy1y2, 2);
	val = -(1.0/(y1*y1))*(1.0-0.5*alpha*(1-sqrt(val)))-0.5*alpha*(1.0+rho)*(1.0/y1+1.0/y2)*R_pow(val,-0.5)*(y2*y2-y2*y1)/(R_pow_di(sumy1y2, 3));	
	return(val);
}


double V2(double y1, double y2, double alpha, double rho)

{
	double  val;	

	
	
	val = V1(y2, y1, alpha,rho);
	
	return(val);
}

double V12(double y1, double y2, double alpha, double rho)

{
	double  sumy1y2, tmp,val;
	

	sumy1y2 = y1+y2;
	tmp = 1.0- 2.0*((1.0+rho)*y1*y2)/(R_pow_di(sumy1y2, 2));
	val = 0.5*(alpha*(1.0+rho))/(y1*y1)*(y1*y1-y1*y2)/R_pow_di(sumy1y2, 3)*R_pow(tmp,-0.5);
	val = val + 0.5*(alpha*(1.0+rho))/(y2*y2)*(y2*y2-y1*y2)/(R_pow_di(sumy1y2, 3))*R_pow(tmp,0.-0.5);
	val = val + 0.5*(alpha*(1.0+rho)*(1.0+rho))*(1.0/y1+1.0/y2)*(y1*y2*(y1-y2)*(y1-y2))/(R_pow_di(sumy1y2, 6))*R_pow(tmp,-1.5);
	val = val + 0.5*(alpha*(1.0+rho))*(1.0/y1+1.0/y2)*(y1*y1+y2*y2-4.0*y1*y2)/(R_pow_di(sumy1y2, 4))*R_pow(tmp,-0.5);

	return(val);
}




double G(double y1, double y2, double alpha, double rho)

{
	double   val;
	

	
	val = exp(-V(y1,y2,alpha,rho));
	return(val);
}

double G1(double y1, double y2, double alpha, double rho)

{
	double   val;	
	val = -G(y1,y2,alpha,rho)*V1(y1,y2,alpha,rho);
	return(val);
}


double G2(double y1, double y2, double alpha, double rho)

{
	double   val;		
	val = -G(y1,y2,alpha,rho)*V2(y1,y2,alpha,rho);
	return(val);
}


double g(double y1, double y2, double alpha, double rho)

{
	double   val;
	

	val = G(y1,y2,alpha,rho)*(V1(y1,y2,alpha,rho)*V2(y1,y2,alpha,rho)-V12(y1,y2,alpha,rho));
	
	
	return(val);
}


void gwrap(double *y1, double *y2, double *alpha,  double *rho, double *val)


{
	*val=g(*y1, *y2, *alpha,  *rho);

}

void Gwrap(double *y1, double *y2, double *alpha,  double *rho, double *val)


{
	*val=G(*y1, *y2, *alpha,  *rho);

}

double B(double y1, double y2, double alpha, double rho)
 
{
  double  l1, l2, val;
  
  
  l1 = e2f(y1);
  l2 = e2f(y2);
  
  val = exp(-1/y1) +exp(-1/y2) -1.0+exp(-V(l1,l2,alpha,rho));
  return(val);
}

double B1(double y1, double y2, double alpha, double rho)
 
  
{
  double  l1, l2, val;
  
  l1 = e2f(y1);
  l2 = e2f(y2);
  val = exp(-1/y1) / (y1*y1); 
  val = val*(1+G(l1, l2, alpha,rho)*V1(l1, l2, alpha, rho)*l1*l1/(1-exp(-1/y1)));
  return(val);
}

double B2(double y1, double y2, double alpha, double rho)
 
{
  double   val;
  
  val = B1(y2, y1, alpha, rho);
  return(val);
}

double b(double y1, double y2, double alpha, double rho)
  
{
  
  double  l1, l2, val;
  
  l1 = e2f(y1);
  l2 = e2f(y2);
  
  val = y1*y1*y2*y2*(1-exp(-1/y1))*(1-exp(-1/y2));
  val = (l1*l1*l2*l2*exp(-1/y1-1/y2))/val;
  
  val = val*(V1(l1,l2,alpha,rho)*V2(l1,l2,alpha,rho)-V12(l1,l2,alpha,rho))*G(l1,l2,alpha,rho);
  
  return(val);
}



void bwrap(double *y1, double *y2, double *alpha,  double *rho, double *val)
  
  
{
  *val=b(*y1, *y2, *alpha,  *rho);
  
}

void Bwrap(double *y1, double *y2, double *alpha,  double *rho, double *val)
  
  
{
  *val=B(*y1, *y2, *alpha,  *rho);
  
}





 
 /*Brown-resniek*/
 
 double over(double d,double *sigma, double *theta)
 {
   double val;
   
   val= (*sigma)*sqrt(2*(1-exp(-R_pow(d/(*theta),1))));
   
   return(val);
 }

 
double VS(double y1, double y2, double tau)
  /*
   evaluate    V where exp(-V) is the bivariate CDF
   of a Brown-resniek model  can also considred as Smith but we have to but the original funcion of  tau:
   with tau mah. distance
   
   */
{
  double val1, val2, val;
  
  val1= (0.5*tau)+(1/tau)*log(y2/y1);
  val2= (0.5*tau)+(1/tau)*log(y1/y2);
  val = (1/y1)*pnorm(val1, 0, 1, 1,0)+(1/y2)*pnorm(val2, 0, 1, 1,0);
  
  return(val);
}

double VS1(double y1, double y2, double tau)
  /*
   evaluate the partial derivative of   V, where exp(-V) is the bivariate CDF
   a Brown-resniek process at point y1,y2.
   */
  
{
  double val1, val2, val;
  
  val1= (0.5*tau)+(1/tau)*log(y2/y1);
  val2= (0.5*tau)+(1/tau)*log(y1/y2);
  
  val = -(pnorm(val1, 0, 1,1, 0)/(y1*y1)) - (dnorm(val1, 0, 1, 0)/(tau*y1*y1)) + (dnorm(val2, 0, 1, 0)/(tau*y1*y2));
  
  return(val);
}


double VS2(double y1, double y2, double tau)
  /*
   evaluate the partial derivative of   V, where exp(-V) is the bivariate CDF
   a Brown-resniek process at point y1,y2.
   */
{
  double val1, val2, val;
  
  val1= (0.5*tau)+(1/tau)*log(y2/y1);
  val2= (0.5*tau)+(1/tau)*log(y1/y2);
  
  val = -(pnorm(val2, 0, 1,1, 0)/(y2*y2)) - (dnorm(val2, 0, 1, 0)/(tau*y2*y2)) + (dnorm(val1, 0, 1, 0)/(tau*y1*y2));
  
  return(val);
}

double VS12(double y1, double y2, double tau)
  /*
   evaluate the partial derivative of   V, where exp(-V) is the bivariate CDF
   a Brown-resniek process at point y1,y2.
   */
{
  double val1, val2, val;
  
  val1= (0.5*tau)+(1/tau)*log(y2/y1);
  val2= (0.5*tau)+(1/tau)*log(y1/y2);
  
  val = -(dnorm(val1, 0, 1, 0)/(tau*y1*y1*y2)) - ((-val1*exp(-0.5*val1*val1)/sqrt(2*M_PI))/(tau*tau*y1*y1*y2)) -(dnorm(val2, 0, 1, 0)/(tau*y1*y2*y2)) - ((-val2*exp(-0.5*val2*val2)/sqrt(2*M_PI))/(tau*tau*y2*y2*y1));
  
  return(val);
}

double K(double y1, double y2, double tau)
  /*
   evaluate the bivariate CDF of
   a Brown-resniek process  at point y1,y2.
   */
  
{
  double   val;
  
  val = exp(-VS(y1,y2,tau));
  
  return(val);
}
double K1(double y1, double y2, double tau)
  /*
   evaluate the partial derivative of   the bivariate CDF of
   a Brown-resniek process at point y1,y2.
   
   */
  
  
{
  double   val;
  
  val = -K(y1,y2,tau)*VS1(y1,y2,tau);
  
  return(val);
}
double K2(double y1, double y2, double tau)
  /*
   evaluate the partial derivative of   the bivariate CDF of
   a Brown-resniek process  at point y1,y2.
   
   */
  
{
  double   val;
  
  val = -K(y1,y2,tau)*VS2(y1,y2,tau);
  
  return(val);
}


double kk(double y1, double y2, double tau)
  /*
   evaluate the bivariate density of
   a Brown-resniek process at point y1,y2.
   */
  
{
  double   val;
  
  val = K(y1,y2,tau)*(VS1(y1,y2,tau)*VS2(y1,y2,tau)-VS12(y1,y2,tau));
  
  return(val);
}

void kwrap(double *y1, double *y2, double *tau, double *val)
  
  
{
  *val=kk(*y1, *y2, *tau );
  
}

void Kwrap(double *y1, double *y2, double *tau, double *val)
  
  
{
  *val=K(*y1, *y2, *tau );
  
}

/*
 Brown-resniek model
 */




/* Inverted Brown-resniek model */  

double Q(double y1, double y2, double tau)
  /*
   evaluate the bivariate CDF of
   a inverse max-stable process Brown-resniek model at point y1,y2.
   
   */
  
{
  double  l1, l2, val;
  
  
  l1 = e2f(y1);
  l2 = e2f(y2);
  
  val = exp(-1/y1) +exp(-1/y2) -1.0+exp(-VS(l1,l2,tau));
  
  return(val);
}

double Q1(double y1, double y2, double tau)
  /*
   evaluate the partial derivative of   the bivariate CDF of
   a inverse max-stable  process Brown-resniek model at point y1,y2 with respect to y1
   */
  
  
{
  double  l1, l2, val;
  
  l1 = e2f(y1);
  l2 = e2f(y2);
  val = exp(-1/y1) / (y1*y1);
  val = val*(1+K(l1, l2, tau)*VS1(l1, l2, tau)*l1*l1/(1-exp(-1/y1)));
  
  return(val);
}

double Q2(double y1, double y2, double tau)
  /*
   evaluate the partial derivative of   the bivariate CDF of
   a inverse max-stable  process Brown-resniek model  at point y1,y2 with respect to y2
   
   */
  
{
  double   val;
  
  val = Q1(y2, y1, tau);
  
  return(val);
}

double q(double y1, double y2, double tau)
  /*
   evaluate the bivariate density of
   a inverse max-stable process Brown-resniek model at point y1,y2.
   
   */
  
{
  
  double  l1, l2, val;
  
  l1 = e2f(y1);
  l2 = e2f(y2);
  
  val = y1*y1*y2*y2*(1-exp(-1/y1))*(1-exp(-1/y2));
  val = (l1*l1*l2*l2*exp(-1/y1-1/y2))/val;
  
  val = val*(VS1(l1,l2,tau)*VS2(l1,l2,tau)-VS12(l1,l2,tau))*K(l1,l2,tau);
  
  return(val);
}

void qwrap(double *y1, double *y2, double *tau, double *val)
  
  
{
  *val=q(*y1, *y2, *tau );
  
}

void Qwrap(double *y1, double *y2, double *tau, double *val)
  
  
{
  *val=Q(*y1, *y2, *tau );
  
}

/*
 END Brown-resniek model
 */




/*
 BEGIN MIXTURE FOR TEG+inverted Brown-resniek
 */
double snew(double y1, double y2, double alpha, double rho, double tau, double mix)
  /*
   evaluate the bivariate density of at point y1,y2.
   */
  
{
  double  val,u1, u2, v1, v2;
  
  
  u1 = y1/mix;
  u2 = y2/mix;
  v1 = y1/(1-mix);
  v2 = y2/(1-mix);
  
  val = 	g(u1, u2,alpha, rho)*Q(v1, v2,tau)/(mix*mix);
  val = val +  G1(u1,u2, alpha, rho)*Q2(v1, v2, tau)/(mix*(1-mix));
  val = val+  G2(u1,u2, alpha, rho)*Q1(v1, v2,tau)/(mix*(1-mix));
  val = val+  G(u1,u2, alpha, rho)*q(v1, v2, tau)/((1-mix)*(1-mix));
  return(val);
}




double Snew(double y1, double y2, double alpha, double rho, double tau, double mix)
  /*
   evaluate the bivariate cumulative distribution function of at point y1,y2.
   */
  
{
  double  val,u1, u2, v1, v2;
  
  
  u1 = y1/mix;
  u2 = y2/mix;
  v1 = y1/(1-mix);
  v2 = y2/(1-mix);
  
  val = 	G(u1, u2,alpha, rho)*Q(v1, v2, tau);
  return(val);
}

void swrap(double *y1, double *y2, double *alpha, double *rho, double *tau, double *mix, double *val)
  
  
{
  *val=snew(*y1, *y2, *alpha, *rho, *tau, *mix);
  
}


void Swrap(double *y1, double *y2, double *alpha, double *rho, double *tau, double *mix, double *val)
  
{
  *val=Snew(*y1, *y2, *alpha, *rho, *tau, *mix);
}
/*
 END MIXTURE TEG+inverted brown-resneik
 */





/* censored likleihood for TEG */
void  davgholpwl_censored(double *y, double *x1, double *x2,  int *n, int *ntimes, double *theta1, double *radius, double *val,  double *delta, int *corrmodel1 ,  double *threshold)
{
  /*
   
   x1, x2 vector of coordinates
   */
  double  alpha, d, indval, rho1;
  int  i, j, k;
  
  
  
  *val = 0;
  
  
  for (i = 0; i < (*n-1); i++) {
    for (j = (i+1); j < *n; j++){   
      d = dist(x1[i], x2[i], x1[j], x2[j]);
      if ( d <= *delta ) {
        rho1 = corr(d,theta1, corrmodel1); 
        
        alpha = overlarea(d, radius);	
        for (k = 0; k < *ntimes; k++) {		
          if (ISNA(y[i+(*n)*k]) || ISNA(y[j+(*n)*k])) {
            *val += 0.0;
          }	
          else {
            
            if ((y[i+(*n)*k] > *threshold) || (y[j+(*n)*k] > *threshold)) {
              indval = g(y[i+(*n)*k], y[j+(*n)*k], alpha,  rho1);
            }
            else {
              
              indval = G(*threshold,*threshold, alpha,  rho1);
            }		
            
            
            *val += log(indval);
          }	
        }
      }
    }	
  }
  
  
  
}





/* censored likleihood for  inverted TEG */
void  invdavgholpwl_censored(double *y, double *x1, double *x2,  int *n, int *ntimes, double *theta1, double *radius, double *val,  double *delta, int *corrmodel1 ,  double *threshold)
{
  /*
   
   x1, x2 vector of coordinates
   */
  double  alpha, d, indval, rho1;
  int  i, j, k;
  
  
  
  *val = 0;
  
  
  for (i = 0; i < (*n-1); i++) {
    for (j = (i+1); j < *n; j++){   
      d = dist(x1[i], x2[i], x1[j], x2[j]);
      if ( d <= *delta ) {
        rho1 = corr(d,theta1, corrmodel1); 
        
        alpha = overlarea(d, radius);	
        for (k = 0; k < *ntimes; k++) {		
          if (ISNA(y[i+(*n)*k]) || ISNA(y[j+(*n)*k])) {
            *val += 0.0;
          }	
          else {
            
            if ((y[i+(*n)*k] > *threshold) || (y[j+(*n)*k] > *threshold)) {
              indval = b(y[i+(*n)*k], y[j+(*n)*k], alpha,  rho1);
            }
            else {
              
              indval = B(*threshold,*threshold, alpha,  rho1);
            }		
            
            
            *val += log(indval);
          }	
        }
      }
    }	
  }
  
  
  
}





/*censored likleihood for brown-resniek */
void  brownpwl_censored(double *y, double *x1, double *x2,  int *n, int *ntimes, double *sigma, double *theta, double *val,  double *delta,  double *threshold)
{
/*

x1, x2 vector of coordinates
*/
	double   d, indval,tau;
	int  i, j, k;

	
    
	*val = 0;
	
	
	for (i = 0; i < (*n-1); i++) {
		for (j = (i+1); j < *n; j++){   
			d = dist(x1[i], x2[i], x1[j], x2[j]);
			if ( d <= *delta ) {
			  tau = over(d, sigma,theta); 
				  for (k = 0; k < *ntimes; k++) {		
				    if (ISNA(y[i+(*n)*k]) || ISNA(y[j+(*n)*k])) {
					*val += 0.0;
				  }	
				  else {

				    	if ((y[i+(*n)*k] > *threshold) || (y[j+(*n)*k] > *threshold)) {
						indval = kk(y[i+(*n)*k], y[j+(*n)*k],tau);
				  	}
					else {

						indval = K(*threshold,*threshold,tau);
					}		


	               			*val += log(indval);
	               		  }	
				}
			}
		}	
	}
	  
	
}



/* censored likleihood for inverted brown-resniek */
void  invbrownpwl_censored(double *y, double *x1, double *x2,  int *n, int *ntimes, double *sigma,  double *theta, double *val,  double *delta,  double *threshold)
{
  /*
   
   x1, x2 vector of coordinates
   */
  double  d, indval,tau;
  int  i, j, k;
  
  *val = 0;
  
  
  for (i = 0; i < (*n-1); i++) {
    for (j = (i+1); j < *n; j++){   
      d = dist(x1[i], x2[i], x1[j], x2[j]);
      if ( d <= *delta ) {
        tau = over(d, sigma, theta);
       
        for (k = 0; k < *ntimes; k++) {		
          if (ISNA(y[i+(*n)*k]) || ISNA(y[j+(*n)*k])) {
            *val += 0.0;
          }	
          else {
            
            if ((y[i+(*n)*k] > *threshold) || (y[j+(*n)*k] > *threshold)) {
              indval = q(y[i+(*n)*k], y[j+(*n)*k],tau);
            }
            else {
              
              indval = Q(*threshold,*threshold, tau);
            }		
            
            
            *val += log(indval);
          }	
        }
      }
    }	
  }
}





/* censored LIKLIHOOOD TRUNCATED WITH inverted Brown-resniek*/
void  TEGbrownpwl_censored(double *y, double *x1, double *x2,  int *n, int *ntimes, double *theta1, double *radius, double *sigma, double *theta2, double *val,  double *delta, int *corrmodel, double *mix, double *threshold)
{
  /*
   
   x1, x2 vector of coordinates
   */
  double  alpha, d, indval, rho, tau;
  int  i, j, k;
  
  
  
  *val = 0;
  
  
  for (i = 0; i < (*n-1); i++) {
    for (j = (i+1); j < *n; j++){
      d = dist(x1[i], x2[i], x1[j], x2[j]);
      if ( d <= *delta ) {
        rho = corr(d,theta1, corrmodel);
        alpha = overlarea(d, radius);
        tau = over(d,sigma,theta2);
        
        for (k = 0; k < *ntimes; k++) {		
          if (ISNA(y[i+(*n)*k]) || ISNA(y[j+(*n)*k])) {
            *val += 0.0;
          }	
          else {
            
            if ((y[i+(*n)*k] > *threshold) || (y[j+(*n)*k] > *threshold)) {
              indval = snew(y[i+(*n)*k], y[j+(*n)*k], alpha,rho,tau, *mix);
              
            }
            else {
              
              indval = Snew(*threshold,*threshold, alpha, rho, tau, *mix);
              
            }		
            
            *val += log(indval);
            
          }	
        }
      }
    }	
  }
  
}
/*
 END LIKLIHOOD TRUNCATED WITH inverted brown-resniek
 */


/*
 SMITH TAU
 */
double over_SMITH(double d,double *sigma)
{
  double val;
  
  val=d/sqrt(*sigma);
  
  return(val);
}

/*
 Censored LIKLIHOOD TRUNCATED WITH inverted smith
 */

void  sgjpwl_censored(double *y, double *x1, double *x2,  int *n, int *ntimes, double *theta, double *radius,double *sigma, double *val,  double *delta, int *corrmodel, double *mix, double *threshold)
{
  /*
   
   x1, x2 vector of coordinates
   */
  double  alpha, d, indval, rho, tau;
  int  i, j, k;
  
  
  
  *val = 0;
  
  
  for (i = 0; i < (*n-1); i++) {
    for (j = (i+1); j < *n; j++){
      d = dist(x1[i], x2[i], x1[j], x2[j]);
      if ( d <= *delta ) {
        rho = corr(d,theta, corrmodel);
        alpha = overlarea(d, radius);
        tau = over_SMITH(d,sigma);
        
        for (k = 0; k < *ntimes; k++) {		
          if (ISNA(y[i+(*n)*k]) || ISNA(y[j+(*n)*k])) {
            *val += 0.0;
          }	
          else {
            
            if ((y[i+(*n)*k] > *threshold) || (y[j+(*n)*k] > *threshold)) {
              indval = snew(y[i+(*n)*k], y[j+(*n)*k], alpha,rho,tau, *mix);
              
            }
            else {
              
              indval = Snew(*threshold,*threshold, alpha, rho, tau, *mix);
              
            }		
            
            
            *val += log(indval);
            
          }	
        }
      }
    }	
  }
  

}
/*
 END LIKLIHOOD TRUNCATED WITH inverted SMITH
 */




/*censored likleihood for smith model */
void  smithpwl_censored(double *y, double *x1, double *x2,  int *n, int *ntimes, double *sigma, double *val,  double *delta,  double *threshold)
{
  /*
   
   x1, x2 vector of coordinates
   */
  double   d, indval,tau;
  int  i, j, k;
  
  
  
  *val = 0;
  
  
  for (i = 0; i < (*n-1); i++) {
    for (j = (i+1); j < *n; j++){   
      d = dist(x1[i], x2[i], x1[j], x2[j]);
      if ( d <= *delta ) {
        tau = over_SMITH(d, sigma); 
        for (k = 0; k < *ntimes; k++) {		
          if (ISNA(y[i+(*n)*k]) || ISNA(y[j+(*n)*k])) {
            *val += 0.0;
          }	
          else {
            
            if ((y[i+(*n)*k] > *threshold) || (y[j+(*n)*k] > *threshold)) {
              indval = kk(y[i+(*n)*k], y[j+(*n)*k],tau);
            }
            else {
              
              indval = K(*threshold,*threshold,tau);
            }		
            
            
            *val += log(indval);
          }	
        }
      }
    }	
  }
  
  
}



/* censored likleihood for inverted smith model */
void  invsmithpwl_censored(double *y, double *x1, double *x2,  int *n, int *ntimes, double *sigma, double *val,  double *delta,  double *threshold)
{
  /*
   
   x1, x2 vector of coordinates
   */
  double  d, indval,tau;
  int  i, j, k;
  
  *val = 0;
  
  
  for (i = 0; i < (*n-1); i++) {
    for (j = (i+1); j < *n; j++){   
      d = dist(x1[i], x2[i], x1[j], x2[j]);
      if ( d <= *delta ) {
        tau = over_SMITH(d, sigma);
        
        for (k = 0; k < *ntimes; k++) {		
          if (ISNA(y[i+(*n)*k]) || ISNA(y[j+(*n)*k])) {
            *val += 0.0;
          }	
          else {
            
            if ((y[i+(*n)*k] > *threshold) || (y[j+(*n)*k] > *threshold)) {
              indval = q(y[i+(*n)*k], y[j+(*n)*k],tau);
            }
            else {
              
              indval = Q(*threshold,*threshold, tau);
            }		
            
            
            *val += log(indval);
          }	
        }
      }
    }	
  }
}





