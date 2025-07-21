/* biv-nt.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b8 = 1.;

/* Selected portion of code taken from: */
/*    http://www.math.wsu.edu/faculty/genz/software/mvtdstpack.f */
/* to compute bivariate normal and Student's t distribution functions. */

/* Author: */
/*          Alan Genz */
/*          Department of Mathematics */
/*          Washington State University */
/*          Pullman, WA 99164-3113 */
/*          Email : alangenz@wsu.edu */

/* except for some auxiliary functions whose authors are indicated */
/* in the respective code below. */



/* *********************************************************************** */

doublereal mvbvn_(doublereal *lower, doublereal *upper, integer *infin, 
	doublereal *correl)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Local variables */
    extern doublereal mvbvu_(doublereal *, doublereal *, doublereal *);


/*     A function for computing bivariate normal probabilities. */

/*  Parameters */

/*     LOWER  REAL, array of lower integration limits. */
/*     UPPER  REAL, array of upper integration limits. */
/*     INFIN  INTEGER, array of integration limits flags: */
/*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
/*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
/*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
/*     CORREL REAL, correlation coefficient. */

    /* Parameter adjustments */
    --infin;
    --upper;
    --lower;

    /* Function Body */
    if (infin[1] == 2 && infin[2] == 2) {
	ret_val = mvbvu_(&lower[1], &lower[2], correl) - mvbvu_(&upper[1], &
		lower[2], correl) - mvbvu_(&lower[1], &upper[2], correl) + 
		mvbvu_(&upper[1], &upper[2], correl);
    } else if (infin[1] == 2 && infin[2] == 1) {
	ret_val = mvbvu_(&lower[1], &lower[2], correl) - mvbvu_(&upper[1], &
		lower[2], correl);
    } else if (infin[1] == 1 && infin[2] == 2) {
	ret_val = mvbvu_(&lower[1], &lower[2], correl) - mvbvu_(&lower[1], &
		upper[2], correl);
    } else if (infin[1] == 2 && infin[2] == 0) {
	d__1 = -upper[1];
	d__2 = -upper[2];
	d__3 = -lower[1];
	d__4 = -upper[2];
	ret_val = mvbvu_(&d__1, &d__2, correl) - mvbvu_(&d__3, &d__4, correl);
    } else if (infin[1] == 0 && infin[2] == 2) {
	d__1 = -upper[1];
	d__2 = -upper[2];
	d__3 = -upper[1];
	d__4 = -lower[2];
	ret_val = mvbvu_(&d__1, &d__2, correl) - mvbvu_(&d__3, &d__4, correl);
    } else if (infin[1] == 1 && infin[2] == 0) {
	d__1 = -upper[2];
	d__2 = -(*correl);
	ret_val = mvbvu_(&lower[1], &d__1, &d__2);
    } else if (infin[1] == 0 && infin[2] == 1) {
	d__1 = -upper[1];
	d__2 = -(*correl);
	ret_val = mvbvu_(&d__1, &lower[2], &d__2);
    } else if (infin[1] == 1 && infin[2] == 1) {
	ret_val = mvbvu_(&lower[1], &lower[2], correl);
    } else if (infin[1] == 0 && infin[2] == 0) {
	d__1 = -upper[1];
	d__2 = -upper[2];
	ret_val = mvbvu_(&d__1, &d__2, correl);
    } else {
	ret_val = 1.;
    }
    return ret_val;
} /* mvbvn_ */

doublereal mvbvu_(doublereal *sh, doublereal *sk, doublereal *r__)
{
    /* Initialized data */

    static struct {
	doublereal e_1[3];
	doublereal fill_2[7];
	doublereal e_3[6];
	doublereal fill_4[4];
	doublereal e_5[10];
	} equiv_21 = { .1713244923791705, .3607615730481384, 
		.4679139345726904, {0}, .04717533638651177, .1069393259953183,
		 .1600783285433464, .2031674267230659, .2334925365383547, 
		.2491470458134029, {0}, .01761400713915212, 
		.04060142980038694, .06267204833410906, .08327674157670475, 
		.1019301198172404, .1181945319615184, .1316886384491766, 
		.1420961093183821, .1491729864726037, .1527533871307259 };

#define w ((doublereal *)&equiv_21)

    static struct {
	doublereal e_1[3];
	doublereal fill_2[7];
	doublereal e_3[6];
	doublereal fill_4[4];
	doublereal e_5[10];
	} equiv_22 = { -.9324695142031522, -.6612093864662647, 
		-.238619186083197, {0}, -.9815606342467191, -.904117256370475,
		 -.769902674194305, -.5873179542866171, -.3678314989981802, 
		-.1252334085114692, {0}, -.9931285991850949, 
		-.9639719272779138, -.9122344282513259, -.8391169718222188, 
		-.7463319064601508, -.636053680726515, -.5108670019508271, 
		-.3737060887154196, -.2277858511416451, -.07652652113349733 };

#define x ((doublereal *)&equiv_22)


    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double asin(doublereal), sin(doublereal), exp(doublereal), sqrt(
	    doublereal);

    /* Local variables */
    static doublereal a, b, c__, d__, h__;
    static integer i__;
    static doublereal k;
    static integer lg;
    static doublereal as;
    static integer ng;
    static doublereal bs, hk, hs, sn, rs, xs, bvn, asr;
    extern doublereal mvphi_(doublereal *);


/*     A function for computing bivariate normal probabilities; */
/*       developed using */
/*         Drezner, Z. and Wesolowsky, G. O. (1989), */
/*         On the Computation of the Bivariate Normal Integral, */
/*         J. Stat. Comput. Simul.. 35 pp. 101-107. */
/*       with extensive modications for double precisions by */
/*         Alan Genz and Yihong Ge */
/*         Department of Mathematics */
/*         Washington State University */
/*         Pullman, WA 99164-3113 */
/*         Email : alangenz@wsu.edu */

/* BVN - calculate the probability that X is larger than SH and Y is */
/*       larger than SK. */

/* Parameters */

/*   SH  REAL, integration limit */
/*   SK  REAL, integration limit */
/*   R   REAL, correlation coefficient */
/*   LG  INTEGER, number of Gauss Rule Points and Weights */

/*     Gauss Legendre Points and Weights, N =  6 */
/*     Gauss Legendre Points and Weights, N = 12 */
/*     Gauss Legendre Points and Weights, N = 20 */
    if (abs(*r__) < .3f) {
	ng = 1;
	lg = 3;
    } else if (abs(*r__) < .75f) {
	ng = 2;
	lg = 6;
    } else {
	ng = 3;
	lg = 10;
    }
    h__ = *sh;
    k = *sk;
    hk = h__ * k;
    bvn = 0.;
    if (abs(*r__) < .925f) {
	hs = (h__ * h__ + k * k) / 2;
	asr = asin(*r__);
	i__1 = lg;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sn = sin(asr * (x[i__ + ng * 10 - 11] + 1) / 2);
	    bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
		    ;
	    sn = sin(asr * (-x[i__ + ng * 10 - 11] + 1) / 2);
	    bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
		    ;
	}
	d__1 = -h__;
	d__2 = -k;
	bvn = bvn * asr / 12.566370614359172 + mvphi_(&d__1) * mvphi_(&d__2);
    } else {
	if (*r__ < 0.) {
	    k = -k;
	    hk = -hk;
	}
	if (abs(*r__) < 1.) {
	    as = (1 - *r__) * (*r__ + 1);
	    a = sqrt(as);
/* Computing 2nd power */
	    d__1 = h__ - k;
	    bs = d__1 * d__1;
	    c__ = (4 - hk) / 8;
	    d__ = (12 - hk) / 16;
	    bvn = a * exp(-(bs / as + hk) / 2) * (1 - c__ * (bs - as) * (1 - 
		    d__ * bs / 5) / 3 + c__ * d__ * as * as / 5);
	    if (hk > -160.) {
		b = sqrt(bs);
		d__1 = -b / a;
		bvn -= exp(-hk / 2) * sqrt(6.283185307179586) * mvphi_(&d__1) 
			* b * (1 - c__ * bs * (1 - d__ * bs / 5) / 3);
	    }
	    a /= 2;
	    i__1 = lg;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = a * (x[i__ + ng * 10 - 11] + 1);
		xs = d__1 * d__1;
		rs = sqrt(1 - xs);
		bvn += a * w[i__ + ng * 10 - 11] * (exp(-bs / (xs * 2) - hk / 
			(rs + 1)) / rs - exp(-(bs / xs + hk) / 2) * (c__ * xs 
			* (d__ * xs + 1) + 1));
/* Computing 2nd power */
		d__1 = -x[i__ + ng * 10 - 11] + 1;
		xs = as * (d__1 * d__1) / 4;
		rs = sqrt(1 - xs);
		bvn += a * w[i__ + ng * 10 - 11] * exp(-(bs / xs + hk) / 2) * 
			(exp(-hk * (1 - rs) / ((rs + 1) * 2)) / rs - (c__ * 
			xs * (d__ * xs + 1) + 1));
	    }
	    bvn = -bvn / 6.283185307179586;
	}
	if (*r__ > 0.) {
	    d__1 = -max(h__,k);
	    bvn += mvphi_(&d__1);
	}
	if (*r__ < 0.) {
/* Computing MAX */
	    d__3 = -h__;
	    d__4 = -k;
	    d__1 = 0., d__2 = mvphi_(&d__3) - mvphi_(&d__4);
	    bvn = -bvn + max(d__1,d__2);
	}
    }
    ret_val = bvn;
    return ret_val;
} /* mvbvu_ */

#undef x
#undef w

doublereal mvphi_(doublereal *z__)
{
    /* Initialized data */

    static doublereal a[44] = { .610143081923200417926465815756,
	    -.434841272712577471828182820888,.176351193643605501125840298123,
	    -.060710795609249414860051215825,.017712068995694114486147141191,
	    -.004321119385567293818599864968,8.54216676887098678819832055e-4,
	    -1.2715509060916274262889394e-4,1.1248167243671189468847072e-5,
	    3.13063885421820972630152e-7,-2.70988068537762022009086e-7,
	    3.0737622701407688440959e-8,2.515620384817622937314e-9,
	    -1.02892992132031912759e-9,2.9944052119949939363e-11,
	    2.605178968726693629e-11,-2.634839924171969386e-12,
	    -6.43404509890636443e-13,1.12457401801663447e-13,
	    1.7281533389986098e-14,-4.264101694942375e-15,
	    -5.45371977880191e-16,1.58697607761671e-16,2.0899837844334e-17,
	    -5.900526869409e-18,-9.41893387554e-19,2.1497735647e-19,
	    4.6660985008e-20,-7.243011862e-21,-2.387966824e-21,1.91177535e-22,
	    1.20482568e-22,-6.72377e-25,-5.747997e-24,-4.28493e-25,
	    2.44856e-25,4.3793e-26,-8.151e-27,-3.089e-27,9.3e-29,1.74e-28,
	    1.6e-29,-8e-30,-2e-30 };

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal b;
    static integer i__;
    static doublereal p, t, bm, bp, xa;


/*     Normal distribution probabilities accurate to 1d-15. */
/*     Reference: J.L. Schonfelder, Math Comp 32(1978), pp 1232-1240. */


    xa = abs(*z__) / 1.414213562373095048801688724209;
    if (xa > 100.) {
	p = 0.;
    } else {
	t = (xa * 8 - 30) / (xa * 4 + 15);
	bm = 0.;
	b = 0.;
	for (i__ = 24; i__ >= 0; --i__) {
	    bp = b;
	    b = bm;
	    bm = t * b - bp + a[i__];
	}
	p = exp(-xa * xa) * (bm - bp) / 4;
    }
    if (*z__ > 0.) {
	p = 1 - p;
    }
    ret_val = p;
    return ret_val;
} /* mvphi_ */


