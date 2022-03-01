// #include <R/R_ext/Applic.h>
#include <R.h>
#include <R_ext/Applic.h>

/*
 lmm is an integer variable.
	 On entry m is the maximum number of variable metric
	    corrections allowed in the limited memory matrix.
	 On exit m is unchanged.

nbd is an integer array of dimension n.
	 On entry nbd represents the type of bounds imposed on the
	   variables, and must be specified as follows:
	   nbd(i)=0 if x(i) is unbounded,
		  1 if x(i) has only a lower bound,
		  2 if x(i) has both lower and upper bounds,
		  3 if x(i) has only an upper bound.
	 On exit nbd is unchanged.
*/

void lbfgsb_f_ (int *n, int *lmm, double *x, double *lower,
		double *upper, int *nbd, double *Fmin, optimfn fn,
		optimgr gr, int *fail, void *ex, double *factr,
		double *pgtol, int *fncount, int *grcount,
		int *maxit)
{
  /* call
     void lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
     double *Fmin, optimfn fn, optimgr gr, int *fail, void *ex,
     double factr, double pgtol, int *fncount, int *grcount,
     int maxit, char *msg, int trace, int nREPORT);
  */
  char msg[60];
  int trace=0;
  int nREPORT=1;

  lbfgsb(*n, *lmm, x, lower, upper, nbd,
	 Fmin, fn, gr, fail, ex,
	 *factr, *pgtol, fncount, grcount,
	 *maxit,
	 msg, trace, nREPORT);

}

