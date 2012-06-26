#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(etarppf)(double*, double*, double*, int*, double*, int*, double*, double*, double*,  double*, int*);

SEXP etarpp(SEXP time, SEXP mag, SEXP refer, SEXP nn, SEXP param, SEXP np, SEXP zts, SEXP zte, SEXP tstart)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8;
    int *i1,*i2,*i3;

    SEXP ans = R_NilValue, x = R_NilValue, ntst = R_NilValue;
    double *zx = NULL;
    int      *zntst = NULL;

    int i, m;

    d1 = NUMERIC_POINTER(time);
    d2 = NUMERIC_POINTER(mag);
    d3 = NUMERIC_POINTER(refer);
    i1 = INTEGER_POINTER(nn);
    d4 = NUMERIC_POINTER(param);
    i2 = INTEGER_POINTER(np);
    d5 = NUMERIC_POINTER(zts);
    d6 = NUMERIC_POINTER(zte);
    d7 = NUMERIC_POINTER(tstart);

    m = *i1;
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, x= allocVector(REALSXP, m));
    SET_VECTOR_ELT(ans, 1, ntst = allocVector(INTSXP, 1));

    d8= NUMERIC_POINTER(x);
    i3 = INTEGER_POINTER(ntst);

    F77_CALL(etarppf) (d1,d2,d3,i1,d4,i2,d5,d6,d7,d8,i3);

    zx = REAL(x);
    zntst = INTEGER(ntst);

    for(i=0; i<m; i++) zx[i] = d8[i];
    *zntst = *i3;

    UNPROTECT(1);

    return ans;
}
