#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(etasimf)(int*, double*, double*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

SEXP etasim(SEXP ic, SEXP bvalue, SEXP tstart, SEXP nd, SEXP thresh, SEXP refer, SEXP a, SEXP b, SEXP c, SEXP d, SEXP p, SEXP mag1, SEXP time1)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13;
    int *i1,*i2;

    SEXP ans = R_NilValue, mag = R_NilValue, x = R_NilValue, probx = R_NilValue;
    double *zmag, *zx, *zprobx = NULL;

    int i, nnd;

    i1 = INTEGER_POINTER(ic);
    d1 = NUMERIC_POINTER(bvalue);
    d2 = NUMERIC_POINTER(tstart);
    i2 = INTEGER_POINTER(nd);
    d3 = NUMERIC_POINTER(thresh);
    d4 = NUMERIC_POINTER(refer);
    d5 = NUMERIC_POINTER(a);
    d6 = NUMERIC_POINTER(b);
    d7 = NUMERIC_POINTER(c);
    d8 = NUMERIC_POINTER(d);
    d9 = NUMERIC_POINTER(p);
    d10 = NUMERIC_POINTER(mag1);
    d11 = NUMERIC_POINTER(time1);

    nnd = *i2;
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, mag = allocVector(REALSXP, nnd));
    SET_VECTOR_ELT(ans, 1, x = allocVector(REALSXP, nnd));
    SET_VECTOR_ELT(ans, 2, probx= allocVector(REALSXP, 1));

    d12 = NUMERIC_POINTER(x);
    d13 = NUMERIC_POINTER(probx);

    F77_CALL(etasimf) (i1,d1,d2,i2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13);

    zmag = REAL(mag);
    zx = REAL(x);
    zprobx = REAL(probx);

    for(i=0; i<nnd; i++) zmag[i] = d10[i];
    for(i=0; i<nnd; i++) zx[i] = d12[i];
    *zprobx = *d13;

    UNPROTECT(1);

    return ans;
}
