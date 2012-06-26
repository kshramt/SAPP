#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(respoif)(double*, double*, double*, double*, double*, int*, double*, int*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, int*);

SEXP respoi(SEXP time, SEXP mag, SEXP dep, SEXP xp, SEXP yp, SEXP nd, SEXP param, SEXP np, SEXP zts, SEXP zte, SEXP tstart, SEXP thresh)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16;
    int *i1,*i2,*i3,*i4;

    SEXP ans = R_NilValue, mag1 = R_NilValue, dep1 = R_NilValue, xp1 = R_NilValue, yp1 = R_NilValue, ntstar = R_NilValue, xx = R_NilValue, x = R_NilValue, nn = R_NilValue;
    double *zmag1, *zdep1, *zxp1, *zyp1, *zxx, *zx = NULL;
    int      *zntstar, *znn = NULL;

    int i, nnd;

    d1 = NUMERIC_POINTER(time);
    d2 = NUMERIC_POINTER(mag);
    d3 = NUMERIC_POINTER(dep);
    d4 = NUMERIC_POINTER(xp);
    d5 = NUMERIC_POINTER(yp);
    i1 = INTEGER_POINTER(nd);
    d6 = NUMERIC_POINTER(param);
    i2 = INTEGER_POINTER(np);
    d7 = NUMERIC_POINTER(zts);
    d8= NUMERIC_POINTER(zte);
    d9= NUMERIC_POINTER(tstart);
    d10= NUMERIC_POINTER(thresh);

    nnd = *i1;
    PROTECT(ans = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(ans, 0, mag1 = allocVector(REALSXP, nnd));
    SET_VECTOR_ELT(ans, 1, dep1 = allocVector(REALSXP, nnd));
    SET_VECTOR_ELT(ans, 2, xp1= allocVector(REALSXP, nnd));
    SET_VECTOR_ELT(ans, 3, yp1= allocVector(REALSXP, nnd));
    SET_VECTOR_ELT(ans, 4, ntstar = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 5, xx = allocVector(REALSXP, nnd));
    SET_VECTOR_ELT(ans, 6, x= allocVector(REALSXP, nnd));
    SET_VECTOR_ELT(ans, 7, nn = allocVector(INTSXP, 1));

    d11= NUMERIC_POINTER(mag1);
    d12= NUMERIC_POINTER(dep1);
    d13= NUMERIC_POINTER(xp1);
    d14= NUMERIC_POINTER(yp1);
    i3 = INTEGER_POINTER(ntstar);
    d15= NUMERIC_POINTER(xx);
    d16= NUMERIC_POINTER(x);
    i4 = INTEGER_POINTER(nn);

    F77_CALL(respoif) (d1,d2,d3,d4,d5,i1,d6,i2,d7,d8,d9,d10,d11,d12,d13,d14,i3,d15,d16,i4);

    zmag1 = REAL(mag1);
    zdep1 = REAL(dep1);
    zxp1 = REAL(xp1);
    zyp1 = REAL(yp1);
    zntstar = INTEGER(ntstar);
    zxx = REAL(xx);
    zx = REAL(x);
    znn = INTEGER(nn);

    for(i=0; i<nnd; i++) zmag1[i] = d11[i];
    for(i=0; i<nnd; i++) zdep1[i] = d12[i];
    for(i=0; i<nnd; i++) zxp1[i] = d13[i];
    for(i=0; i<nnd; i++) zyp1[i] = d14[i];
    *zntstar = *i3;
    for(i=0; i<nnd; i++) zxx[i] = d15[i];
    for(i=0; i<nnd; i++) zx[i] = d16[i];
    *znn = *i4;

    UNPROTECT(1);

    return ans;
}
