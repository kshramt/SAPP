#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(ptspecf)( double*, int*, double*, double*, double*, double*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

SEXP ptspec(SEXP data, SEXP n, SEXP intval, SEXP pprd, SEXP prdm, SEXP prd, SEXP nfre, SEXP nt, SEXP  is)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15;
    int *i1,*i2,*i3,*i4;

    SEXP ans = R_NilValue, prb = R_NilValue, r1 = R_NilValue, rwx = R_NilValue, rwy = R_NilValue, phs = R_NilValue, wt = R_NilValue, ht = R_NilValue, w = R_NilValue, h = R_NilValue, g = R_NilValue;

    double *zprb, *zr1, *zrwx, *zrwy, *zphs, *zwt, *zht, *zw, *zh, *zg = NULL;

    int nh1, nnt;
    int i;

    d1 = NUMERIC_POINTER(data);
    i1 = INTEGER_POINTER(n);
    d2 = NUMERIC_POINTER(intval);
    d3 = NUMERIC_POINTER(pprd);
    d4 = NUMERIC_POINTER(prdm);
    d5 = NUMERIC_POINTER(prd);
    i2 = INTEGER_POINTER(nfre);
    i3 = INTEGER_POINTER(nt);
    i4 = INTEGER_POINTER(is);

    nh1 = *i2 + 1;
    nnt = *i3;
    PROTECT(ans = allocVector(VECSXP, 10));
    SET_VECTOR_ELT(ans, 0, prb = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, r1 = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, rwx = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 3, rwy = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 4, phs = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 5, wt = allocVector(REALSXP, nnt));
    SET_VECTOR_ELT(ans, 6, ht = allocVector(REALSXP, nnt));
    SET_VECTOR_ELT(ans, 7, w = allocVector(REALSXP, nh1));
    SET_VECTOR_ELT(ans, 8, h = allocVector(REALSXP, nh1));
    SET_VECTOR_ELT(ans, 9, g = allocVector(REALSXP, nh1));

    d6 = NUMERIC_POINTER(prb);
    d7 = NUMERIC_POINTER(r1);
    d8 = NUMERIC_POINTER(rwx);
    d9 = NUMERIC_POINTER(rwy);
    d10 = NUMERIC_POINTER(phs);
    d11 = NUMERIC_POINTER(wt);
    d12 = NUMERIC_POINTER(ht);
    d13 = NUMERIC_POINTER(w);
    d14 = NUMERIC_POINTER(h);
    d15 = NUMERIC_POINTER(g);

    F77_CALL(ptspecf) (d1,i1,d2,d3,d4,d5,i2,i3,i4,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15);

    zprb = REAL(prb);
    zr1 = REAL(r1);
    zrwx = REAL(rwx);
    zrwy = REAL(rwy);
    zphs = REAL(phs);
    zwt = REAL(wt);
    zht = REAL(ht);
    zw = REAL(w);
    zh = REAL(h);
    zg = REAL(g);

    *zprb = *d6;
    *zr1 = *d7;
    *zrwx = *d8;
    *zrwy = *d9;
    *zphs = *d10;
    for(i=0; i<nnt; i++) zwt[i] = d11[i];
    for(i=0; i<nnt; i++) zht[i] = d12[i];
    for(i=0; i<nh1; i++) zw[i] = d13[i];
    for(i=0; i<nh1; i++) zh[i] = d14[i];
    for(i=0; i<nh1; i++) zg[i] = d15[i];

    UNPROTECT(1);

    return ans;
}
