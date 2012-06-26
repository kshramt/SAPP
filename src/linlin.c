#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(linlinf)(int*, double*, int*, double*, int*, int*, double*, double*, int*, int*, int*, int*, int*,  int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, int*, int*);

SEXP linlin(SEXP n, SEXP x, SEXP opt, SEXP t, SEXP nn, SEXP mm, SEXP xx, SEXP yy, SEXP kkx, SEXP kky, SEXP kmax, SEXP kkc, SEXP kkt, SEXP nlmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13;

    SEXP ans = R_NilValue, x1 = R_NilValue, x2 = R_NilValue, aic = R_NilValue, f = R_NilValue, prb = R_NilValue,  r1= R_NilValue, rwx = R_NilValue, rwy = R_NilValue, px = R_NilValue, pg = R_NilValue, phs = R_NilValue, id = R_NilValue, rmd = R_NilValue, eee = R_NilValue, nl = R_NilValue, ier = R_NilValue;

    double *xx1, *xx2, *xaic, *xf, *xprb, *xr1, *xrwx, *xrwy, *xphs, *xpx, *xpg, *xrmd, *xeee = NULL;
    int      *xid, *xnl, *xier = NULL;

    int np, nlm;
    int i;

    i1 = INTEGER_POINTER(n);
    d1 = NUMERIC_POINTER(x);
    i2 = INTEGER_POINTER(opt);
    d2 = NUMERIC_POINTER(t);
    i3 = INTEGER_POINTER(nn);
    i4 = INTEGER_POINTER(mm);
    d3 = NUMERIC_POINTER(xx);
    d4 = NUMERIC_POINTER(yy);
    i5 = INTEGER_POINTER(kkx);
    i6 = INTEGER_POINTER(kky);
    i7 = INTEGER_POINTER(kmax);
    i8 = INTEGER_POINTER(kkc);
    i9 = INTEGER_POINTER(kkt);
    i10 = INTEGER_POINTER(nlmax);

    np = *i1;
    nlm = *i10;

    PROTECT(ans = allocVector(VECSXP, 16));
    SET_VECTOR_ELT(ans, 0, x1 = allocVector(REALSXP, np));
    SET_VECTOR_ELT(ans, 1, x2 = allocVector(REALSXP, np));
    SET_VECTOR_ELT(ans, 2, aic = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 3, f = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 4, prb = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 5, r1 = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 6, rwx = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 7, rwy = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 8, phs = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 9, px = allocVector(REALSXP, 5*np));
    SET_VECTOR_ELT(ans, 10, pg = allocVector(REALSXP, 5*np));
    SET_VECTOR_ELT(ans, 11, id = allocVector(INTSXP, nlm));
    SET_VECTOR_ELT(ans, 12, rmd = allocVector(REALSXP, nlm));
    SET_VECTOR_ELT(ans, 13, eee = allocVector(REALSXP, nlm));
    SET_VECTOR_ELT(ans, 14, nl = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 15, ier = allocVector(INTSXP, 1));

    d5 = NUMERIC_POINTER(x1);
    d6 = NUMERIC_POINTER(x2);
    d7 = NUMERIC_POINTER(aic);
    d8 = NUMERIC_POINTER(f);
    d9 = NUMERIC_POINTER(prb);
    d10 = NUMERIC_POINTER(r1);
    d11 = NUMERIC_POINTER(rwx);
    d12 = NUMERIC_POINTER(rwy);
    d13 = NUMERIC_POINTER(phs);
    d14 = NUMERIC_POINTER(px);
    d15 = NUMERIC_POINTER(pg);
    i11 = INTEGER_POINTER(id);
    d16 = NUMERIC_POINTER(rmd);
    d17 = NUMERIC_POINTER(eee);
    i12 = INTEGER_POINTER(nl);
    i13 = INTEGER_POINTER(ier);

    F77_CALL(linlinf) (i1,d1,i2,d2,i3,i4,d3,d4,i5,i6,i7,i8,i9,i10,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,i11,d16,d17,i12,i13);

    xx1 = REAL(x1);
    xx2 = REAL(x2);
    xaic = REAL(aic);
    xf = REAL(f);
    xprb = REAL(prb);
    xr1 = REAL(r1);
    xrwx = REAL(rwx);
    xrwy = REAL(rwy);
    xphs = REAL(phs);
    xpx = REAL(px);
    xpg = REAL(pg);
    xid = INTEGER(id);
    xrmd = REAL(rmd);
    xeee = REAL(eee);
    xnl = INTEGER(nl);
    xier = INTEGER(ier);

    for(i=0; i<np; i++) xx1[i] = d5[i];
    for(i=0; i<np; i++) xx2[i] = d6[i];
    *xaic = *d7;
    *xf = *d8;
    *xprb = *d9;
    *xr1 = *d10;
    *xrwx = *d11;
    *xrwy = *d12;
    *xphs = *d13;
    for(i=0; i<5*np; i++) xpx[i] = d14[i];
    for(i=0; i<5*np; i++) xpg[i] = d15[i];
    for(i=0; i<nlm; i++) xid[i] = i11[i];
    for(i=0; i<nlm; i++) xrmd[i] = d16[i];
    for(i=0; i<nlm; i++) xeee[i] = d17[i];
    *xnl = *i12;
    *xier = *i13;

    UNPROTECT(1);

    return ans;
}
