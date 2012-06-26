#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(eptrenf)(double*, double*, int*, int*, int*, int*, double*, double*, double*, double*, int*,  double*, double*, double*, double*, int*, double*, double*, int*, int*, int*, int*);

SEXP eptren(SEXP y, SEXP t, SEXP n, SEXP nfunct, SEXP npara, SEXP nsub, SEXP cycle, SEXP nmax, SEXP np, SEXP nlmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10;

    SEXP ans = R_NilValue, xa = R_NilValue, aic = R_NilValue, aicm = R_NilValue, morder = R_NilValue, xval = R_NilValue, fval = R_NilValue, x = R_NilValue, g = R_NilValue, id = R_NilValue, rmd = R_NilValue, eee = R_NilValue, nl = R_NilValue;
    double *xxa, *xaic, *xaicm, *xxval, *xfval, *xx, *xg, *xrmd, *xeee = NULL;
    int      *xmorder, *xid, *xnl = NULL;

    int npa, nm, nnp, nlm;
    int i;

    d1 = NUMERIC_POINTER(y);
    d2 = NUMERIC_POINTER(t);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(nfunct);
    i3 = INTEGER_POINTER(npara);
    i4 = INTEGER_POINTER(nsub);
    d3 = NUMERIC_POINTER(cycle);
    i8 = INTEGER_POINTER(nmax);
    i9 = INTEGER_POINTER(np);
    i10 = INTEGER_POINTER(nlmax);

    npa = *i3;
    nm = *i8;
    nnp = *i9;
    nlm = *i10;

    PROTECT(ans = allocVector(VECSXP, 12));
    SET_VECTOR_ELT(ans, 0, xa = allocVector(REALSXP, nm*npa));
    SET_VECTOR_ELT(ans, 1, aic = allocVector(REALSXP, npa));
    SET_VECTOR_ELT(ans, 2, aicm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 3, morder = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 4, xval = allocVector(REALSXP, nnp));
    SET_VECTOR_ELT(ans, 5, fval = allocVector(REALSXP, nnp));
    SET_VECTOR_ELT(ans, 6, x = allocVector(REALSXP, nm*npa));
    SET_VECTOR_ELT(ans, 7, g = allocVector(REALSXP, nm*npa));
    SET_VECTOR_ELT(ans, 8, id = allocVector(INTSXP, nlm));
    SET_VECTOR_ELT(ans, 9, rmd = allocVector(REALSXP, nlm));
    SET_VECTOR_ELT(ans, 10, eee = allocVector(REALSXP, nlm));
    SET_VECTOR_ELT(ans, 11, nl = allocVector(INTSXP, 1));

    d4 = NUMERIC_POINTER(xa);
    d5 = NUMERIC_POINTER(aic);
    d6 = NUMERIC_POINTER(aicm);
    i5 = INTEGER_POINTER(morder);
    d7 = NUMERIC_POINTER(xval);
    d8 = NUMERIC_POINTER(fval);
    d9 = NUMERIC_POINTER(x);
    d10 = NUMERIC_POINTER(g);
    i6 = INTEGER_POINTER(id);
    d11 = NUMERIC_POINTER(rmd);
    d12 = NUMERIC_POINTER(eee);
    i7 = INTEGER_POINTER(nl);

    F77_CALL(eptrenf) (d1,d2,i1,i2,i3,i4,d3,d4,d5,d6,i5,d7,d8,d9,d10,i6,d11,d12,i7,i8,i9,i10);

    xxa = REAL(xa);
    xaic = REAL(aic);
    xaicm = REAL(aicm);
    xmorder = INTEGER(morder);
    xxval = REAL(xval);
    xfval = REAL(fval);
    xx = REAL(x);
    xg = REAL(g);
    xid = INTEGER(id);
    xrmd = REAL(rmd);
    xeee = REAL(eee);
    xnl = INTEGER(nl);

    for(i=0; i<nm*npa; i++) xxa[i] = d4[i];
    for(i=0; i<npa; i++) xaic[i] = d5[i];
    *xaicm = *d6;
    *xmorder = *i5;
    for(i=0; i<nnp; i++) xxval[i] = d7[i];
    for(i=0; i<nnp; i++) xfval[i] = d8[i];
    for(i=0; i<nm*npa; i++) xx[i] = d9[i];
    for(i=0; i<nm*npa; i++) xg[i] = d10[i];
    for(i=0; i<nlm; i++) xid[i] = i6[i];
    for(i=0; i<nlm; i++) xrmd[i] = d11[i];
    for(i=0; i<nlm; i++) xeee[i] = d12[i];
    *xnl = *i7;

    UNPROTECT(1);

    return ans;
}
