#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(momorif)(double*, int*, double*, int*, double*, double*, double*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, double*, double*, int*, int*);

SEXP momori(SEXP y, SEXP n, SEXP pai, SEXP np, SEXP zts, SEXP zte, SEXP tstart, SEXP ncount, SEXP nfunct, SEXP nlmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17,*d18, *d19, *d20;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7;

    SEXP ans = R_NilValue, ff = R_NilValue, x = R_NilValue, g = R_NilValue, pa = R_NilValue, ahaic = R_NilValue;
    SEXP t0 = R_NilValue, ti = R_NilValue, ak = R_NilValue, c = R_NilValue, p = R_NilValue, cls = R_NilValue, id = R_NilValue, rmd = R_NilValue, x1 = R_NilValue, h = R_NilValue, hf = R_NilValue, nl = R_NilValue;
    double *xff, *xx, *xg, *xpa, *xahaic, *xt0, *xti, *xak, *xc, *xp, *xcls, *xrmd, *xx1, *xh, *xhf = NULL;
    int      *xid, *xnl;

    int nnp, kn, nc, nlm;
    int i;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    d2 = NUMERIC_POINTER(pai);
    i2 = INTEGER_POINTER(np);
    d3 = NUMERIC_POINTER(zts);
    d4 = NUMERIC_POINTER(zte);
    d5 = NUMERIC_POINTER(tstart);
    i3 = INTEGER_POINTER(ncount);
    i4 = INTEGER_POINTER(nfunct);
    i7 = INTEGER_POINTER(nlmax);

    nnp = *i2;
    kn = (nnp-1)/3;
    nc = *i3;
    nlm = *i7;

    PROTECT(ans = allocVector(VECSXP, 17));
    SET_VECTOR_ELT(ans, 0, ff = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, x = allocVector(REALSXP, nnp*2));
    SET_VECTOR_ELT(ans, 2, g = allocVector(REALSXP, nnp*2));
    SET_VECTOR_ELT(ans, 3, pa = allocVector(REALSXP, nnp));
    SET_VECTOR_ELT(ans, 4, ahaic = allocVector(REALSXP, nc));
    SET_VECTOR_ELT(ans, 5, t0 = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 6, ti = allocVector(REALSXP, kn));
    SET_VECTOR_ELT(ans, 7, ak = allocVector(REALSXP, kn));
    SET_VECTOR_ELT(ans, 8, c = allocVector(REALSXP, kn));
    SET_VECTOR_ELT(ans, 9, p = allocVector(REALSXP, kn));
    SET_VECTOR_ELT(ans, 10, cls = allocVector(REALSXP, kn));
    SET_VECTOR_ELT(ans, 11, id = allocVector(INTSXP, nlm));
    SET_VECTOR_ELT(ans, 12, rmd = allocVector(REALSXP, nlm));
    SET_VECTOR_ELT(ans, 13, x1 = allocVector(REALSXP, nnp*nlm));
    SET_VECTOR_ELT(ans, 14, h = allocVector(REALSXP, nnp*nnp*2));
    SET_VECTOR_ELT(ans, 15, hf = allocVector(REALSXP, nnp*nnp*4));
    SET_VECTOR_ELT(ans, 16, nl = allocVector(INTSXP, 1));


    d6 = NUMERIC_POINTER(ff);
    d7 = NUMERIC_POINTER(x);
    d8 = NUMERIC_POINTER(g);
    d9 = NUMERIC_POINTER(pa);
    d10 = NUMERIC_POINTER(ahaic);
    d11 = NUMERIC_POINTER(t0);
    d12 = NUMERIC_POINTER(ti);
    d13 = NUMERIC_POINTER(ak);
    d14 = NUMERIC_POINTER(c);
    d15 = NUMERIC_POINTER(p);
    d16 = NUMERIC_POINTER(cls);
    i5 = INTEGER_POINTER(id);
    d17 = NUMERIC_POINTER(rmd);
    d18 = NUMERIC_POINTER(x1);
    d19 = NUMERIC_POINTER(h);
    d20 = NUMERIC_POINTER(hf);
    i6 = INTEGER_POINTER(nl);

    F77_CALL(momorif) (d1,i1,d2,i2,d3,d4,d5,i3,i4,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,i5,d17,d18,d19,d20,i6,i7);

    xff = REAL(ff);
    xx = REAL(x);
    xg = REAL(g);
    xpa = REAL(pa);
    xahaic = REAL(ahaic);
    xt0 = REAL(t0);
    xti = REAL(ti);
    xak = REAL(ak);
    xc = REAL(c);
    xp = REAL(p);
    xcls = REAL(cls);
    xid = INTEGER(id);
    xrmd = REAL(rmd);
    xx1 = REAL(x1);
    xh = REAL(h);
    xhf = REAL(hf);
    xnl = INTEGER(nl);

    *xff = *d6;
    for(i=0; i<nnp*2; i++) xx[i] = d7[i];
    for(i=0; i<nnp*2; i++) xg[i] = d8[i];
    for(i=0; i<nnp; i++) xpa[i] = d9[i];
    for(i=0; i<nc; i++) xahaic[i] = d10[i];
    *xt0 = *d11;
    for(i=0; i<kn; i++) xti[i] = d12[i];
    for(i=0; i<kn; i++) xak[i] = d13[i];
    for(i=0; i<kn; i++) xc[i] = d14[i];
    for(i=0; i<kn; i++) xp[i] = d15[i];
    for(i=0; i<kn; i++) xcls[i] = d16[i];
    for(i=0; i<nlm; i++) xid[i] = i5[i];
    for(i=0; i<nlm; i++) xrmd[i] = d17[i];
    for(i=0; i<nnp*nlm; i++) xx1[i] = d18[i];
    for(i=0; i<nnp*nnp*2; i++) xh[i] = d19[i];
    for(i=0; i<nnp*nnp*4; i++) xhf[i] = d20[i];
    *xnl = *i6;

    UNPROTECT(1);

    return ans;
}
