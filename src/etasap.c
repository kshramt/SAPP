#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(etasapf)(double*, double*, int*, double*, double*, double*, int*, double*, double*, double*, int*, int*, double*, double*, double*, double*, int*, double*, double*, int*, int*);

SEXP etasap(SEXP xx, SEXP mag, SEXP nn, SEXP refer, SEXP thresh, SEXP param, SEXP np, SEXP zts, SEXP zte, SEXP tstart, SEXP nfunct, SEXP app, SEXP nlmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7;

    SEXP ans = R_NilValue, f = R_NilValue, x = R_NilValue, g = R_NilValue, aic2 = R_NilValue, id = R_NilValue, ee = R_NilValue, x1 = R_NilValue, nl = R_NilValue;
    double *zf, *zx, *zg, *zaic2, *zee, *zx1 = NULL;
    int      *zid, *znl;

    int i, nnp, nlm;

    d1 = NUMERIC_POINTER(xx);
    d2 = NUMERIC_POINTER(mag);
    i1 = INTEGER_POINTER(nn);
    d3 = NUMERIC_POINTER(refer);
    d4 = NUMERIC_POINTER(thresh);
    d5 = NUMERIC_POINTER(param);
    i2 = INTEGER_POINTER(np);
    d6 = NUMERIC_POINTER(zts);
    d7 = NUMERIC_POINTER(zte);
    d8= NUMERIC_POINTER(tstart);
    i3 = INTEGER_POINTER(nfunct);
    i4 = INTEGER_POINTER(app);
    i7 = INTEGER_POINTER(nlmax); 

    nnp = *i2;
    nlm = *i7;
    PROTECT(ans = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(ans, 0, f = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, x = allocVector(REALSXP, nnp));
    SET_VECTOR_ELT(ans, 2, g = allocVector(REALSXP, nnp));
    SET_VECTOR_ELT(ans, 3, aic2= allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 4, id = allocVector(INTSXP, nlm));
    SET_VECTOR_ELT(ans, 5, ee = allocVector(REALSXP, nlm));
    SET_VECTOR_ELT(ans, 6, x1 = allocVector(REALSXP, nlm*nnp));
    SET_VECTOR_ELT(ans, 7, nl = allocVector(INTSXP, 1));

    d9= NUMERIC_POINTER(f);
    d10= NUMERIC_POINTER(x);
    d11= NUMERIC_POINTER(g);
    d12= NUMERIC_POINTER(aic2);
    i5 = INTEGER_POINTER(id); 
    d13= NUMERIC_POINTER(ee);
    d14= NUMERIC_POINTER(x1);
    i6 = INTEGER_POINTER(nl); 

    F77_CALL(etasapf) (d1,d2,i1,d3,d4,d5,i2,d6,d7,d8,i3,i4,d9,d10,d11,d12,i5,d13,d14,i6,i7);

    zf = REAL(f);
    zx = REAL(x);
    zg = REAL(g);
    zaic2 = REAL(aic2);
    zid = INTEGER(id);
    zee = REAL(ee);
    zx1 = REAL(x1);
    znl = INTEGER(nl);

    *zf = *d9;
    for(i=0; i<nnp; i++) zx[i] = d10[i];
    for(i=0; i<nnp; i++) zg[i] = d11[i];
    *zaic2 = *d12;
    for(i=0; i<nlm; i++) zid[i] = i5[i];
    for(i=0; i<nlm; i++) zee[i] = d13[i];
    for(i=0; i<nlm*nnp; i++) zx1[i] = d14[i];
    *znl = *i6;

    UNPROTECT(1);

    return ans;
}
