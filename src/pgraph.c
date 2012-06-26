#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(pgraphf)(int*, int*, double*, int*, int*, double*, double*, double*, double*, int*, double*, double*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,  double*, double*, int*, int*);

SEXP pgraph(SEXP nfunct, SEXP isi, SEXP zd, SEXP nn, SEXP npoint, SEXP days, SEXP h, SEXP delta, SEXP  dmax, SEXP kmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17,*d18,*d19;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;

    SEXP ans = R_NilValue, xtau = R_NilValue, y = R_NilValue, kn = R_NilValue, xl = R_NilValue, xx = R_NilValue, ydev  = R_NilValue, ui = R_NilValue, ncn = R_NilValue, sui = R_NilValue, xp = R_NilValue, xrate = R_NilValue, dlt = R_NilValue, xtime = R_NilValue, yvar = R_NilValue, sigma = R_NilValue, k = R_NilValue, ier = R_NilValue;

    double *zxtau, *zy, *zxl, *zxx, *zydev, *zui, *zncn, *zsui, *zxp, *zxrate, *zdlt, *zxtime, *zyvar, *zsigma = NULL;
    int      *zkn, *zk, *zier = NULL;

    int nn1, nn2, npt, kmx;
    int i;

    i1 = INTEGER_POINTER(nfunct);
    i2 = INTEGER_POINTER(isi);
    d1 = NUMERIC_POINTER(zd);
    i3 = INTEGER_POINTER(nn);
    i4 = INTEGER_POINTER(npoint);
    d2 = NUMERIC_POINTER(days);
    d3 = NUMERIC_POINTER(h);
    d4 = NUMERIC_POINTER(delta);
    d5 = NUMERIC_POINTER(dmax);
    i5 = INTEGER_POINTER(kmax);

    nn2 = 2*(*i3);
    nn1 = *i3 -1;
    npt = *i4;
    kmx = *i5;

    PROTECT(ans = allocVector(VECSXP, 17));
    SET_VECTOR_ELT(ans, 0, xtau = allocVector(REALSXP, nn2));
    SET_VECTOR_ELT(ans, 1, y = allocVector(REALSXP, nn2));
    SET_VECTOR_ELT(ans, 2, kn = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 3, xl = allocVector(REALSXP, nn1));
    SET_VECTOR_ELT(ans, 4, xx = allocVector(REALSXP, nn1*6));
    SET_VECTOR_ELT(ans, 5, ydev = allocVector(REALSXP, nn1));
    SET_VECTOR_ELT(ans, 6, ui = allocVector(REALSXP, nn1));
    SET_VECTOR_ELT(ans, 7, ncn = allocVector(REALSXP, nn1));
    SET_VECTOR_ELT(ans, 8, sui = allocVector(REALSXP, nn1));
    SET_VECTOR_ELT(ans, 9, xp = allocVector(REALSXP, 4));
    SET_VECTOR_ELT(ans, 10, xrate = allocVector(REALSXP, npt+1));
    SET_VECTOR_ELT(ans, 11, dlt = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 12, xtime = allocVector(REALSXP, kmx));
    SET_VECTOR_ELT(ans, 13, yvar = allocVector(REALSXP, kmx*5));
    SET_VECTOR_ELT(ans, 14, sigma = allocVector(REALSXP, kmx));
    SET_VECTOR_ELT(ans, 15, k = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 16, ier = allocVector(INTSXP, 1));

    d6 = NUMERIC_POINTER(xtau);
    d7 = NUMERIC_POINTER(y);
    i6 = INTEGER_POINTER(kn);
    d8 = NUMERIC_POINTER(xl);
    d9 = NUMERIC_POINTER(xx);
    d10 = NUMERIC_POINTER(ydev);
    d11 = NUMERIC_POINTER(ui);
    d12 = NUMERIC_POINTER(ncn);
    d13 = NUMERIC_POINTER(sui);
    d14 = NUMERIC_POINTER(xp);
    d15 = NUMERIC_POINTER(xrate);
    d16 = NUMERIC_POINTER(dlt);
    d17 = NUMERIC_POINTER(xtime);
    d18 = NUMERIC_POINTER(yvar);
    d19 = NUMERIC_POINTER(sigma);
    i7 = INTEGER_POINTER(k);
    i8 = INTEGER_POINTER(ier);

    F77_CALL(pgraphf) (i1,i2,d1,i3,i4,d2,d3,d4,d5,i5,d6,d7,i6,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,i7,i8);

    zxtau = REAL(xtau);
    zy = REAL(y);
    zkn = INTEGER(kn);
    zxl = REAL(xl);
    zxx = REAL(xx);
    zydev = REAL(ydev);
    zui = REAL(ui);
    zncn = REAL(ncn);
    zsui = REAL(sui);
    zxp = REAL(xp);
    zxrate = REAL(xrate);
    zdlt = REAL(dlt);
    zxtime = REAL(xtime);
    zyvar = REAL(yvar);
    zsigma = REAL(sigma);
    zk = INTEGER(k);
    zier = INTEGER(ier);

    for(i=0; i<nn2; i++) zxtau[i] = d6[i];
    for(i=0; i<nn2; i++) zy[i] = d7[i];
    *zkn = *i6;
    for(i=0; i<nn1; i++) zxl[i] = d8[i];
    for(i=0; i<nn1*6; i++) zxx[i] = d9[i];
    for(i=0; i<nn1; i++) zydev[i] = d10[i];
    for(i=0; i<nn1; i++) zui[i] = d11[i];
    for(i=0; i<nn1; i++) zncn[i] = d12[i];
    for(i=0; i<nn1; i++) zsui[i] = d13[i];
    for(i=0; i<4; i++) zxp[i] = d14[i];
    for(i=0; i<npt+1; i++) zxrate[i] = d15[i];
    *zdlt = *d16;
    for(i=0; i<kmx; i++) zxtime[i] = d17[i];
    for(i=0; i<kmx*5; i++) zyvar[i] = d18[i];
    for(i=0; i<kmx; i++) zsigma[i] = d19[i];
    *zk = *i7;
    *zier = *i8;

    UNPROTECT(1);

    return ans;
}
