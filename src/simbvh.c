#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(simbvhf)( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, int*, int*, double*, int*, int*);

SEXP simbvh(SEXP kxx, SEXP kxy, SEXP kxz, SEXP kyx, SEXP kyy, SEXP kyz, SEXP t, SEXP coxx, SEXP coxy, SEXP coyx, SEXP coyy, SEXP axx, SEXP axy, SEXP axz, SEXP ayx, SEXP ayy, SEXP ayz, SEXP ptxmax, SEXP ptymax, SEXP kmax, SEXP nnmax, SEXP mmmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11;

    SEXP ans = R_NilValue, x = R_NilValue, y = R_NilValue, ii1 = R_NilValue, jj1 = R_NilValue, err = R_NilValue;

    double *xx, *yy, *zerr = NULL;
    int     *ii, *jj = NULL;
    int nnn, mmm;
    int i;

    i1 = INTEGER_POINTER(kxx);
    i2 = INTEGER_POINTER(kxy);
    i3 = INTEGER_POINTER(kxz);
    i4 = INTEGER_POINTER(kyx);
    i5 = INTEGER_POINTER(kyy);
    i6 = INTEGER_POINTER(kyz);
    d1 = NUMERIC_POINTER(t);
    d2 = NUMERIC_POINTER(coxx);
    d3 = NUMERIC_POINTER(coxy);
    d4 = NUMERIC_POINTER(coyx);
    d5 = NUMERIC_POINTER(coyy);
    d6 = NUMERIC_POINTER(axx);
    d7 = NUMERIC_POINTER(axy);
    d8 = NUMERIC_POINTER(axz);
    d9 = NUMERIC_POINTER(ayx);
    d10 = NUMERIC_POINTER(ayy);
    d11 = NUMERIC_POINTER(ayz);
    d12 = NUMERIC_POINTER(ptxmax);
    d13 = NUMERIC_POINTER(ptymax);
    i7 = INTEGER_POINTER(kmax);
     i10 = INTEGER_POINTER(nnmax);
     i11 = INTEGER_POINTER(mmmax);

    nnn = *i10;
    mmm = *i11;
    PROTECT(ans = allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ans, 0, x = allocVector(REALSXP, nnn));
    SET_VECTOR_ELT(ans, 1, y = allocVector(REALSXP, mmm));
    SET_VECTOR_ELT(ans, 2, ii1 = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 3, jj1 = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 4, err = allocVector(REALSXP, 1));

    d14 = NUMERIC_POINTER(x);
    d15 = NUMERIC_POINTER(y);
    i8 = INTEGER_POINTER(ii1);
    i9 = INTEGER_POINTER(jj1);
    d16 = NUMERIC_POINTER(err);

    F77_CALL(simbvhf) (i1,i2,i3,i4,i5,i6,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,i7,d14,d15,i8,i9,d16,i10,i11);

    xx = REAL(x);
    yy = REAL(y);
    ii = INTEGER(ii1);
    jj = INTEGER(jj1);
    zerr = REAL(err);

    for(i=0; i<nnn; i++) xx[i] = d14[i];
    for(i=0; i<mmm; i++) yy[i] = d15[i];
    *ii = *i8;
    *jj = *i9;
    *zerr = *d16;

    UNPROTECT(1);

    return ans;
}
