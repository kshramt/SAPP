#include <R.h>
#include <Rdefines.h>
#include "sapp.h"

extern void F77_NAME(linsimf)(int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, double*, int*, double*, int*, int*, double*);

SEXP linsim(SEXP kxx, SEXP kxy, SEXP kxz, SEXP t, SEXP c, SEXP d, SEXP ax, SEXP ay, SEXP at, SEXP yy, SEXP mm, SEXP ptmax, SEXP kmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7;

    SEXP ans = R_NilValue, xx = R_NilValue, ii1= R_NilValue, jj1= R_NilValue, err = R_NilValue;

    double *zxx, *zerr = NULL;
    int      *ii, *jj = NULL;

    int mm2;
    int i;

    i1 = INTEGER_POINTER(kxx);
    i2 = INTEGER_POINTER(kxy);
    i3 = INTEGER_POINTER(kxz);
    d1 = NUMERIC_POINTER(t);
    d2 = NUMERIC_POINTER(c);
    d3 = NUMERIC_POINTER(d);
    d4 = NUMERIC_POINTER(ax);
    d5 = NUMERIC_POINTER(ay);
    d6 = NUMERIC_POINTER(at);
    d7 = NUMERIC_POINTER(yy);
    i4 = INTEGER_POINTER(mm);
    d8 = NUMERIC_POINTER(ptmax);
    i5 = INTEGER_POINTER(kmax);

    mm2 = (*i4) * 2;

    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, xx = allocVector(REALSXP, mm2));
    SET_VECTOR_ELT(ans, 1, ii1 = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 2, jj1 = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 3, err = allocVector(REALSXP, 1));

    d9 = NUMERIC_POINTER(xx);
    i6 = INTEGER_POINTER(ii1);
    i7 = INTEGER_POINTER(jj1);
    d10 = NUMERIC_POINTER(err);

    F77_CALL(linsimf) (i1,i2,i3,d1,d2,d3,d4,d5,d6,d7,i4,d8,i5,d9,i6,i7,d10);

    zxx = REAL(xx);
    ii = INTEGER(ii1);
    jj = INTEGER(jj1);
    zerr = REAL(err);

    for(i=0; i<mm2; i++) zxx[i] = d9[i];
    *ii = *i6;
    *jj = *i7;
    *zerr = *d10;

    UNPROTECT(1);

    return ans;
}
