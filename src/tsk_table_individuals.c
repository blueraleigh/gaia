#include <R.h>
#include <Rinternals.h>

#include <tskit.h>
#include "error.h"


SEXP
C_treeseq_individuals(SEXP treeseq)
{
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);
    tsk_individual_table_t *ind = &ts->tables->individuals;
    int i;
    int j;
    int m;
    int n = (int)ind->num_rows;
    SEXP df = PROTECT(Rf_allocMatrix(VECSXP, n, 4));
    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(colnames, 0, Rf_mkChar("individual_id"));
    SET_STRING_ELT(colnames, 1, Rf_mkChar("location"));
    SET_STRING_ELT(colnames, 2, Rf_mkChar("parents"));
    SET_STRING_ELT(colnames, 3, Rf_mkChar("metadata"));
    Rf_setAttrib(df, R_DimNamesSymbol, Rf_list2(R_NilValue, colnames));
    for (i = 0; i < n; ++i)
    {
        SET_VECTOR_ELT(df, i+0*n, Rf_ScalarInteger(i));
        m = ind->location_offset[i+1] - ind->location_offset[i];
        SET_VECTOR_ELT(df, i+1*n, Rf_allocVector(REALSXP, m));
        for (j = 0; j < m; ++j)
            REAL(VECTOR_ELT(df, i+1*n))[j] = ind->location[ind->location_offset[i] + j];
        m = ind->parents_offset[i+1] - ind->parents_offset[i];
        SET_VECTOR_ELT(df, i+2*n, Rf_allocVector(INTSXP, m));
        for (j = 0; j < m; ++j)
            INTEGER(VECTOR_ELT(df, i+2*n))[j] = ind->parents[ind->parents_offset[i] + j];
        m = ind->metadata_offset[i+1] - ind->metadata_offset[i];
        SET_VECTOR_ELT(df, i+3*n, Rf_allocVector(RAWSXP, m));
        memcpy(RAW(VECTOR_ELT(df, i+3*n)),
            ind->metadata + ind->metadata_offset[i], m);
    }
    UNPROTECT(2);
    return df;
}
