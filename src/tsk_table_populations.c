#include <R.h>
#include <Rinternals.h>

#include <tskit.h>
#include "error.h"


SEXP
C_treeseq_populations(SEXP treeseq)
{
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);
    tsk_population_table_t *pop = &ts->tables->populations;
    int i;
    int m;
    int n = (int)pop->num_rows;
    SEXP df = PROTECT(Rf_allocMatrix(VECSXP, n, 2));
    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(colnames, 0, Rf_mkChar("population_id"));
    SET_STRING_ELT(colnames, 1, Rf_mkChar("metadata"));
    Rf_setAttrib(df, R_DimNamesSymbol, Rf_list2(R_NilValue, colnames));
    for (i = 0; i < n; ++i)
    {
        SET_VECTOR_ELT(df, i+0*n, Rf_ScalarInteger(i));
        m = pop->metadata_offset[i+1] - pop->metadata_offset[i];
        SET_VECTOR_ELT(df, i+1*n, Rf_allocVector(RAWSXP, m));
        memcpy(RAW(VECTOR_ELT(df, i+1*n)),
            pop->metadata + pop->metadata_offset[i], m);
    }
    UNPROTECT(2);
    return df;
}
