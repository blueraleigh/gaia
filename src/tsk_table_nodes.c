#include <R.h>
#include <Rinternals.h>

#include <tskit.h>
#include "error.h"


SEXP
C_treeseq_nodes(SEXP treeseq)
{
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);
    tsk_node_table_t *nodes = &ts->tables->nodes;
    int i;
    int m;
    int n = (int)nodes->num_rows;
    SEXP df = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP rownames = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 6));
    SET_STRING_ELT(colnames, 0, Rf_mkChar("node_id"));
    SET_STRING_ELT(colnames, 1, Rf_mkChar("is_sample"));
    SET_STRING_ELT(colnames, 2, Rf_mkChar("time"));
    SET_STRING_ELT(colnames, 3, Rf_mkChar("population_id"));
    SET_STRING_ELT(colnames, 4, Rf_mkChar("individual_id"));
    SET_STRING_ELT(colnames, 5, Rf_mkChar("metadata"));
    // node id
    SET_VECTOR_ELT(df, 0, Rf_allocVector(INTSXP, n));
    // is_sample
    SET_VECTOR_ELT(df, 1, Rf_allocVector(INTSXP, n));
    // time
    SET_VECTOR_ELT(df, 2, Rf_allocVector(REALSXP, n));
    // population id
    SET_VECTOR_ELT(df, 3, Rf_allocVector(INTSXP, n));
    // individual id
    SET_VECTOR_ELT(df, 4, Rf_allocVector(INTSXP, n));
    // metadata
    SET_VECTOR_ELT(df, 5, Rf_allocVector(VECSXP, n));
    for (i = 0; i < n; ++i)
    {
        INTEGER(VECTOR_ELT(df, 0))[i] = i;
        INTEGER(VECTOR_ELT(df, 1))[i] = nodes->flags[i] & TSK_NODE_IS_SAMPLE;
        REAL(VECTOR_ELT(df, 2))[i] = nodes->time[i];
        INTEGER(VECTOR_ELT(df, 3))[i] = (int)nodes->population[i];
        INTEGER(VECTOR_ELT(df, 4))[i] = (int)nodes->individual[i];
        m = nodes->metadata_offset[i+1] - nodes->metadata_offset[i];
        SET_VECTOR_ELT(VECTOR_ELT(df, 5), i, Rf_allocVector(RAWSXP, m));
        memcpy(RAW(VECTOR_ELT(VECTOR_ELT(df, 5), i)),
            nodes->metadata + nodes->metadata_offset[i], m);
        INTEGER(rownames)[i] = i+1;
    }
    Rf_setAttrib(df, R_NamesSymbol, colnames);
    Rf_setAttrib(df, R_ClassSymbol, Rf_mkString("data.frame"));
    Rf_setAttrib(df, Rf_mkString("row.names"), rownames);
    UNPROTECT(3);
    return df;
}
