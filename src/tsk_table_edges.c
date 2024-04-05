#include <R.h>
#include <Rinternals.h>
#include <tskit.h>


SEXP
C_treeseq_edges(SEXP treeseq)
{
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);
    tsk_edge_table_t *edges = &ts->tables->edges;
    int i;
    int n = (int)edges->num_rows;
    SEXP df = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP rownames = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(colnames, 0, Rf_mkChar("edge_id"));
    SET_STRING_ELT(colnames, 1, Rf_mkChar("left"));
    SET_STRING_ELT(colnames, 2, Rf_mkChar("right"));
    SET_STRING_ELT(colnames, 3, Rf_mkChar("parent_id"));
    SET_STRING_ELT(colnames, 4, Rf_mkChar("child_id"));
    // edge id
    SET_VECTOR_ELT(df, 0, Rf_allocVector(INTSXP, n));
    // left
    SET_VECTOR_ELT(df, 1, Rf_allocVector(REALSXP, n));
    // right
    SET_VECTOR_ELT(df, 2, Rf_allocVector(REALSXP, n));
    // parent id
    SET_VECTOR_ELT(df, 3, Rf_allocVector(INTSXP, n));
    // child id
    SET_VECTOR_ELT(df, 4, Rf_allocVector(INTSXP, n));
    for (i = 0; i < n; ++i)
    {
        INTEGER(VECTOR_ELT(df, 0))[i] = i;
        REAL(VECTOR_ELT(df, 1))[i] = edges->left[i];
        REAL(VECTOR_ELT(df, 2))[i] = edges->right[i];
        INTEGER(VECTOR_ELT(df, 3))[i] = (int)edges->parent[i];
        INTEGER(VECTOR_ELT(df, 4))[i] = (int)edges->child[i];
        INTEGER(rownames)[i] = i+1;
    }
    Rf_setAttrib(df, R_NamesSymbol, colnames);
    Rf_setAttrib(df, R_ClassSymbol, Rf_mkString("data.frame"));
    Rf_setAttrib(df, Rf_mkString("row.names"), rownames);
    UNPROTECT(3);
    return df;
}
