#include <R.h>
#include <Rinternals.h>

#include <tskit.h>
#include "error.h"

void 
tsx_treeseq_exptr_free(SEXP treeseq) 
{
  tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);
  if (ts != NULL) {
    tsk_treeseq_free(ts);
    free(ts);
    ts = NULL;
  }
}

void 
tsx_tree_exptr_free(SEXP tr) 
{
  tsk_tree_t *tree = (tsk_tree_t *)R_ExternalPtrAddr(tr);
  if (tree != NULL) {
    tsk_tree_free(tree);
    free(tree);
    tree = NULL;
  }
}

SEXP 
C_treeseq_load(SEXP filename) 
{
  tsk_treeseq_t *ts = malloc(sizeof(*ts));
  tsk_tree_t *tree = malloc(sizeof(*tree));
  if (!ts || !tree)
    TSX_ERROR("memory allocation failed");
  
  int rc = tsk_treeseq_load(ts, CHAR(STRING_ELT(filename, 0)), 0);
  check_tsk_error(rc);
  
  rc = tsk_tree_init(tree, ts, TSK_SAMPLE_LISTS);
  check_tsk_error(rc);
  
  SEXP exptr1 = PROTECT(R_MakeExternalPtr(ts, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(exptr1, tsx_treeseq_exptr_free, TRUE);
  SEXP exptr2 = PROTECT(R_MakeExternalPtr(tree, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(exptr2, tsx_tree_exptr_free, TRUE);
  UNPROTECT(2);
  
  return Rf_list2(exptr1, exptr2);
}
