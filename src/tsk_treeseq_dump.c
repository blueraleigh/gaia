#include <R.h>
#include <Rinternals.h>

#include <tskit.h>
#include "error.h"

SEXP
C_treeseq_write(SEXP treeseq, SEXP filename)
{
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);
    check_tsk_error(tsk_treeseq_dump(ts, CHAR(STRING_ELT(filename, 0)), 0));
    return R_NilValue;
}
