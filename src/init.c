#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>


#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

SEXP C_treeseq_load(SEXP);
SEXP C_treeseq_write(SEXP,SEXP);
SEXP C_treeseq_to_phylo(SEXP);
SEXP C_treeseq_simplify(SEXP,SEXP,SEXP);
SEXP C_treeseq_sample(SEXP);
SEXP C_treeseq_nodes(SEXP);
SEXP C_treeseq_edges(SEXP);
SEXP C_treeseq_individuals(SEXP);
SEXP C_treeseq_intervals(SEXP);
SEXP C_treeseq_populations(SEXP);
SEXP C_treeseq_indexes(SEXP);
SEXP C_treeseq_drop_edges(SEXP,SEXP,SEXP);
SEXP C_treeseq_fixup_edge_ids(SEXP, SEXP, SEXP);

SEXP C_treeseq_discrete_mpr(SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_treeseq_discrete_mpr_minimize(SEXP,SEXP);
SEXP C_treeseq_discrete_mpr_edge_history(SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_treeseq_discrete_mpr_ancestry_flux(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,
    SEXP,SEXP);
SEXP C_treeseq_discrete_mpr_ancestry(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,
    SEXP,SEXP);

SEXP C_treeseq_quadratic_mpr(SEXP,SEXP,SEXP);
SEXP C_treeseq_quadratic_mpr_minimize(SEXP);
SEXP C_treeseq_quadratic_mpr_minimize_discrete(SEXP,SEXP);

SEXP C_treeseq_linear_mpr(SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_treeseq_linear_mpr_minimize(SEXP);
SEXP C_treeseq_linear_mpr_minimize_discrete(SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(C_treeseq_load, 1),
    CALLDEF(C_treeseq_write, 2),
    CALLDEF(C_treeseq_to_phylo, 1),
    CALLDEF(C_treeseq_simplify, 3),
    CALLDEF(C_treeseq_sample, 1),
    CALLDEF(C_treeseq_nodes, 1),
    CALLDEF(C_treeseq_edges, 1),
    CALLDEF(C_treeseq_individuals, 1),
    CALLDEF(C_treeseq_intervals, 1),
    CALLDEF(C_treeseq_populations, 1),
    CALLDEF(C_treeseq_indexes, 1),
    CALLDEF(C_treeseq_drop_edges, 3),
    CALLDEF(C_treeseq_fixup_edge_ids, 3),

    CALLDEF(C_treeseq_discrete_mpr, 5),
    CALLDEF(C_treeseq_discrete_mpr_minimize, 2),
    CALLDEF(C_treeseq_discrete_mpr_edge_history, 5),
    CALLDEF(C_treeseq_discrete_mpr_ancestry_flux, 9),
    CALLDEF(C_treeseq_discrete_mpr_ancestry, 10),

    CALLDEF(C_treeseq_quadratic_mpr, 3),
    CALLDEF(C_treeseq_quadratic_mpr_minimize, 1),
    CALLDEF(C_treeseq_quadratic_mpr_minimize_discrete, 2),

    CALLDEF(C_treeseq_linear_mpr, 5),
    CALLDEF(C_treeseq_linear_mpr_minimize, 1),
    CALLDEF(C_treeseq_linear_mpr_minimize_discrete, 2),

    {NULL, NULL, 0}
};


void attribute_visible R_init_gaia(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
