#include "coin.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef callMethods[] = {
    CALLDEF(coin_init, 0),
    CALLDEF(R_blocksetup, 1),
    CALLDEF(R_blockperm, 1),
    CALLDEF(R_MonteCarloIndependenceTest, 4),
    CALLDEF(R_maxstattrafo, 2),
    CALLDEF(R_kronecker, 2),
    CALLDEF(R_outersum, 2),
    CALLDEF(R_ExpectCovarInfluence, 2),
    CALLDEF(R_ExpectCovarLinearStatistic, 4),
    CALLDEF(R_LinearStatistic, 3),
    CALLDEF(R_PermutedLinearStatistic, 4),
    CALLDEF(R_cpermdist2, 5),
    CALLDEF(R_cpermdist1, 1),
    CALLDEF(R_split_up_2sample, 4),
    {NULL, NULL, 0}
};

void attribute_visible R_init_coin(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
