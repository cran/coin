
#include "coin_common.h"

/* Classes.c */
extern SEXP coin_init
(
    void
);


/* Helpers.c */
extern SEXP R_blocksetup
(
    SEXP block
);

extern SEXP R_blockperm
(
    SEXP block
);

extern SEXP R_MonteCarloIndependenceTest
(
    SEXP x,
    SEXP y,
    SEXP block,
    SEXP B
);

extern SEXP R_maxstattrafo
(
    SEXP x,
    SEXP cutpoints
);


/* LinearStatistic.c */
extern SEXP R_kronecker
(
    SEXP A,
    SEXP B
);

extern SEXP R_outersum
(
    SEXP A,
    SEXP B
);

extern SEXP R_ExpectCovarInfluence
(
    SEXP y,
    SEXP weights
);

extern SEXP R_ExpectCovarLinearStatistic
(
    SEXP x,
    SEXP y,
    SEXP weights,
    SEXP expcovinf
);

extern SEXP R_LinearStatistic
(
    SEXP x,
    SEXP y,
    SEXP weights
);

extern SEXP R_PermutedLinearStatistic
(
    SEXP x,
    SEXP y,
    SEXP indx,
    SEXP perm
);


/* StreitbergRoehmel.c */
extern SEXP R_cpermdist2
(
    SEXP score_a,
    SEXP score_b,
    SEXP m_a,
    SEXP m_b,
    SEXP retProb
);

extern SEXP R_cpermdist1
(
    SEXP scores
);


/* vandeWiel.c */
extern SEXP R_split_up_2sample
(
    SEXP scores,
    SEXP m,
    SEXP obs,
    SEXP tol
);
