#include "coin_common.h"

/* Helpers.c */
extern SEXP R_ExpectationCovarianceStatistic
(
    SEXP x,
    SEXP y,
    SEXP weights,
    SEXP subset,
    SEXP block,
    SEXP varonly,
    SEXP tol
);

extern SEXP R_PermutedLinearStatistic
(
    SEXP x,
    SEXP y,
    SEXP weights,
    SEXP subset,
    SEXP block,
    SEXP nresample
);

extern SEXP R_quadform
(
    SEXP linstat,
    SEXP expect,
    SEXP MPinv_sym
);

extern SEXP R_kronecker
(
    SEXP A,
    SEXP B
);

extern SEXP R_MPinv_sym
(
    SEXP x,
    SEXP n,
    SEXP tol
);

extern SEXP R_unpack_sym
(
    SEXP x,
    SEXP names,
    SEXP diagonal
);

extern SEXP R_maxstattrafo
(
    SEXP x,
    SEXP cutpoints
);

extern SEXP R_outersum
(
    SEXP A,
    SEXP B
);


/* StreitbergRoehmel.c */
extern SEXP R_cpermdist2
(
    SEXP score_b,
    SEXP m_a
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
