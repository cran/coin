\name{SurvivalTests}
\alias{surv_test}
\alias{logrank_test}
\alias{logrank_test.formula}
\alias{logrank_test.IndependenceProblem}
\concept{Logrank test}
\concept{Gehan-Breslow test}
\concept{Tarone-Ware test}
\concept{Prentice test}
\concept{Prentice-Marek test}
\concept{Andersen-Borgan-Gill-Keiding test}
\concept{Fleming-Harrington test}
\concept{Gaugler-Kim-Liao test}
\concept{Self test}
\encoding{UTF-8}
\title{Two- and \eqn{K}-Sample Tests for Censored Data}
\description{
  Testing the equality of the survival distributions in two or more independent
  groups.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\method{logrank_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{logrank_test}{IndependenceProblem}(object, ties.method = c("mid-ranks", "Hothorn-Lausen",
                                     "average-scores"),
             type = c("logrank", "Gehan-Breslow", "Tarone-Ware", "Prentice",
                      "Prentice-Marek", "Andersen-Borgan-Gill-Keiding",
                      "Fleming-Harrington", "Gaugler-Kim-Liao", "Self"),
             rho = NULL, gamma = NULL, \dots)
}
\arguments{
  \item{formula}{
    a formula of the form \code{y ~ x | block} where \code{y} is a survival
    object (see \code{\link[survival]{Surv}} in package \pkg{survival}),
    \code{x} is a factor and \code{block} is an optional factor for
    stratification.
  }
  \item{data}{
    an optional data frame containing the variables in the model formula.
  }
  \item{subset}{
    an optional vector specifying a subset of observations to be used.  Defaults
    to \code{NULL}.
  }
  \item{weights}{
    an optional formula of the form \code{~ w} defining integer valued case
    weights for each observation.  Defaults to \code{NULL}, implying equal
    weight for all observations.
  }
  \item{object}{
    an object inheriting from class \code{"\linkS4class{IndependenceProblem}"}.
  }
  \item{ties.method}{
    a character, the method used to handle ties: the score generating function
    either uses mid-ranks (\code{"mid-ranks"}, default), the Hothorn-Lausen
    method (\code{"Hothorn-Lausen"}) or averages the scores of randomly broken
    ties (\code{"average-scores"}); see \sQuote{Details}.
  }
  \item{type}{
    a character, the type of test: either \code{"logrank"} (default),
    \code{"Gehan-Breslow"}, \code{"Tarone-Ware"}, \code{"Prentice"},
    \code{"Prentice-Marek"}, \code{"Andersen-Borgan-Gill-Keiding"},
    \code{"Fleming-Harrington"}, \code{"Gaugler-Kim-Liao"} or \code{"Self"}; see
    \sQuote{Details}.
  }
  \item{rho}{
    a numeric, the \eqn{\rho} constant when \code{type} is \code{"Tarone-Ware"},
    \code{"Fleming-Harrington"}, \code{"Gaugler-Kim-Liao"} or \code{"Self"}; see
    \sQuote{Details}.  Defaults to \code{NULL}, implying \code{0.5} for
    \code{type = "Tarone-Ware"} and \code{0} otherwise.
  }
  \item{gamma}{
    a numeric, the \eqn{\gamma} constant when \code{type} is
    \code{"Fleming-Harrington"}, \code{"Gaugler-Kim-Liao"} or \code{"Self"}; see
    \sQuote{Details}.  Defaults to \code{NULL}, implying \code{0}.
  }
  \item{\dots}{
    further arguments to be passed to \code{\link{independence_test}}.
  }
}
\details{
  \code{logrank_test} provides the weighted logrank test reformulated as a
  linear rank test.  The family of weighted logrank tests encompasses a large
  collection of tests commonly used in the analysis of survival data including,
  but not limited to, the standard (unweighted) logrank test, the Gehan-Breslow
  test, the Tarone-Ware class of tests, the Prentice test, the Prentice-Marek
  test, the Andersen-Borgan-Gill-Keiding test, the Fleming-Harrington class of
  tests and the Self class of tests.  A general description of these methods is
  given by Klein and Moeschberger (2003, Ch. 7).  See \enc{Letón}{Leton} and
  Zuluaga (2001) for the linear rank test formulation.

  The null hypothesis of equality, or conditional equality given \code{block},
  of the survival distribution of \code{y} in the groups defined by \code{x} is
  tested.  In the two-sample case, the two-sided null hypothesis is \eqn{H_0\!:
    \theta = 1}{H_0: theta = 1}, where \eqn{\theta = \lambda_2 / \lambda_1}
  and \eqn{\lambda_s} is the hazard rate in the \eqn{s}th sample.  In case
  \code{alternative = "less"}, the null hypothesis is \eqn{H_0\!: \theta \ge
    1}{H_0: theta >= 1}, i.e., the survival is lower in sample 1 than in sample
  2.  When \code{alternative = "greater"}, the null hypothesis is \eqn{H_0\!:
    \theta \le 1}{H_0: theta <= 1}, i.e., the survival is higher in sample 1
  than in sample 2.

  If \code{x} is an ordered factor, the default scores, \code{1:nlevels(x)}, can
  be altered using the \code{scores} argument (see
  \code{\link{independence_test}}); this argument can also be used to coerce
  nominal factors to class \code{"ordered"}.  In this case, a linear-by-linear
  association test is computed and the direction of the alternative hypothesis
  can be specified using the \code{alternative} argument.  This type of
  extension of the standard logrank test was given by Tarone (1975) and later
  generalized to general weights by Tarone and Ware (1977).

  Let \eqn{(t_i, \delta_i)}, \eqn{i = 1, 2, \ldots, n}, represent a
  right-censored random sample of size \eqn{n}, where \eqn{t_i} is the observed
  survival time and \eqn{\delta_i} is the status indicator (\eqn{\delta_i} is 0
  for right-censored observations and 1 otherwise).  To allow for ties in the
  data, let \eqn{t_{(1)} < t_{(2)} < \cdots < t_{(m)}}{t_(1) < t_(2) < \dots <
    t_(m)} represent the \eqn{m}, \eqn{m \le n}, ordered distinct event times.
  At time \eqn{t_{(k)}}{t_(k)}, \eqn{k = 1, 2, \ldots, m}, the number of events
  and the number of subjects at risk are given by \eqn{d_k = \sum_{i = 1}^n
    I\!\left(t_i = t_{(k)}\,|\, \delta_i = 1\right)}{d_k = sum(i = 1, \dots, n)
    I(t_i = t_(k) | delta_i = 1)} and \eqn{n_k = n - r_k}, respectively, where
  \eqn{r_k} depends on the ties handling method.

  Three different methods of handling ties are available using
  \code{ties.method}: mid-ranks (\code{"mid-ranks"}, default), the
  Hothorn-Lausen method (\code{"Hothorn-Lausen"}) and average-scores
  (\code{"average-scores"}).  The first and last method are discussed and
  contrasted by Callaert (2003), whereas the second method is defined in Hothorn
  and Lausen (2003).  The mid-ranks method leads to
  \deqn{
    r_k = \sum_{i = 1}^n I\!\left(t_i < t_{(k)}\right)
  }{
    r_k = sum(i = 1, \dots, n) I(t_i < t_(k))
  }
  whereas the Hothorn-Lausen method uses
  \deqn{
    r_k = \sum_{i = 1}^n I\!\left(t_i \le t_{(k)}\right) - 1.
  }{
    r_k = sum(i = 1, \dots, n) I(t_i <= t_(k)) - 1.
  }
  The scores assigned to right-censored and uncensored observations at the
  \eqn{k}th event time are given by
  \deqn{
    C_k = \sum_{j = 1}^k w_j \frac{d_j}{n_j}
    \quad \mbox{and} \quad
    c_k = C_k - w_k,
  }{
    C_k = sum(j=1,\dots,k) w_j * (d_j / n_j)   and   c_k = C_k - w_k,
  }
  respectively, where \eqn{w} is the logrank weight.  For the average-scores
  method, used by, e.g., the software package StatXact, the \eqn{d_k} events
  observed at the \eqn{k}th event time are arbitrarily ordered by assigning them
  distinct values \eqn{t_{(k_l)}}{t_(k_l)}, \eqn{l = 1, 2, \ldots, d_k},
  infinitesimally to the left of \eqn{t_{(k)}}{t_(k)}.  Then scores
  \eqn{C_{k_l}}{C_k_l} and \eqn{c_{k_l}}{c_k_l} are computed as indicated above,
  effectively assuming that no event times are tied.  The scores \eqn{C_k} and
  \eqn{c_k} are assigned the average of the scores \eqn{C_{k_l}}{C_k_l} and
  \eqn{c_{k_l}}{c_k_l} respectively.  It then follows that the score for the
  \eqn{i}th subject is
  \deqn{
    a_i = \left\{
      \begin{array}{ll}
        C_{k'} & \mbox{if } \delta_i = 0 \\
        c_{k'} & \mbox{otherwise}
      \end{array}
    \right.
  }{
                           C_k'   if delta_i = 0
                     a_i =
                           c_k'   otherwise
  }
  where \eqn{k' = \max \{k: t_i \ge t_{(k)}\}}{k' = max\{k : t_i >= t_(k)\}}.

  The \code{type} argument allows for a choice between some of the most
  well-known members of the family of weighted logrank tests, each corresponding
  to a particular weight function.  The standard logrank test (\code{"logrank"},
  default) was suggested by Mantel (1966), Peto and Peto (1972) and Cox (1972)
  and has \eqn{w_k = 1}.  The Gehan-Breslow test (\code{"Gehan-Breslow"})
  proposed by Gehan (1965) and later extended to \eqn{K} samples by Breslow
  (1970) is a generalization of the Wilcoxon rank-sum test, where \eqn{w_k =
    n_k}.  The Tarone-Ware class of tests (\code{"Tarone-Ware"}) discussed by
  Tarone and Ware (1977) has \eqn{w_k = n_k^\rho}, where \eqn{\rho} is a
  constant; \eqn{\rho = 0.5} (default) was suggested by Tarone and Ware (1977),
  but note that \eqn{\rho = 0} and \eqn{\rho = 1} lead to the the standard
  logrank test and Gehan-Breslow test respectively.  The Prentice test
  (\code{"Prentice"}) is another generalization of the Wilcoxon rank-sum test
  proposed by Prentice (1978), where
  \deqn{
    w_k = \prod_{j = 1}^k \frac{n_j}{n_j + d_j}.
  }{
    w_k = prod(j = 1, \dots, k) n_j / n_j + d_j).
  }
  The Prentice-Marek test (\code{"Prentice-Marek"}) is yet another
  generalization of the Wilcoxon rank-sum test discussed by Prentice and Marek
  (1979), with
  \deqn{
    w_k = \tilde{S}_k = \prod_{j = 1}^k \frac{n_j + 1 - d_j}{n_j + 1}.
  }{
    w_k = Stilde_k = prod(j = 1, \dots, k) (n_j + 1 - d_j) / (n_j + 1).
  }
  The Andersen-Borgan-Gill-Keiding test (\code{"Andersen-Borgan-Gill-Keiding"})
  suggested by Andersen \emph{et al.} (1982) is a modified version of the
  Prentice-Marek test using
  \deqn{
    w_k = \frac{n_k}{n_k + 1}
            \prod_{j = 0}^{k - 1} \frac{n_j + 1 - d_j}{n_j + 1}
  }{
    w_k = (n_k / (n_k + 1)) prod(j = 0, \dots, k - 1) (n_j + 1 - d_j) / (n_j + 1)
  }
  where \eqn{n_0 \equiv n}{n_0 := n} and \eqn{d_0 \equiv 0}{d_0 := 0}.  The
  Fleming-Harrington class of tests (\code{"Fleming-Harrington"}) proposed by
  Fleming and Harrington (1991) uses \eqn{w_k = \hat{S}_k^\rho (1 -
    \hat{S}_k)^\gamma}{w_k = Shat_k^rho * (1 - Shat_k)^gamma}, where \eqn{\rho}
  and \eqn{\gamma} are constants and
  \deqn{
    \hat{S}_k = \prod_{j = 0}^{k - 1} \frac{n_j - d_j}{n_j},
    \quad
    \hat{S}_0 \equiv 1
  }{
    Shat_k = prod(j = 0, \dots, k - 1) (n_j - d_j) / n_j,   Shat_0 := 1
  }
  is the \emph{left-continuous} Kaplan-Meier estimator of the survival function;
  \eqn{\rho = 0} and \eqn{\gamma = 0} lead to the standard logrank test.  The
  Gaugler-Kim-Liao class of tests (\code{"Gaugler-Kim-Liao"}) discussed by
  Gaugler \emph{et al.} (2007) is a modified version of the Fleming-Harrington
  class of tests, replacing \eqn{\hat{S}_k}{Shat_k} with
  \eqn{\tilde{S}_k}{Stilde_k} so that \eqn{w_k = \tilde{S}_k^\rho (1 -
    \tilde{S}_k)^\gamma}{w_k = Stilde_k^rho * (1 - Stilde_k)^gamma}, where
  \eqn{\rho} and \eqn{\gamma} are constants; \eqn{\rho = 0} and \eqn{\gamma = 0}
  lead to the standard logrank test.  The Self class of tests (\code{"Self"})
  suggested by Self (1991) has \eqn{w_k = v_k^\rho (1 - v_k)^\gamma}{w_k =
    v_k^rho * (1 - v_k)^gamma}, where
  \deqn{
    v_k = \frac{1}{2} \frac{t_{(k-1)} + t_{(k)}}{t_{(m)}},
    \quad
    t_{(0)} \equiv 0
  }{
    v_k = 1 / 2 * (t_(k - 1) + t_(k)) / t_(m),   t_(0) := 1
  }
  is the standardized mid-point between the \eqn{(k - 1)}th and the \eqn{k}th
  event time.  (This is a slight generalization of Self's original proposal in
  order to allow for non-integer follow-up times.)  Again, \eqn{\rho} and
  \eqn{\gamma} are constants and \eqn{\rho = 0} and \eqn{\gamma = 0} lead to
  the standard logrank test.

  The conditional null distribution of the test statistic is used to obtain
  \eqn{p}-values and an asymptotic approximation of the exact distribution is
  used by default (\code{distribution = "asymptotic"}).  Alternatively, the
  distribution can be approximated via Monte Carlo resampling or computed
  exactly for univariate two-sample problems by setting \code{distribution} to
  \code{"approximate"} or \code{"exact"} respectively.  See
  \code{\link{asymptotic}}, \code{\link{approximate}} and \code{\link{exact}}
  for details.
}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
}
\note{
  Peto and Peto (1972) proposed the test statistic implemented in
  \code{logrank_test} and named it the \emph{logrank test}.  However, the
  Mantel-Cox test (Mantel, 1966; Cox, 1972), as implemented in
  \code{\link[survival]{survdiff}} (in package \pkg{survival}), is also known as
  the logrank test.  These tests are similar, but differ in the choice of
  probability model: the (Peto-Peto) logrank test uses the permutational
  variance, whereas the Mantel-Cox test is based on the hypergeometric variance.

  Combining \code{\link{independence_test}} or \code{\link{symmetry_test}} with
  \code{\link{logrank_trafo}} offers more flexibility than \code{logrank_test}
  and allows for, among other things, maximum-type versatile test procedures
  (e.g., Lee, 1996; see \sQuote{Examples}) and user-supplied logrank weights
  (see \code{\link{GTSG}} for tests against Weibull-type or crossing-curve
  alternatives).

  Starting with version 1.1-0, \code{logrank_test} replaced \code{surv_test}
  which was made \strong{defunct} in version 1.2-0.  Furthermore,
  \code{logrank_trafo} is now an increasing function for all choices of
  \code{ties.method}, implying that the test statistic has the same sign
  irrespective of the ties handling method.  Consequently, the sign of the test
  statistic will now be the opposite of what it was in earlier versions unless
  \code{ties.method = "average-scores"}.  (In versions of \pkg{coin} prior to
  1.1-0, \code{logrank_trafo} was a decreasing function when \code{ties.method}
  was other than \code{"average-scores"}.)

  Starting with version 1.2-0, mid-ranks and the Hothorn-Lausen method can no
  longer be specified with \code{ties.method = "logrank"} and
  \code{ties-method = "HL"} respectively.
}
\references{
  Andersen, P. K., Borgan, \enc{Ø}{O}., Gill, R. and Keiding, N.  (1982).
  Linear nonparametric tests for comparison of counting processes, with
  applications to censored survival data (with discussion).  \emph{International
    Statistical Review} \bold{50}(3), 219--258.  \doi{10.2307/1402489}

  Breslow, N.  (1970).  A generalized Kruskal-Wallis test for comparing \eqn{K}
  samples subject to unequal patterns of censorship.  \emph{Biometrika}
  \bold{57}(3), 579--594.  \doi{10.1093/biomet/57.3.579}

  Callaert, H.  (2003).  Comparing statistical software packages: The case of
  the logrank test in StatXact.  \emph{The American Statistician} \bold{57}(3),
  214--217.  \doi{10.1198/0003130031900}

  Cox, D. R. (1972).  Regression models and life-tables (with discussion).
  \emph{Journal of the Royal Statistical Society} B \bold{34}(2), 187--220.

  Fleming, T. R. and Harrington, D. P.  (1991).  \emph{Counting Processes and
    Survival Analysis}.  New York: John Wiley & Sons.

  Gaugler, T., Kim, D. and Liao, S.  (2007).  Comparing two survival time
  distributions: An investigation of several weight functions for the weighted
  logrank statistic.  \emph{Communications in Statistics -- Simulation and
    Computation} \bold{36}(2), 423--435.  \doi{10.1080/03610910601161272}

  Gehan, E. A.  (1965).  A generalized Wilcoxon test for comparing arbitrarily
  single-censored samples.  \emph{Biometrika} \bold{52}(1--2), 203--223.
  \doi{10.1093/biomet/52.1-2.203}

  Hothorn, T. and Lausen, B.  (2003).  On the exact distribution of maximally
  selected rank statistics.  \emph{Computational Statistics & Data Analysis}
  \bold{43}(2), 121--137.  \doi{10.1016/S0167-9473(02)00225-6}

  Klein, J. P. and Moeschberger, M. L.  (2003).  \emph{Survival Analysis:
    Techniques for Censored and Truncated Data}, Second Edition.  New York:
  Springer.

  Lee, J. W.  (1996).  Some versatile tests based on the simultaneous use of
  weighted log-rank statistics.  \emph{Biometrics} \bold{52}(2), 721--725.
  \doi{10.2307/2532911}

  \enc{Letón}{Leton}, E. and Zuluaga, P.  (2001).  Equivalence between score
  and weighted tests for survival curves.  \emph{Communications in Statistics --
    Theory and Methods} \bold{30}(4), 591--608.  \doi{10.1081/STA-100002138}

  Mantel, N.  (1966).  Evaluation of survival data and two new rank order
  statistics arising in its consideration.  \emph{Cancer Chemotherapy Reports}
  \bold{50}(3), 163--170.

  Peto, R. and Peto, J.  (1972).  Asymptotic efficient rank invariant test
  procedures (with discussion).  \emph{Journal of the Royal Statistical Society}
  A \bold{135}(2), 185--207.  \doi{10.2307/2344317}

  Prentice, R. L.  (1978).  Linear rank tests with right censored data.
  \emph{Biometrika} \bold{65}(1), 167--179.  \doi{10.1093/biomet/65.1.167}

  Prentice, R. L. and Marek, P.  (1979).  A qualitative discrepancy between
  censored data rank tests.  \emph{Biometrics} \bold{35}(4), 861--867.
  \doi{10.2307/2530120}

  Self, S. G.  (1991).  An adaptive weighted log-rank test with application to
  cancer prevention and screening trials.  \emph{Biometrics} \bold{47}(3),
  975--986.  \doi{10.2307/2532653}

  Tarone, R. E.  (1975).  Tests for trend in life table analysis.
  \emph{Biometrika} \bold{62}(3), 679--682.  \doi{10.1093/biomet/62.3.679}

  Tarone, R. E. and Ware, J.  (1977).  On distribution-free tests for equality
  of survival distributions.  \emph{Biometrika} \bold{64}(1), 156--160.
  \doi{10.1093/biomet/64.1.156}
}
\examples{
## Example data (Callaert, 2003, Tab. 1)
callaert <- data.frame(
    time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
    group = factor(rep(0:1, c(7, 8)))
)

## Logrank scores using mid-ranks (Callaert, 2003, Tab. 2)
with(callaert,
     logrank_trafo(Surv(time)))

## Asymptotic Mantel-Cox test (p = 0.0523)
survdiff(Surv(time) ~ group, data = callaert)

## Exact logrank test using mid-ranks (p = 0.0505)
logrank_test(Surv(time) ~ group, data = callaert, distribution = "exact")

## Exact logrank test using average-scores (p = 0.0468)
logrank_test(Surv(time) ~ group, data = callaert, distribution = "exact",
             ties.method = "average-scores")


## Lung cancer data (StatXact 9 manual, p. 213, Tab. 7.19)
lungcancer <- data.frame(
    time = c(257, 476, 355, 1779, 355,
             191, 563, 242, 285, 16, 16, 16, 257, 16),
    event = c(0, 0, 1, 1, 0,
              1, 1, 1, 1, 1, 1, 1, 1, 1),
    group = factor(rep(1:2, c(5, 9)),
                   labels = c("newdrug", "control"))
)

## Logrank scores using average-scores (StatXact 9 manual, p. 214)
with(lungcancer,
     logrank_trafo(Surv(time, event), ties.method = "average-scores"))

## Exact logrank test using average-scores (StatXact 9 manual, p. 215)
logrank_test(Surv(time, event) ~ group, data = lungcancer,
             distribution = "exact", ties.method = "average-scores")

## Exact Prentice test using average-scores (StatXact 9 manual, p. 222)
logrank_test(Surv(time, event) ~ group, data = lungcancer,
             distribution = "exact", ties.method = "average-scores",
             type = "Prentice")


## Approximative (Monte Carlo) versatile test (Lee, 1996)
rho.gamma <- expand.grid(rho = seq(0, 2, 1), gamma = seq(0, 2, 1))
lee_trafo <- function(y)
    logrank_trafo(y, ties.method = "average-scores",
                  type = "Fleming-Harrington",
                  rho = rho.gamma["rho"], gamma = rho.gamma["gamma"])

it <- independence_test(Surv(time, event) ~ group, data = lungcancer,
                        distribution = approximate(nresample = 10000),
                        ytrafo = function(data)
                            trafo(data, surv_trafo = lee_trafo))
pvalue(it, method = "step-down")
}
\keyword{htest}
\keyword{survival}
