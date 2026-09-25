#' Detection of instruments and control variables and validity testing with double machine learning
#'
#' @description Tests the validity of a pre-specified instrument and pre-specified control variables (mode 1) or
#'   learns the partition of instruments and control variables from the data
#'   (mode 2), following the two-step procedure of Apfel, Hatamyar, Huber, and
#'   Kueck (2025).
#'
#'   \strong{Mode 1 — user-specified \code{x} and \code{z}:} The user provides
#'   the covariates \code{x} and a single candidate instrument \code{z}
#'   directly.  The function skips instrument selection and runs only the
#'   validity test: it tests whether \code{z} is conditionally
#'   independent of \code{y} given \code{d} and \code{x}.
#'
#'   \strong{Mode 2 — data-driven selection via \code{q}:} The user provides a
#'   matrix \code{q} of candidate variables without specifying which are
#'   instruments and which are controls.  The function runs the full two-stage
#'   procedure: (1) a first-stage DML relevance test selects candidate instruments
#'   associated with \code{d}, and (2) the validity test retains
#'   those that are conditionally independent of \code{y}.
#'
#' @param y Outcome variable, numeric vector.  Must not contain missings.
#' @param d Treatment variable, numeric vector.  Must not contain missings.
#' @param z \strong{Mode 1 only.}  Candidate instrument, numeric vector.  Must
#'   not contain missings.  Provide together with \code{x}; leave \code{NULL}
#'   (default) to use mode 2 with \code{q}.
#' @param x \strong{Mode 1 only.}  Covariate matrix.  Must not contain
#'   missings.  Provide together with \code{z}; leave \code{NULL} (default)
#'   to use mode 2 with \code{q}.
#' @param q \strong{Mode 2 only.}  Matrix of candidate variables whose columns
#'   are iteratively considered as potential instruments (with the remaining columns
#'   serving as control variables).  Must not contain missings and must have at least
#'   two columns.  Leave \code{NULL} (default) to use mode 1 with \code{x}
#'   and \code{z}.
#' @param alpha \strong{Mode 2 only.}  Significance level for the first-stage
#'   relevance test.  A candidate passes if its DML p-value is below
#'   \code{alpha}.  The paper recommends \code{0.1 / log(n)}.  Default is
#'   \code{0.1 / log(length(y))}.
#' @param critval Significance level for the exclusion-restriction test.
#'   Candidates whose p-value \emph{exceeds} \code{critval} are classified as
#'   valid IVs (null of conditional independence not rejected).  Note that a higher value is a more conservative choice; Apfel et al. recommend 0.3. Default is \code{0.3}.
#' @param MLmethod Machine learning method for estimating nuisance parameters
#'   via the \code{SuperLearner} package.  Must be one of \code{"lasso"}
#'   (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"},
#'   \code{"ensemble"}, or \code{"parametric"}.
#' @param k Number of folds in k-fold cross-fitting.  Default is \code{3}.
#' @param trim Trimming threshold for propensity scores: observations whose
#'   estimated propensity score falls below \code{trim} or above
#'   \code{1 - trim} are discarded.  Default is \code{0.01}.
#' @param L Number of partition bins for non-binary instruments with 15 or
#'   more unique values.  Variables with fewer than 15 unique values use their
#'   natural categories.  Default is \code{4}.
#' @param seed Random seed.  Default is \code{123}.
#'
#' @details
#'   \strong{Mode 1} runs the validity test of Apfel et al.
#'   (2025) directly on the user-supplied \code{z} and \code{x}, testing
#'   \eqn{H_0: E[Y \mid D, X] = E[Y \mid D, X, Z]}.  This is appropriate
#'   when the researcher has a single pre-specified candidate IV and wants to
#'   assess its validity (and that of the control variables) without selection.
#'
#'   \strong{Mode 2} implements the full data-driven algorithm of Apfel et al.
#'   (2025).  For a matrix \eqn{Q} of \eqn{p} candidate variables it
#'   constructs partitions \eqn{P_j = \{Z = Q_j,\, X = Q_{[j]}\}} and:
#'   \enumerate{
#'     \item Fits a DML partially linear regression (PLR) of \code{d} on each
#'       column \eqn{Q_j} controlling for \eqn{Q_{[j]}}, retaining columns
#'       whose p-value is below \code{alpha} as first-stage candidates
#'       \eqn{\hat{S}}.
#'     \item Applies the validity test to each candidate in
#'       \eqn{\hat{S}}, classifying those with p-value above \code{critval} as
#'       valid IVs \eqn{\hat{V}}.
#'   }
#'   The detected set is \eqn{\hat{P}_{\text{pass}} = \hat{S} \cap \hat{V}}.
#'
#'   In both modes the validity test uses the doubly robust score
#'   of Apfel et al. (2025): eq. (7) for binary instruments and the
#'   multi-partition score of eq. (18) / Appendix C for non-binary instruments.
#'
#' @return
#'   \strong{Mode 1} returns a list with:
#'   \describe{
#'     \item{\code{teststat}}{Test statistic (should be near zero under the null).}
#'     \item{\code{se}}{Standard error of the test statistic.}
#'     \item{\code{pval}}{P-value of the exclusion-restriction test.}
#'     \item{\code{n_eff}}{Number of observations retained after trimming.}
#'     \item{\code{critval}}{Significance level used for the exclusion test.}
#'   }
#'
#'   \strong{Mode 2} returns a list with:
#'   \describe{
#'     \item{\code{valid_IVs}}{Integer vector of column indices of \code{q}
#'       detected as valid IVs.  Empty if none found.}
#'     \item{\code{firststage_pass}}{Column indices passing the first-stage
#'       test (\eqn{\hat{S}}).}
#'     \item{\code{firststage_pvals}}{Named numeric vector of first-stage
#'       p-values for all columns of \code{q}.}
#'     \item{\code{exclusion_pvals}}{Named numeric vector of
#'       exclusion-restriction p-values.  \code{NA} for columns not in
#'       \eqn{\hat{S}}.}
#'     \item{\code{exclusion_results}}{Named list of full test output objects
#'       for each first-stage candidate.}
#'     \item{\code{alpha}}{Significance level used for the first-stage test.}
#'     \item{\code{critval}}{Significance level used for the exclusion test.}
#'   }
#'
#' @references Apfel, N., Hatamyar, J., Huber, M., Kueck, J. (2025): "Learning control variables and instruments for causal analysis in observational data," arXiv:2407.04448.
#' @references Huber, M., Kloiber, K.,  Lafférs, L. (2026): "Testing Full Mediation of Treatment Effects and the Identifiability of Causal Mechanisms,"  arXiv:2603.04109.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 2000; p <- 10
#' Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
#' Q_raw <- mvtnorm::rmvnorm(n, rep(0, p), Sigma)
#' pis   <- 1 / (1 + exp(-2 * Q_raw))
#' Q     <- matrix(rbinom(n * p, 1, pis), nrow = n)
#' beta  <- c(0.8 / (1:4), rep(0, 5))
#' W <- rnorm(n); U <- rnorm(n); V <- rnorm(n)
#' D <- as.numeric(Q[, 1:9] %*% beta + Q[, p] + V > 0)
#' Y <- D + Q[, 1:9] %*% beta + W + U   # Q[, p] has no direct effect on Y
#'
#' # Mode 1: test a pre-specified instrument (last column) directly
#' res1 <- detectIV(y = Y, d = D, z = Q[, p], x = Q[, 1:9])
#' cat("Mode 1 p-value:", round(res1$pval, 3), "\n")
#'
#' # Mode 2: learn instruments and controls from Q
#' res2 <- detectIV(y = Y, d = D, q = Q, alpha = 0.1 / log(n))
#' cat("Mode 2 column indices of detected valid IVs:", res2$valid_IVs, "\n")
#' }
#'
#' @importFrom stats pnorm rnorm printCoefmat setNames
#' @import checkmate
#' @export

detectIV <- function(y, d,
                     z = NULL, x = NULL,   # Mode 1
                     q = NULL,             # Mode 2
                     alpha   = 0.1 / log(length(y)),
                     critval = 0.3,
                     MLmethod = "lasso",
                     k    = 3,
                     trim = 0.01,
                     L    = 4,
                     seed = 123) {



  ## ---- Determine mode and validate arguments --------------------------------

  mode1 <- !is.null(z) && !is.null(x)
  mode2 <- !is.null(q)

  if (mode1 && mode2)
    stop("Provide either (z, x) for mode 1 OR q for mode 2, not both.")
  if (!mode1 && !mode2)
    stop("Provide either (z, x) for mode 1 or q for mode 2.")
  if (!is.null(z) && is.null(x))
    stop("'x' must be provided alongside 'z' in mode 1.")
  if (is.null(z) && !is.null(x))
    stop("'z' must be provided alongside 'x' in mode 1.")

  checkmate::assertNumeric(y, any.missing = FALSE)
  checkmate::assertNumeric(d, any.missing = FALSE)
  checkmate::assertNumber(critval, lower = 0, upper = 1)
  checkmate::assertChoice(MLmethod,
                          c("lasso", "randomforest", "xgboost", "svm", "ensemble", "parametric"))
  checkmate::assertIntegerish(k, lower = 2)
  checkmate::assertNumber(trim, lower = 0, upper = 0.5)
  checkmate::assertCount(L, positive = TRUE)
  checkmate::assertNumber(seed)

  set.seed(seed)

  n <- length(y)
  if(3*k > n) stop("'k' cannot exceed a third of the sample size.")
  if (length(d) != n) stop("'y' and 'd' must have the same length.")

  ## ==========================================================================
  ## MODE 1: user supplies z and x — run exclusion-restriction test directly
  ## ==========================================================================
  if (mode1) {

    checkmate::assertNumeric(z, any.missing = FALSE)
    checkmate::assertMatrix(as.matrix(x), any.missing = FALSE)
    if (length(z) != n) stop("'z' must have the same length as 'y'.")
    if (nrow(as.matrix(x)) != n) stop("'x' must have the same number of rows as 'y'.")

    res <- test_conditional_independence(
      Y        = y,
      Z        = z,
      D        = d,
      X        = as.matrix(x),
      MLmethod = MLmethod,
      K        = k,
      epsilon  = trim,
      L        = L
    )

    return(c(res, list(critval = critval)))
  }

  ## ==========================================================================
  ## MODE 2: user supplies q — full two-stage data-driven procedure
  ## ==========================================================================

  checkmate::assertMatrix(q, any.missing = FALSE, min.cols = 2)
  checkmate::assertNumber(alpha, lower = 0, upper = 1)
  if (nrow(q) != n) stop("'q' must have the same number of rows as 'y'.")

  if (is.null(colnames(q))) colnames(q) <- paste0("Q", seq_len(ncol(q)))
  qnames <- colnames(q)
  p <- ncol(q)

  ## ---- STEP 1: First-stage relevance (DML-PLR) ----------------------------
  ##
  ## For each column j, fit PLR: D ~ Q_j | Q_{-j}.
  ## PLR estimate: theta_hat = sum(e * V) / sum(e^2), SE via sandwich.

  firststage_pvals <- setNames(rep(NA_real_, p), qnames)

  for (j in seq_len(p)) {

    z_j  <- q[, j]
    x_mj <- q[, -j, drop = FALSE]
    colnames(x_mj) <- paste0("V", seq_len(ncol(x_mj)))

    d_hat <- MLmean2(y = d,   x = x_mj, MLmethod = MLmethod, k = k, seed = seed)
    z_hat <- MLmean2(y = z_j, x = x_mj, MLmethod = MLmethod, k = k, seed = seed)

    V_res <- d   - d_hat
    e_res <- z_j - z_hat

    denom <- sum(e_res^2)
    if (denom < .Machine$double.eps) next

    theta_hat <- sum(e_res * V_res) / denom
    score_i   <- e_res * (V_res - theta_hat * e_res)
    se_hat    <- sqrt(sum(score_i^2) / denom^2)

    firststage_pvals[j] <- 2 * pnorm(-abs(theta_hat / se_hat))
  }

  firststage_pass <- which(firststage_pvals < alpha & !is.na(firststage_pvals))
  names(firststage_pass) <- qnames[firststage_pass]

  ## ---- STEP 2: Exclusion-restriction test ----------------------------------

  exclusion_pvals   <- setNames(rep(NA_real_, p), qnames)
  exclusion_results <- setNames(vector("list", p), qnames)

  if (length(firststage_pass) == 0) {
    warning("No candidate variable passed the first-stage relevance test. ",
            "Returning empty valid_IVs.")
  } else {

    candidates_df <- as.data.frame(q[, firststage_pass, drop = FALSE])
    non_cand_cols <- setdiff(seq_len(p), firststage_pass)
    controls_base <- if (length(non_cand_cols) > 0)
      as.data.frame(q[, non_cand_cols, drop = FALSE])
    else NULL

    for (i in seq_along(firststage_pass)) {

      j   <- firststage_pass[i]
      z_j <- q[, j]

      other_cands <- candidates_df[, -i, drop = FALSE]
      controls_i  <- if (!is.null(controls_base))
        cbind(other_cands, controls_base)
      else other_cands

      result_i <- tryCatch(
        test_conditional_independence(
          Y        = y,
          Z        = z_j,
          D        = d,
          X        = as.matrix(controls_i),
          MLmethod = MLmethod,
          K        = k,
          epsilon  = trim,
          L        = L
        ),
        error = function(e) {
          warning("Exclusion-restriction test failed for column '", qnames[j],
                  "': ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(result_i)) {
        exclusion_results[[qnames[j]]] <- result_i
        pval_i <- result_i$pval
        if (is.na(pval_i) && isTRUE(result_i$teststat == 0)) pval_i <- 1
        exclusion_pvals[j] <- pval_i
      }
    }
  }

  excl_pass <- which(exclusion_pvals > critval & !is.na(exclusion_pvals))
  valid_IVs <- intersect(firststage_pass, excl_pass)

  list(
    valid_IVs         = valid_IVs,
    firststage_pass   = firststage_pass,
    firststage_pvals  = firststage_pvals,
    exclusion_pvals   = exclusion_pvals,
    exclusion_results = exclusion_results[qnames[firststage_pass]],
    alpha             = alpha,
    critval           = critval
  )
}
