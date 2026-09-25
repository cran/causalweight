#' Testing effect homogeneity across studies with double machine learning
#'
#' @param y Dependent variable, must not contain missings.
#' @param d Treatment variable, must be binary (0/1), must not contain missings.
#' @param x Covariates, must not contain missings.
#' @param z Study indicator across which effect homogeneity is tested, must not contain missings. If \code{z} is a factor, a character vector, or a numeric vector taking only integer values, each distinct value defines one study and contrasts are formed for all but the last of them. If \code{z} is numeric and not integer-valued, it is discretised into \code{L} bins of equal size, each of which is then treated as one study.
#' @param MLmethod Machine learning method for estimating the nuisance parameters based on the \code{SuperLearner} package. Must be either \code{"lasso"} (default) for lasso estimation, \code{"randomforest"} for random forests, \code{"xgboost"} for xg boosting, \code{"svm"} for support vector machines, \code{"ensemble"} for using an ensemble algorithm based on all previously mentioned machine learners, or \code{"parametric"} for linear or logit regression.
#' @param k Number of folds in k-fold cross-fitting. Default is 5.
#' @param L Number of bins into which \code{z} is discretised if \code{z} is numeric and not integer-valued. Ignored otherwise, in which case the number of studies is determined by the distinct values of \code{z}. Default is 4.
#' @param trim Trimming rule for discarding observations whose estimated joint probabilities of receiving a treatment state and belonging to a study given \code{x} are smaller than \code{trim} or larger than \code{1-trim} (to avoid too small denominators in weighting by the inverse of these probabilities). Default is 0.05.
#' @param normalized If set to \code{TRUE}, then the inverse probability-based weights are normalized such that they average to one within each combination of treatment state and study. Default is \code{TRUE}.
#' @param seed Default is 123.
#' @details Tests whether the conditional average treatment effect (CATE) of a binary treatment \code{d} on an outcome \code{y} given covariates \code{x} is homogeneous across the studies indexed by \code{z}. The null hypothesis is that for each study, the CATE within that study coincides with the CATE in the remaining studies at any value of \code{x}. If all studies are randomized experiments, a rejection indicates that treatment effects vary with characteristics that differ across studies but are not contained in \code{x}, which speaks to the external validity of the experimental effects. If some studies are experimental while others are observational, a rejection may in addition reflect unobserved confounding in the observational studies, such that differences between the estimands may arise from unobserved confounding, effect heterogeneity, or both. See Armendariz and Huber (2026). Estimation is based on a Neyman-orthogonal score that is quadratic in the doubly robust double difference of conditional mean outcomes across treatment states and studies, in the spirit of Apfel, Hatamyar, Huber and Kueck (2024), combined with \code{k}-fold cross-fitting. The nuisance parameters, namely the conditional mean outcomes \eqn{E[Y|D=d,X,Z=z]} and the joint probabilities \eqn{Pr(D=d,Z=z|X)}, are estimated by the machine learner selected in \code{MLmethod}. The test statistic is asymptotically normal under the null hypothesis.
#' @return A list with the following elements:
#' @return \code{teststat}: estimate of the target parameter, which is zero under the null hypothesis of effect homogeneity across studies.
#' @return \code{se}: standard error of \code{teststat}.
#' @return \code{pval}: p-value of the two-sided test of the null hypothesis, based on the t-statistic \code{teststat/se}.
#' @return \code{n_eff}: effective sample size after trimming and inverse probability weighting.
#' @return \code{n_eff_z}: effective sample size for each study-specific contrast.
#' @return \code{ntrimmed}: number of observations discarded by trimming.
#' @return \code{studies}: labels of the studies for which contrasts are formed.
#' @return \code{psi}: individual contributions to the score function.
#' @return \code{w_full}: inverse probability-based weights, zero for trimmed observations.
#' @references Armendariz, A., Huber, M. (2026): Testing Effect Homogeneity and Confounding in High-Dimensional Experimental and Observational Studies. arXiv:2602.19703.
#' @references Apfel, N., Hatamyar, J., Huber, M., Kueck, J. (2024): Learning control variables and instruments for causal analysis in observational data. arXiv:2407.04448.
#' @examples \dontrun{
#' # Two studies with homogeneous CATEs, such that the test should not reject.
#' set.seed(123)
#' n <- 2000
#' x <- matrix(rnorm(n * 5), ncol = 5)
#' z <- sample(1:2, n, replace = TRUE)                     # two studies
#' d <- rbinom(n, 1, plogis(0.5 * x[, 1] + 0.5 * (z == 2)))
#' y <- d + x[, 2] + rnorm(n)                              # CATE is 1 in both studies
#' output <- catehetDML(y = y, d = d, x = x, z = z)
#' output$teststat
#' output$pval
#'
#' # Same design, but with CATEs differing across studies, such that the test should reject.
#' y <- d + 1 * (z == 2) * d + x[, 2] + rnorm(n)           # CATE is 1 vs. 2
#' output <- catehetDML(y = y, d = d, x = x, z = z)
#' output$teststat
#' output$pval
#' }
#' @importFrom stats binomial gaussian model.matrix pnorm predict quantile sd
#' @import glmnet ranger xgboost nnls
#' @export

catehetDML <- function(y, d, x, z,
                       MLmethod   = "lasso",
                       k          = 5,
                       L          = 4,
                       trim       = 0.05,
                       normalized = TRUE,
                       seed       = 123) {

  #### checking the arguments ####
  if (!MLmethod %in% c("lasso", "randomforest", "xgboost", "svm", "ensemble", "parametric"))
    stop("'MLmethod' must be one of 'lasso', 'randomforest', 'xgboost', 'svm', 'ensemble', 'parametric'.")
  if (!(is.numeric(k) && length(k) == 1 && k >= 2))
    stop("'k' must be a single number of at least 2.")
  if (!(is.numeric(trim) && length(trim) == 1 && trim >= 0 && trim < 0.5))
    stop("'trim' must be a single number in [0, 0.5).")

  y <- as.numeric(y)
  d <- as.numeric(d)
  if (!all(d %in% c(0, 1)))   stop("'d' must be binary (0/1).")
  if (length(y) != length(d)) stop("'y' and 'd' must have the same length.")
  if (length(z) != length(d)) stop("'z' must have the same length as 'd'.")
  if (NROW(x) != length(d))   stop("'x' must have the same number of rows as 'd'.")
  if (anyNA(y) || anyNA(d) || anyNA(x) || anyNA(z))
    stop("Inputs must not contain missing values.")

  set.seed(seed)

  n   <- length(y)
  Xdf <- data.frame(x)

  #### encoding the study indicator ####
  Zb <- catehet_zdummies(z, L = L)
  L  <- ncol(Zb)

  #### estimating the nuisance parameters by k-fold cross-fitting ####
  mu11_mat <- mu10_mat <- mu01_mat <- mu00_mat <- matrix(NA, n, L)
  p_1_z    <- p_0_z    <- p_1_zm   <- p_0_zm   <- matrix(NA, n, L)

  for (l in seq_len(L)) {

    mu11_mat[, l] <- catehet_MLmean(y, data.frame(Xdf, d = rep(1, n), z_ind = rep(1, n)),
                                    MLmethod = MLmethod, k = k)
    mu10_mat[, l] <- catehet_MLmean(y, data.frame(Xdf, d = rep(0, n), z_ind = rep(1, n)),
                                    MLmethod = MLmethod, k = k)
    mu01_mat[, l] <- catehet_MLmean(y, data.frame(Xdf, d = rep(1, n), z_ind = rep(0, n)),
                                    MLmethod = MLmethod, k = k)
    mu00_mat[, l] <- catehet_MLmean(y, data.frame(Xdf, d = rep(0, n), z_ind = rep(0, n)),
                                    MLmethod = MLmethod, k = k)

    p_1_z[, l]  <- catehet_MLmean(as.integer(d == 1 & Zb[, l] == 1), Xdf,
                                  MLmethod = MLmethod, k = k)
    p_0_z[, l]  <- catehet_MLmean(as.integer(d == 0 & Zb[, l] == 1), Xdf,
                                  MLmethod = MLmethod, k = k)
    p_1_zm[, l] <- catehet_MLmean(as.integer(d == 1 & Zb[, l] == 0), Xdf,
                                  MLmethod = MLmethod, k = k)
    p_0_zm[, l] <- catehet_MLmean(as.integer(d == 0 & Zb[, l] == 0), Xdf,
                                  MLmethod = MLmethod, k = k)
  }

  #### building the score ####
  tmp <- catehet_psi(y, mu11_mat, mu10_mat, mu01_mat, mu00_mat,
                     p_1_z, p_0_z, p_1_zm, p_0_zm,
                     d, Zb, trim = trim, normalized = normalized)

  psi     <- tmp$psi
  w_full  <- tmp$w_full
  n_eff_z <- tmp$n_eff_z

  #### assembling the test ####
  keep  <- which(w_full > 0)
  w_k   <- w_full[keep]
  n_eff <- if (length(w_k)) sum(w_k)^2 / sum(w_k^2) else 0

  theta_hat <- if (length(keep))     mean(psi[keep]) else NA
  sigma_hat <- if (length(keep) > 1) sd(psi[keep])   else NA
  se        <- if (!is.na(sigma_hat) && n_eff > 0) sigma_hat / sqrt(n_eff) else NA
  p_value   <- if (!is.na(se) && !is.na(theta_hat)) 2 * pnorm(-abs(theta_hat / se)) else NA

  list(teststat = theta_hat,
       se       = se,
       pval     = p_value,
       n_eff    = n_eff,
       n_eff_z  = n_eff_z,
       ntrimmed = n - length(keep),
       studies  = colnames(Zb),
       psi      = psi,
       w_full   = w_full)
}
