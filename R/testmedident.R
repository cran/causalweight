#' Test for identification in causal mediation and dynamic treatment models
#'
#' This function tests for identification in causal mediation and dynamic treatment models based on covariates and instrumental variables using machine learning methods.
#'
#' @param y Outcome variable.
#' @param d Treatment variable.
#' @param m Mediator variable (optional).
#' @param x Baseline covariates (prior to treatment assignment).
#' @param w Post-treatment covariates (prior to mediator assignment, optional).
#' @param z1 Instrument for the treatment.
#' @param z2 Instrument for the mediator (optional).
#' @param testmediator Logical indicating if the mediator should be used as dependent variable (in addition to outcome \code{y}) when testing if the effect of treatment \code{d} is identified. Default is TRUE.
#' @param seed Random seed for sample splitting in cross-fitting. Default is 123.
#' @param MLmethod Machine learning method for estimating conditional outcome/mediator means required for testing. Default is "lasso".
#' @param k Number of cross-fitting folds. Default is 3.
#' @param zeta_sigma Tuning parameter defining the standard deviation of a random, mean zero, and normal variable that is added to the test statistic to avoid a degenerate distribution of test statistic under the null hypothesis. \code{zeta_sigma} gauges the trade-off between power and size of the test. Default is the minimum of 0.5 and 500/(# of observations).
#' @details This function implements a hypothesis test for identifying causal effects in mediation and dynamic treatment models involving sequential assignment of treatment and mediator variables. The test jointly verifies the exogeneity/ignorability of treatment and mediator variables conditional on covariates and the validity of (distinct) instruments for the treatment and mediator (ignorability of instrument assignment and exclusion restriction). If the null hypothesis holds, dynamic and pathwise causal effects may be identified based on the sequential exogeneity/ignorability of the treatment and the mediator given the covariates. The function employs machine learning techniques to control for covariates in a data-driven manner.  
#' @return A list with the following components:
#' \item{teststat}{Test statistic.}
#' \item{se}{Standard error of the test statistic.}
#' \item{pval}{Two-sided p-value of the test.}
#' @references Huber, M., Kloiber, K., and Lafférs, L. (2024): "Testing identification in mediation and dynamic treatment models", arXiv preprint 2406.13826.
#' @examples
#' \dontrun{
#' # Example with simulated data in which null hypothesis holds
#' n=2000
#' x=rnorm(n)
#' z1=rnorm(n)
#' z2=rnorm(n)
#' d=1*(0.5*x+0.5*z1+rnorm(n)>0)      # Treatment equation
#' m=0.5*x+0.5*d+0.5*z2+rnorm(n)      # Mediator equation
#' y=0.5*x+d+0.5*m+rnorm(n)           # Outcome equation
#' # Run test and report p-value
#' testmedident(y=y, d=d, m=m, x=x, z1=z1, z2=z2)$pval
#' }
#' @importFrom SuperLearner SuperLearner
#' @import sandwich
#' @export
testmedident=function(y, d, m=NULL, x, w=NULL, z1, z2=NULL, testmediator=TRUE, seed = 123, MLmethod ="lasso", k = 3, zeta_sigma = min(0.5,500/length(y))){
  zeta=rnorm(length(y),0,sd = zeta_sigma) # sample from a normal distribution to avoid degenerate distribution of test statistic under the null
  scoresinstr=MLmean(y = y, x = cbind(x,d,z1), d=d, MLmethod = MLmethod, k = k, zeta = zeta, seed = seed)
  scoresnoinstr=MLmean(y = y, x = cbind(x,d), d=d, MLmethod = MLmethod, k = k ,zeta = zeta, seed = seed)
  if (testmediator==TRUE & !is.null(m)){   # tests with mediator as outcome
    scoresinstr=rbind(scoresinstr, MLmean(y = m, x = cbind(x,d,z1), d=d, MLmethod = MLmethod, k = k, zeta = zeta, seed = seed))
    scoresnoinstr=rbind(scoresnoinstr, MLmean(y = m, x = cbind(x,d), d=d, MLmethod = MLmethod, k = k ,zeta = zeta, seed = seed))
  }
  if (!is.null(m) & !is.null(z2)){   # tests with mediator as treatment
    if (!is.null(w)) x=cbind(x,w)
    scoresinstr=rbind(scoresinstr, MLmean(y = y, x = cbind(x,d,m,z2), d=d, MLmethod = MLmethod, k = k, zeta = zeta, seed = seed))
    scoresnoinstr=rbind(scoresnoinstr, MLmean(y = y, x = cbind(x,d,m), d=d, MLmethod = MLmethod, k = k ,zeta = zeta, seed = seed))
  }
  outinstr=scoresinstr[,1]
  outnoinstr=scoresnoinstr[,1]
  cluster=c()
  if (testmediator==TRUE & (is.null(m) | is.null(z2))) cluster=rep(c(1:length(y)),2)   # clustering because same observations are used twice when also using mediator as outcome
  if (testmediator==TRUE & !is.null(m) & !is.null(z2)) cluster=rep(c(1:length(y)),3)   # same observations used three times when using mediator as outcome and treatment
  if (testmediator==FALSE & !is.null(m) & !is.null(z2)) cluster=rep(c(1:length(y)),2)  # some observations used twice when also using mediator as treatment
  stat=(outinstr-outnoinstr)^2+scoresinstr[,2]; teststat=mean(stat)
  if (is.null(cluster)) {
    se=summary(lm(stat~1))$coefficients[, "Std. Error"]}  else {  # se not accounting for clustering
      se=sqrt(vcovCL(lm(stat~1), cluster = cluster))[1,1]}        # se accounting for clustering
  pval=2*pnorm((-1)*abs(teststat/se))
  list(teststat=teststat, se=se, pval=pval)
}   