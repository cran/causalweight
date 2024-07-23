#' Difference-in-differences based on inverse probability weighting
#' @description Difference-in-differences-based estimation of the average treatment effect on the treated in the post-treatment period, given a binary treatment with one pre- and one post-treatment period. Permits controlling for differences in observed covariates across treatment groups and/or time periods based on inverse probability weighting.
#' @param y Dependent variable, must not contain missings.
#' @param d Treatment, must be binary (either 1 or 0), must not contain missings.
#' @param t Time period, must be binary, 0 for pre-treatment and 1 for post-treatment, must not contain missings.
#' @param x Covariates to be controlled for by inverse probability weighting. Default is \code{NULL}.
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @param trim Trimming rule for discarding observations with extreme propensity scores in the 3 reweighting steps, which reweight (1) treated in the pre-treatment period, (2) non-treated in the post-treatment period, and (3) non-treated in the pre-treatment period according to the covariate distribution of the treated in the post-treatment period. Default is 0.05, implying that observations with a probability lower than 5 percent of not being treated in some weighting step are discarded.
#' @param cluster A cluster ID for block or cluster bootstrapping when units are clustered rather than iid. Must be numerical. Default is \code{NULL} (standard bootstrap without clustering).
#' @details Estimation of the average treatment effect on the treated in the post-treatment period based Difference-in-differences. Inverse probability weighting is used to control for differences in covariates across treatment groups and/or over time. That is, (1) treated observations in the pre-treatment period, (2) non-treated observations in the post-treatment period, and (3) non-treated observations in the pre-treatment period are reweighted according to the covariate distribution of the treated observations in the post-treatment period. The respective propensity scores are obtained by probit regressions.
#' @return A didweight object contains 4 components, \code{eff}, \code{se}, \code{pvalue}, and \code{ntrimmed}.
#' @return \code{eff}: estimate of the average treatment effect on the treated in the post-treatment period.
#' @return \code{se}: standard error obtained by bootstrapping the effect.
#' @return \code{pvalue}: p-value based on the t-statistic.
#' @return \code{ntrimmed}: total number of discarded (trimmed) observations in any of the 3 reweighting steps due to extreme propensity score values.
#' @references Abadie, A. (2005): "Semiparametric Difference-in-Differences Estimators", The Review of Economic Studies, 72, 1-19.
#' @references Lechner, M. (2011): "The Estimation of Causal Effects by Difference-in-Difference Methods", Foundations and Trends in Econometrics, 4, 165-224.
#' @examples # A little example with simulated data (4000 observations)
#' \dontrun{
#' n=4000                            # sample size
#' t=1*(rnorm(n)>0)                  # time period
#' u=rnorm(n)                        # time constant unobservable
#' x=0.5*t+rnorm(n)                  # time varying covariate
#' d=1*(x+u+rnorm(n)>0)              # treatment
#' y=d*t+t+x+u+rnorm(n)              # outcome
#' # The true effect equals 1
#' didweight(y=y,d=d,t=t,x=x, boot=199)}
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm dnorm quantile
#' @export

didweight<-function(y,d,t,x=NULL, boot=1999, trim=0.05, cluster=NULL){
  results=ipw.did(y=y,d=d,t=t, x=x, trim=trim)
  beff=bootstrap.did(y=y,d=d,t=t, x=x, boot=boot,trim=trim, cluster=cluster)
  se=sd(beff)
  list (effect=results[1], se=se, pvalue=2*pnorm(-abs(results[1]/se)), ntrimmed=results[2])
}
