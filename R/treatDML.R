#' Binary or multiple  discrete treatment effect evaluation with double machine learning
#' @description Treatment effect estimation for assessing the average effects of discrete (multiple or binary) treatments. Combines estimation based on (doubly robust) efficient score functions with double machine learning to control for confounders in a data-driven way.
#' @param y Dependent variable, must not contain missings.
#' @param d Treatment variable, must be discrete, must not contain missings.
#' @param x Covariates, must not contain missings.
#' @param s  Indicator function for defining a subpopulation for whom the treatment effect is estimated as a function of the subpopulation's distribution of \code{x}. Default is \code{NULL} (estimation of the average treatment effect in the total population).
#' @param dtreat Value of the treatment in the treatment group. Default is 1.
#' @param dcontrol Value of the treatment in the control group. Default is 0.
#' @param trim Trimming rule for discarding observations with treatment propensity scores that are smaller than \code{trim} or larger than \code{1-trim} (to avoid too small denominators in weighting by the inverse of the propensity scores). Default is 0.01.
#' @param MLmethod Machine learning method for estimating the nuisance parameters based on the \code{SuperLearner} package. Must be either  \code{"lasso"} (default) for lasso estimation,  \code{"randomforest"} for random forests, \code{"xgboost"} for xg boosting,  \code{"svm"} for support vector machines, \code{"ensemble"} for using an ensemble algorithm based on all previously mentioned machine learners, or \code{"parametric"} for linear or logit regression.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @param normalized If set to \code{TRUE}, then the inverse probability-based weights are normalized such that they add up to 1 within treatment groups. Default is \code{TRUE}.
#' @details Estimation of the causal effects of binary or multiple discrete treatments under conditional independence, assuming that confounders jointly affecting the treatment and the outcome can be controlled for by observed covariates. Estimation is based on the (doubly robust) efficient score functions for potential outcomes in combination with double machine learning with cross-fitting, see Chernozhukov et al (2018). To this end, one part of the data is used for estimating the model parameters of the treatment and outcome equations based machine learning. The other part of the data is used for predicting the efficient score functions. The roles of the data parts are swapped (using k-fold cross-fitting) and the average treatment effect is estimated based on averaging the predicted efficient score functions in the total sample.
#' Standard errors are based on asymptotic approximations using the estimated variance of the (estimated) efficient score functions.
#' @return A \code{treatDML} object contains eight components, \code{effect}, \code{se}, \code{pval}, \code{ntrimmed}, \code{meantreat}, \code{meancontrol}, \code{pstreat}, and \code{pscontrol}:
#' @return \code{effect}: estimate of the average treatment effect.
#' @return \code{se}: standard error of the effect.
#' @return \code{pval}: p-value of the effect estimate.
#' @return \code{ntrimmed}: number of discarded (trimmed) observations due to extreme propensity scores.
#' @return \code{meantreat}: Estimate of the mean potential outcome under treatment.
#' @return \code{meancontrol}: Estimate of the mean potential outcome under control.
#' @return \code{pstreat}: P-score estimates for treatment in treatment group.
#' @return \code{pscontrol}: P-score estimates for treatment in control group.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @references van der Laan, M., Polley, E., Hubbard, A. (2007): "Super Learner", Statistical Applications in Genetics and Molecular Biology, 6.
#' @examples # A little example with simulated data (2000 observations)
#' \dontrun{
#' n=2000                            # sample size
#' p=100                             # number of covariates
#' s=2                               # number of covariates that are confounders
#' x=matrix(rnorm(n*p),ncol=p)       # covariate matrix
#' beta=c(rep(0.25,s), rep(0,p-s))   # coefficients determining degree of confounding
#' d=(x%*%beta+rnorm(n)>0)*1         # treatment equation
#' y=x%*%beta+0.5*d+rnorm(n)       # outcome equation
#' # The true ATE is equal to 0.5
#' output=treatDML(y,d,x)
#' cat("ATE: ",round(c(output$effect),3),", standard error: ",
#'     round(c(output$se),3), ", p-value: ",round(c(output$pval),3))
#' output$ntrimmed}

#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm dnorm quantile coef fitted gaussian median
#' @import SuperLearner glmnet ranger xgboost e1071 mvtnorm
#' @export

treatDML=function(y,d,x, s=NULL, dtreat=1, dcontrol=0, trim=0.01, MLmethod="lasso", k=3, normalized=TRUE){
  dtre=1*(d==dtreat); dcon=1*(d==dcontrol)
  scorestreat=hdtreat(y=y,d=dtre, x=x, s=s, trim=trim, MLmethod=MLmethod, k=k)
  scorescontrol=hdtreat(y=y,d=dcon, x=x, s=s, trim=trim, MLmethod=MLmethod,k=k)
  trimmed=1*(scorescontrol[,7]+scorestreat[,7]>0)        #number of trimmed observations
  scorestreat=scorestreat[trimmed==0,]
  scorescontrol=scorescontrol[trimmed==0,]
  if (normalized==FALSE){
    tscores=(scorestreat[,1]*scorestreat[,2]*(scorestreat[,3]-scorestreat[,4])/(scorestreat[,5])+scorestreat[,6]*scorestreat[,4])/mean(scorestreat[,6])
    cscores=(scorescontrol[,1]*scorescontrol[,2]*(scorescontrol[,3]-scorescontrol[,4])/(scorescontrol[,5])+scorescontrol[,6]*scorescontrol[,4])/mean(scorescontrol[,6])
  }
  if (normalized!=FALSE){
    ntreat=nrow(scorestreat)
    weightsumtreat=sum(scorestreat[,1]*scorestreat[,2]/(scorestreat[,5]))
    tscores=(ntreat*scorestreat[,1]*scorestreat[,2]*(scorestreat[,3]-scorestreat[,4])/(scorestreat[,5]))/weightsumtreat+(scorestreat[,6]*scorestreat[,4])/mean(scorestreat[,6])
    ncontrol=nrow(scorescontrol)
    weightsumcontrol=sum(scorescontrol[,1]*scorescontrol[,2]/(scorescontrol[,5]))
    cscores=(ncontrol*scorescontrol[,1]*scorescontrol[,2]*(scorescontrol[,3]-scorescontrol[,4])/(scorescontrol[,5]))/weightsumcontrol+(scorescontrol[,6]*scorescontrol[,4])/mean(scorescontrol[,6])
  }
    meantreat=mean(tscores)
    meancontrol=mean(cscores)
    effect=meantreat - meancontrol
    se=sqrt(mean((tscores-cscores-effect)^2)/length(tscores))
    pval= 2*pnorm((-1)*abs(effect/se))
  list(effect=effect, se=se, pval=pval, ntrimmed=sum(trimmed), meantreat=meantreat, meancontrol=meancontrol, pstreat=scorestreat[,5], pscontrol=scorescontrol[,5])
}

