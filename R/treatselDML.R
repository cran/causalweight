#' Binary or multiple  treatment effect evaluation with double machine learning under sample selection/outcome attrition
#' @description Average treatment effect (ATE) estimation for assessing the average effects of discrete (multiple or binary) treatments under sample selection/outcome attrition. Combines estimation based on Neyman-orthogonal score functions with double machine learning to control for confounders in a data-driven way.
#' @param y Dependent variable, may contain missings.
#' @param d Treatment variable, must be discrete, must not contain missings.
#' @param x Covariates, must not contain missings.
#' @param s Selection indicator. Must be 1 if \code{y} is observed (non-missing) and 0 if \code{y} is not observed (missing).
#' @param z Optional instrumental variable(s) for selection \code{s}. If \code{NULL}, outcome selection based on observables (\code{x},\code{d}) - known as "missing at random" - is assumed. If \code{z} is defined, outcome selection based on unobservables - known as "non-ignorable missingness" - is assumed. Default is \code{NULL}.
#' @param selected Must be 1 if ATE is to be estimated for the selected population without missing outcomes. Must be 0 if the ATE is to be estimated for the total population. Default is 0 (ATE for total population). This parameter is ignored if \code{z} is \code{NULL} (under MAR, the ATE in the total population is estimated).
#' @param dtreat Value of the treatment in the treatment group. Default is 1.
#' @param dcontrol Value of the treatment in the control group. Default is 0.
#' @param trim Trimming rule for discarding observations with (products of) propensity scores that are smaller than \code{trim} (to avoid too small denominators in weighting by the inverse of the propensity scores). If \code{selected} is 0 (ATE estimation for the total population), observations with products of the treatment and selection propensity scores that are smaller than \code{trim} are discarded.  If \code{selected} is 1 (ATE estimation for the subpopulation with observed outcomes), observations with treatment propensity scores smaller than \code{trim} are discarded. Default for \code{trim} is 0.01.
#' @param MLmethod Machine learning method for estimating the nuisance parameters based on the \code{SuperLearner} package. Must be either  \code{"lasso"} (default) for lasso estimation,  \code{"randomforest"} for random forests, \code{"xgboost"} for xg boosting,  \code{"svm"} for support vector machines, \code{"ensemble"} for using an ensemble algorithm based on all previously mentioned machine learners, or \code{"parametric"} for linear or logit regression.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @param normalized If set to \code{TRUE}, then the inverse probability-based weights are normalized such that they add up to 1 within treatment groups. Default is \code{TRUE}.
#' @details Estimation of the causal effects of binary or multiple discrete treatments under conditional independence, assuming that confounders jointly affecting the treatment and the outcome can be controlled for by observed covariates, and sample selection/outcome attrition. The latter might either be related to observables, which implies a missing at random assumption, or in addition also to unobservables, if an instrument for sample selection is available. Estimation is based on Neyman-orthogonal score functions for potential outcomes in combination with double machine learning with cross-fitting, see Chernozhukov et al (2018). To this end, one part of the data is used for estimating the model parameters of the treatment and outcome equations based machine learning. The other part of the data is used for predicting the efficient score functions. The roles of the data parts are swapped (using k-fold cross-fitting) and the average treatment effect is estimated based on averaging the predicted efficient score functions in the total sample.  Standard errors are based on asymptotic approximations using the estimated variance of the (estimated) efficient score functions.
#' @return A \code{treatDML} object contains eight components, \code{effect}, \code{se}, \code{pval}, \code{ntrimmed}, \code{meantreat}, \code{meancontrol}, \code{pstreat}, and \code{pscontrol}:
#' @return \code{effect}: estimate of the average treatment effect.
#' @return \code{se}: standard error of the effect.
#' @return \code{pval}: p-value of the effect estimate.
#' @return \code{ntrimmed}: number of discarded (trimmed) observations due to extreme propensity scores.
#' @return \code{meantreat}: Estimate of the mean potential outcome under treatment.
#' @return \code{meancontrol}: Estimate of the mean potential outcome under control.
#' @return \code{pstreat}: P-score estimates for treatment in treatment group.
#' @return \code{pscontrol}: P-score estimates for treatment in control group.
#' @references Bia, M., Huber, M., Laffers, L. (2020): "Double machine learning for sample selection models", working paper, University of Fribourg.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @references van der Laan, M., Polley, E., Hubbard, A. (2007): "Super Learner", Statistical Applications in Genetics and Molecular Biology, 6.
#' @examples # A little example with simulated data (2000 observations)
#' \dontrun{
#' n=2000                            # sample size
#' p=100                             # number of covariates
#' s=2                               # number of covariates that are confounders
#' sigma=matrix(c(1,0.5,0.5,1),2,2)
#' e=(2*rmvnorm(n,rep(0,2),sigma))
#' x=matrix(rnorm(n*p),ncol=p)       # covariate matrix
#' beta=c(rep(0.25,s), rep(0,p-s))   # coefficients determining degree of confounding
#' d=(x%*%beta+rnorm(n)>0)*1         # treatment equation
#' z=rnorm(n)
#' s=(x%*%beta+0.25*d+z+e[,1]>0)*1   # selection equation
#' y=x%*%beta+0.5*d+e[,2]            # outcome equation
#'  y[s==0]=0
#' # The true ATE is equal to 0.5
#' output=treatselDML(y,d,x,s,z)
#' cat("ATE: ",round(c(output$effect),3),", standard error: ",
#'     round(c(output$se),3), ", p-value: ",round(c(output$pval),3))
#' output$ntrimmed}
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm dnorm quantile coef fitted gaussian median
#' @import SuperLearner glmnet ranger xgboost e1071 mvtnorm
#' @export

treatselDML=function(y,d,x,s,z=NULL, selected=0, dtreat=1, dcontrol=0, trim=0.01, MLmethod="lasso", k=3, normalized=TRUE){
  y[s==0]=0
  dtre=1*(d==dtreat); dcon=1*(d==dcontrol)
  scorestreat=hdseltreat(y=y,d=dtre, x=x, s=s, z=z, trim=trim, MLmethod=MLmethod, k=k, selected=selected)
  scorescontrol=hdseltreat(y=y,d=dcon, x=x, s=s, z=z, trim=trim, MLmethod=MLmethod,k=k, selected=selected)
  trimmed=1*(scorescontrol[,7]+scorestreat[,7]>0)        #number of trimmed observations
  scorestreat=scorestreat[trimmed==0,]
  scorescontrol=scorescontrol[trimmed==0,]

  if (selected!=1) {
    if (normalized==FALSE){
      tscores=scorestreat[,2]*scorestreat[,6]*(scorestreat[,3]-scorestreat[,4])/(scorestreat[,5]*scorestreat[,1])+scorestreat[,4]
      cscores=scorescontrol[,2]*scorescontrol[,6]*(scorescontrol[,3]-scorescontrol[,4])/(scorescontrol[,5]*scorescontrol[,1])+scorescontrol[,4]
    }
    if (normalized!=FALSE){
      tscores=(nrow(scorestreat)*scorestreat[,2]*scorestreat[,6]*(scorestreat[,3]-scorestreat[,4])/(scorestreat[,5]*scorestreat[,1]))/(sum(scorestreat[,2]*scorestreat[,6]/(scorestreat[,5]*scorestreat[,1])))+scorestreat[,4]
      cscores=(nrow(scorescontrol)*scorescontrol[,2]*scorescontrol[,6]*(scorescontrol[,3]-scorescontrol[,4])/(scorescontrol[,5]*scorescontrol[,1]))/(sum(scorescontrol[,2]*scorescontrol[,6]/(scorescontrol[,5]*scorescontrol[,1])))+scorescontrol[,4]
    }
  }
  if (selected==1) {
    if (normalized==FALSE){
      tscores=scorestreat[,6]*(scorestreat[,2]*(scorestreat[,3]-scorestreat[,4])/(scorestreat[,5])+scorestreat[,4])/mean(scorestreat[,6])
      cscores=scorescontrol[,6]*(scorescontrol[,2]*(scorescontrol[,3]-scorescontrol[,4])/(scorescontrol[,5])+scorescontrol[,4])/mean(scorescontrol[,6])
    }
    if (normalized!=FALSE){
      tscores=(nrow(scorestreat)*scorestreat[,6]*scorestreat[,2]*(scorestreat[,3]-scorestreat[,4])/scorestreat[,5]) / (sum(scorestreat[,6]*scorestreat[,2]/scorestreat[,5]))   +scorestreat[,6]*scorestreat[,4]/mean(scorestreat[,6])
      cscores=(nrow(scorescontrol)*scorescontrol[,6]*scorescontrol[,2]*(scorescontrol[,3]-scorescontrol[,4])/scorescontrol[,5]) / (sum(scorescontrol[,6]*scorescontrol[,2]/scorescontrol[,5]))  +scorescontrol[,6]*scorescontrol[,4]/mean(scorescontrol[,6])
    }

  }
  meantreat=mean(tscores)
  meancontrol=mean(cscores)
  effect=meantreat - meancontrol
  se=sqrt(mean((tscores-cscores-effect)^2)/sum(tscores))
  pval= 2*pnorm((-1)*abs(effect/se))
  list(effect=effect, se=se, pval=pval, ntrimmed=sum(trimmed), meantreat=meantreat, meancontrol=meancontrol, pstreat=scorestreat[,5], pscontrol=scorescontrol[,5])
}
