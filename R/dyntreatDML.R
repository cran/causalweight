#' Dynamic treatment effect evaluation with double machine learning
#' @description Dynamic treatment effect estimation for assessing the average effects of sequences of treatments (consisting of two sequential treatments). Combines estimation based on (doubly robust) efficient score functions with double machine learning to control for confounders in a data-driven way.
#' @param y2 Dependent variable in the second period (=outcome period), must not contain missings.
#' @param d1 Treatment in the first period, must be discrete, must not contain missings.
#' @param d2 Treatment in the second period, must be discrete, must not contain missings.
#' @param x0 Covariates in the baseline period (prior to the treatment in the first period), must not contain missings.
#' @param x1 Covariates in the first period (prior to the treatment in the second period), must not contain missings.
#' @param s  Indicator function for defining a subpopulation for whom the treatment effect is estimated as a function of the subpopulation's distribution of \code{x0}. Default is \code{NULL} (estimation of the treatment effect in the total population).
#' @param d1treat Value of the first treatment in the treatment sequence. Default is 1.
#' @param d2treat Value of the second treatment in the treatment sequence. Default is 1.
#' @param d1control Value of the first treatment in the control sequence. Default is 0.
#' @param d2control Value of the second treatment in the control sequence. Default is 0.
#' @param trim Trimming rule for discarding observations with products of treatment propensity scores in the first and second period that are smaller than \code{trim} (to avoid too small denominators in weighting by the inverse of the propensity scores). Default is 0.01.
#' @param MLmethod Machine learning method for estimating the nuisance parameters based on the \code{SuperLearner} package. Must be either  \code{"lasso"} (default) for lasso estimation,  \code{"randomforest"} for random forests, \code{"xgboost"} for xg boosting,  \code{"svm"} for support vector machines, \code{"ensemble"} for using an ensemble algorithm based on all previously mentioned machine learners, or \code{"parametric"} for linear or logit regression.
#' @param fewsplits If set to \code{TRUE}, the same training data are used for estimating a nested model of conditional mean outcomes, namely \code{E[E[y2|d1,d2,x0,x1]|d1,x0]}. If \code{fewsplits} is \code{FALSE}, the training data are split for the sequential estimation of the nested model. Default of \code{fewsplits} is \code{FALSE}.
#' @param normalized If set to \code{TRUE}, then the inverse probability-based weights are normalized such that they add up to 1 within treatment groups. Default is \code{TRUE}.
#' @details Estimation of the causal effects of sequences of two treatments under sequential conditional independence, assuming that all confounders of the treatment in either period and the outcome of interest are observed. Estimation is based on the (doubly robust) efficient score functions for potential outcomes, see e.g. Bodory, Huber, and Laffers (2020), in combination with double machine learning with cross-fitting, see Chernozhukov et al (2018). To this end, one part of the data is used for estimating the model parameters of the treatment and outcome equations based machine learning. The other part of the data is used for predicting the efficient score functions. The roles of the data parts are swapped (using 3-fold cross-fitting) and the average dynamic treatment effect is estimated based on averaging the predicted efficient score functions in the total sample.
#' Standard errors are based on asymptotic approximations using the estimated variance of the (estimated) efficient score functions.
#' @return A \code{dyntreatDML} object contains ten components, \code{effect}, \code{se}, \code{pval}, \code{ntrimmed}, \code{meantreat}, \code{meancontrol}, \code{psd1treat}, \code{psd2treat}, \code{psd1control}, and \code{psd2control} :
#' @return \code{effect}: estimate of the average effect of the treatment sequence.
#' @return \code{se}: standard error of the effect estimate.
#' @return \code{pval}: p-value of the effect estimate.
#' @return \code{ntrimmed}: number of discarded (trimmed) observations due to low products of propensity scores.
#' @return \code{meantreat}: Estimate of the mean potential outcome under the treatment sequence.
#' @return \code{meancontrol}: Estimate of the mean potential outcome under the control sequence.
#' @return \code{psd1treat}: P-score estimates for first treatment in treatment sequence.
#' @return \code{psd2treat}: P-score estimates for second treatment in treatment sequence.
#' @return \code{psd1control}: P-score estimates for first treatment in control sequence.
#' @return \code{psd2control}: P-score estimates for second treatment in control sequence.
#' @references Bodory, H., Huber, M., Laffers, L. (2020): "Evaluating (weighted) dynamic treatment effects by double machine learning", working paper, arXiv preprint arXiv:2012.00370.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @references van der Laan, M., Polley, E., Hubbard, A. (2007): "Super Learner", Statistical Applications in Genetics and Molecular Biology, 6.
#' @examples # A little example with simulated data (2000 observations)
#' \dontrun{
#' n=2000
#' # sample size
#' p0=10
#' # number of covariates at baseline
#' s0=5
#' # number of covariates that are confounders at baseline
#' p1=10
#' # number of additional covariates in period 1
#' s1=5
#' # number of additional covariates that are confounders in period 1
#' x0=matrix(rnorm(n*p0),ncol=p0)
#' # covariate matrix at baseline
#' beta0=c(rep(0.25,s0), rep(0,p0-s0))
#' # coefficients determining degree of confounding for baseline covariates
#' d1=(x0%*%beta0+rnorm(n)>0)*1
#' # equation of first treatment in period 1
#' x1=matrix(rnorm(n*p1),ncol=p1)
#' # covariate matrix for covariates of period 1
#' beta1=c(rep(0.25,s1), rep(0,p1-s1))
#' # coefficients determining degree of confounding for additonal covariates of period 1
#' d2=(x0%*%beta0+x1%*%beta1+0.5*d1+rnorm(n)>0)*1
#' # equation of second treatment in period 2
#' y2=x0%*%beta0+x1%*%beta1+1*d1+0.5*d2+rnorm(n)
#' # outcome equation in period 2
#' output=dyntreatDML(y2=y2,d1=d1,d2=d2,x0=x0,x1=x1,
#'        d1treat=1,d2treat=1,d1control=0,d2control=0)
#' cat("dynamic ATE: ",round(c(output$effect),3),", standard error: ",
#'     round(c(output$se),3), ", p-value: ",round(c(output$pval),3))
#' output$ntrimmed
#' # The true effect of the treatment sequence is 1.5}

#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm dnorm quantile coef fitted gaussian median
#' @import SuperLearner glmnet ranger xgboost e1071 mvtnorm
#' @export

dyntreatDML=function(y2,d1,d2,x0,x1, s=NULL, d1treat=1, d2treat=1, d1control=0, d2control=0,  trim=0.01, MLmethod="lasso", fewsplits=FALSE, normalized=TRUE){
  if (length(d1treat)==1) {d1tre=1*(d1==d1treat)} else {d1tre=d1treat}
  if (length(d2treat)==1) {d2tre=1*(d2==d2treat)} else {d2tre=d2treat}
  if (length(d1control)==1) {d1con=1*(d1==d1control)} else {d1con=d1control}
  if (length(d2control)==1) {d2con=1*(d2==d2control)} else {d2con=d2control}
  scorestreat=hddyntreat(y2=y2,d1=d1tre,d2=d2tre,x0=x0,x1=x1, s=s, trim=trim, MLmethod=MLmethod, fewsplits=fewsplits)
  scorescontrol=hddyntreat(y2=y2,d1=d1con,d2=d2con,x0=x0,x1=x1, s=s, trim=trim, MLmethod=MLmethod, fewsplits=fewsplits)
  trimmed=1*(scorescontrol[,10]+scorestreat[,10]>0)        #number of trimmed observations
  scorestreat=scorestreat[trimmed==0,]
  scorescontrol=scorescontrol[trimmed==0,]
  if (normalized==FALSE){
    tscores=(scorestreat[,1]*scorestreat[,2]*scorestreat[,3]*(scorestreat[,4]-scorestreat[,5])/(scorestreat[,6]*scorestreat[,7])+scorestreat[,1]*scorestreat[,2]*(scorestreat[,5]-scorestreat[,8])/scorestreat[,6]+scorestreat[,9]*scorestreat[,8])/mean(scorestreat[,9])
    cscores=(scorescontrol[,1]*scorescontrol[,2]*scorescontrol[,3]*(scorescontrol[,4]-scorescontrol[,5])/(scorescontrol[,6]*scorescontrol[,7])+scorescontrol[,1]*scorescontrol[,2]*(scorescontrol[,5]-scorescontrol[,8])/scorescontrol[,6]+scorescontrol[,9]*scorescontrol[,8])/mean(scorescontrol[,9])
  }
  if (normalized!=FALSE){
    ntreat=nrow(scorestreat)
    weightsumtreat1=sum(scorestreat[,1]*scorestreat[,2]*scorestreat[,3]/(scorestreat[,6]*scorestreat[,7]))
    weightsumtreat2=sum(scorestreat[,1]*scorestreat[,2]/scorestreat[,6])
    tscores=(ntreat*scorestreat[,1]*scorestreat[,2]*scorestreat[,3]*(scorestreat[,4]-scorestreat[,5])/(scorestreat[,6]*scorestreat[,7]))/weightsumtreat1+(ntreat*scorestreat[,1]*scorestreat[,2]*(scorestreat[,5]-scorestreat[,8])/scorestreat[,6])/weightsumtreat2+(scorestreat[,9]*scorestreat[,8])/mean(scorestreat[,9])
    ncontrol=nrow(scorescontrol)
    weightsumcontrol1=sum(scorescontrol[,1]*scorescontrol[,2]*scorescontrol[,3]/(scorescontrol[,6]*scorescontrol[,7]))
    weightsumcontrol2=sum(scorescontrol[,1]*scorescontrol[,2]/scorescontrol[,6])
    cscores=(ncontrol*scorescontrol[,1]*scorescontrol[,2]*scorescontrol[,3]*(scorescontrol[,4]-scorescontrol[,5])/(scorescontrol[,6]*scorescontrol[,7]))/weightsumcontrol1+(ncontrol*scorescontrol[,1]*scorescontrol[,2]*(scorescontrol[,5]-scorescontrol[,8])/scorescontrol[,6])/weightsumcontrol2+(scorescontrol[,9]*scorescontrol[,8])/mean(scorescontrol[,9])
  }
  meantreat=mean(tscores)
  meancontrol=mean(cscores)
  effect=meantreat - meancontrol
  se=sqrt(mean((tscores-cscores-effect)^2)/length(tscores))
  pval= 2*pnorm((-1)*abs(effect/se))
  list(effect=effect, se=se, pval=pval, ntrimmed=sum(trimmed), meantreat=meantreat, meancontrol=meancontrol, psd1treat=scorestreat[,6], psd2treat=scorestreat[,7], psd1control=scorescontrol[,6], psd2control=scorescontrol[,7])
}
