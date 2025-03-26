#' ATET Estimation for Binary Treatments using Double Machine Learning
#' @description This function estimates the average treatment effect on the treated (ATET) for a binary treatment. Combines estimation based on (doubly robust) efficient score functions with double machine learning to control for confounders in a data-driven way.
#' @param y Outcome variable. Should not contain missing values.
#' @param d Treatment group indicator (binary). Should not contain missing values.
#' @param x Covariates to be controlled for. Should not contain missing values.
#' @param MLmethod Machine learning method for estimating nuisance parameters using the \code{SuperLearner} package. Must be one of \code{"lasso"} (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"}, \code{"ensemble"}, or \code{"parametric"}.
#' @param est Estimation method. Must be one of \code{"dr"} (default) for doubly robust, \code{"ipw"} for inverse probability weighting (not doubly robust!), or \code{"reg"} for regression (not doubly robust!).
#' @param trim Trimming threshold (in percentage) for discarding observations with propensity scores too close to one. Default is 0.05, implying that treatment propensity scores larger than 1-0.05=0.95 (such that the probability to be not treated is below 0.05) are trimmed.
#' @param cluster Optional clustering variable for calculating cluster-robust standard errors.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @details This function estimates the Average Treatment Effect on the Treated (ATET) under conditional independence, assuming that confounders jointly affecting the treatment and the outcome can be controlled for by observed covariates. Estimation is based on the (doubly robust) efficient score functions for potential outcomes in combination with double machine learning with cross-fitting, see Chernozhukov et al (2018). To this end, one part of the data is used for estimating the model parameters of the treatment and outcome equations based machine learning. The other part of the data is used for predicting the efficient score functions. The roles of the data parts are swapped (using k-fold cross-fitting) and the ATET is estimated based on averaging the predicted efficient score functions in the total sample. Besides double machine learning, the function also provides inverse probability weighting and regression adjustment methods (which are, however, not doubly robust).
#' @return A list with the following components:
#' @return \code{ATET}: Estimate of the Average Treatment Effect on the Treated (ATET) in the post-treatment period.
#' @return \code{se}: Standard error of the ATET estimate.
#' @return \code{pval}: P-value of the ATET estimate.
#' @return \code{trimmed}: Number of discarded (trimmed) observations.
#' @return \code{treat}: Treatments of untrimmed observations.
#' @return \code{outcome}: Outcomes of untrimmed observations.
#' @return \code{pscores}: Treatment propensity scores of untrimmed observations.
#' @return \code{outcomepred}: Conditional outcome predictions under nontreatment of untrimmed observations.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @examples
#' \dontrun{
#' n=2000                            # sample size
#' p=50                              # number of covariates
#' s=2                               # number of covariates that are confounders
#' x=matrix(rnorm(n*p),ncol=p)       # covariate matrix
#' beta=c(rep(0.25,s), rep(0,p-s))   # coefficients determining degree of confounding
#' d=(x%*%beta+rnorm(n)>0)*1         # treatment equation
#' y=x%*%beta+0.5*d+rnorm(n)         # outcome equation
#' # The true ATE is equal to 0.5
#' output=ATETDML(y,d,x)
#' cat("ATET: ",round(c(output$ATET),3),", standard error: ",
#'     round(c(output$se),3), ", p-value: ",round(c(output$pval),3))
#' output$ntrimmed
#' }
#' @importFrom stats rnorm lm predict sd dnorm
#' @importFrom SuperLearner SuperLearner
#' @import sandwich
#' @export


ATETDML<-function(y, d, x, MLmethod="lasso", est="dr",  trim=0.05, cluster=NULL, k=3){
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  controls=data.frame(x)
  stepsize=ceiling((1/k)*length(d))                               # sample size in folds
  set.seed(1); idx= sample(length(d), replace=FALSE)              # shuffle data
  param=c();
  for (i in 1:k){                                                         # start of cross-fitting loop
    tesample=idx[((i-1)*stepsize+1):(min((i)*stepsize,length(d)))]
    trsample=idx[-tesample]                                                                # cross-fitting loop
    ytr=y[trsample]; dtr=d[trsample] 
    controlstr=data.frame(1,controls)[trsample,]; 
    mud0=MLfunct(y=ytr, x=controlstr, d1=1*(dtr==0), MLmethod=MLmethod, ybin=ybin) # outcome model under non-treatment 
    controlstest=data.frame(1,controls)[tesample, ]     # test data 
    mud0=predict(mud0, controlstest, onlySL = TRUE)$pred #outcome prediction under non-treatment
    rhod1=MLfunct(y=dtr, x=controlstr, MLmethod=MLmethod,  ybin=1)  # now estimation of and prediction based on propensity score models
    rhod1=predict(rhod1, controlstest, onlySL = TRUE)$pred
    param=rbind(param, cbind(rhod1,mud0))
  }  # close cross-fitting loop
  param = param[order(idx),]                              # sort nuisance parameters according to original order of observations
  param=cbind(d, y, param) 
  trimmed=1*( param[,3]>(1-trim) )     # identify observations with small propensity scores in subgroups which are dropped
  param=param[trimmed==0,]
  resd1=param[,1]/sum(param[,1])    # normalize weights to add up to one
  resd0=( (1-param[,1])*param[,3]/(1-param[,3]))/sum( (1-param[,1])*param[,3]/(1-param[,3]))
  reg=resd1*(param[,2]-param[,4])      
  if (est=="dr"){
    resd0=resd0*(param[,2]-param[,4])
    score=sum(1-trimmed)*(reg-resd0) # estimate the score function
  }
  if (est=="reg") score=sum(1-trimmed)*(reg)
  if (est=="ipw") score=sum(1-trimmed)*(resd1*param[,2]-resd0*param[,2])
  ATET=mean(score)        # ATET
  if (is.null(cluster)){
    se=summary(lm(score~1))$coefficients[, "Std. Error"]     # se without clustering
  }else{
    se=sqrt(vcovCL(lm(score~1), cluster = cluster[trimmed==0]))[1,1]}    # se with clustering 
  pval= 2*pnorm((-1)*abs(ATET/se))                           # p-value
  list(ATET=ATET, se=se, pval=pval, ntrimmed=sum(trimmed), treat=param[,1], outcome=param[,2], pscores=param[,3], outcomepred=param[,4])
}