#' Difference-in-Differences in Repeated Cross-Sections for Binary Treatments using Double Machine Learning
#' @description This function estimates the average treatment effect on the treated (ATET) in the post-treatment period for a binary treatment using a doubly robust Difference-in-Differences (DiD) approach for repeated cross-sections that is combined with double machine learning. It controls for (possibly time-varying) confounders in a data-driven manner and supports various machine learning methods for estimating nuisance parameters through k-fold cross-fitting.
#' @param y Outcome variable. Should not contain missing values.
#' @param d Treatment group indicator (binary). Should not contain missing values.
#' @param t Time period indicator (binary). Should be 1 for post-treatment period and 0 for pre-treatment period. Should not contain missing values.
#' @param x Covariates to be controlled for. Should not contain missing values.
#' @param MLmethod Machine learning method for estimating nuisance parameters using the \code{SuperLearner} package. Must be one of \code{"lasso"} (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"}, \code{"ensemble"}, or \code{"parametric"}.
#' @param est Estimation method. Must be one of \code{"dr"} (default) for doubly robust, \code{"ipw"} for inverse probability weighting (not doubly robust!), or \code{"reg"} for regression (not doubly robust!).
#' @param trim Trimming threshold (in percentage) for discarding observations with too small propensity scores within any subgroup defined by the treatment group and time. Default is 0.05.
#' @param cluster Optional clustering variable for calculating cluster-robust standard errors.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @details This function estimates the Average Treatment Effect on the Treated (ATET) in the post-treatment period based on Difference-in-Differences in repeated cross-sections when controlling for confounders in a data-adaptive manner using double machine learning. The function supports different machine learning methods to estimate nuisance parameters (conditional mean outcomes and propensity scores) as well as cross-fitting to mitigate overfitting. Besides double machine learning, the function also provides inverse probability weighting and regression adjustment methods (which are, however, not doubly robust).
#' @return A list with the following components:
#' @return \code{ATET}: Estimate of the Average Treatment Effect on the Treated (ATET) in the post-treatment period.
#' @return \code{se}: Standard error of the ATET estimate.
#' @return \code{pval}: P-value of the ATET estimate.
#' @return \code{trimmed}: Number of discarded (trimmed) observations.
#' @return \code{pscores}: Propensity scores of untrimmed observations (4 columns): under treatment in period 1, under treatment in period 0, under control in period 1, under control in period 0.
#' @return \code{outcomepred}: Conditional outcome predictions of untrimmed observations (3 columns): in treatment group in period 0, in control group in period 1, in control group in period 0.
#' @return \code{treat}: Treatment status of untrimmed observations.
#' @return \code{time}: Time period of untrimmed observations.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @references Zimmert, M. (2020): "Efficient difference-in-differences estimation with high-dimensional common trend confounding", arXiv preprint 1809.01643.
#' @examples
#' \dontrun{
#' # Example with simulated data
#' n=4000                            # sample size
#' t=1*(rnorm(n)>0)                  # time period
#' u=runif(n,0,1)                    # time constant unobservable
#' x= 0.25*t+runif(n,0,1)            # time varying covariate
#' d=1*(x+u+2*rnorm(n)>0)            # treatment
#' y=d*t+t+x+u+2*rnorm(n)            # outcome
#' # true effect is equal to 1
#' results=didDML(y=y, d=d, t=t, x=x)
#' cat("ATET: ", round(results$ATET, 3), ", Standard error: ", round(results$se, 3))
#' }
#' @importFrom stats rnorm lm predict sd dnorm
#' @importFrom SuperLearner SuperLearner
#' @import sandwich
#' @export
didDML<-function(y, d, t, x, MLmethod="lasso", est="dr",  trim=0.05, cluster=NULL, k=3){
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  controls=data.frame(x)
  stepsize=ceiling((1/k)*length(d))                               # sample size in folds
  set.seed(1); idx= sample(length(d), replace=FALSE)              # shuffle data
  param=c();
  for (i in 1:k){                                                         # start of cross-fitting loop
    tesample=idx[((i-1)*stepsize+1):(min((i)*stepsize,length(d)))]
    trsample=idx[-tesample]                                                                # cross-fitting loop
    ytr=y[trsample]; dtr=d[trsample]; ttr=t[trsample]
    controlstr=data.frame(1,controls)[trsample,];
    mud0t1=MLfunct(y=ytr, x=controlstr, d1=1*(dtr==0 & ttr==1), MLmethod=MLmethod, ybin=ybin) # outcome model under non-treatment in post-treatment period
    mud1t0=MLfunct(y=ytr, x=controlstr, d1=1*(dtr==1 & ttr==0), MLmethod=MLmethod, ybin=ybin) # outcome model under treatment in pre-treatment period
    mud0t0=MLfunct(y=ytr, x=controlstr, d1=1*(dtr==0 & ttr==0), MLmethod=MLmethod, ybin=ybin) # outcome model under treatment in pre-treatment period
    controlstest=data.frame(1,controls)[tesample, ]     # test data
    mud0t1=predict(mud0t1, controlstest, onlySL = TRUE)$pred #outcome prediction under non-treatment, post-treatment period
    mud1t0=predict(mud1t0, controlstest, onlySL = TRUE)$pred  #outcome prediction under treatment, pre-treatment period
    mud0t0=predict(mud0t0, controlstest, onlySL = TRUE)$pred   #outcome prediction under non-treatment, pre-treatment period
    rhod1t1=MLfunct(y=1*(dtr==1 & ttr==1), x=controlstr, MLmethod=MLmethod,  ybin=1)  # now estimation of and prediction based on propensity score models
    rhod1t1=predict(rhod1t1, controlstest, onlySL = TRUE)$pred
    rhod1t0=MLfunct(y=1*(dtr==1 & ttr==0), x=controlstr, MLmethod=MLmethod,  ybin=1)
    rhod1t0=predict(rhod1t0, controlstest, onlySL = TRUE)$pred
    rhod0t1=MLfunct(y=1*(dtr==0 & ttr==1), x=controlstr, MLmethod=MLmethod,  ybin=1)
    rhod0t1=predict(rhod0t1, controlstest, onlySL = TRUE)$pred
    rhod0t0=MLfunct(y=1*(dtr==0 & ttr==0), x=controlstr, MLmethod=MLmethod,  ybin=1)
    rhod0t0=predict(rhod0t0, controlstest, onlySL = TRUE)$pred
    param=rbind(param, cbind(rhod1t1,rhod1t0,rhod0t1,rhod0t0,mud1t0,mud0t1,mud0t0))
  }  # close cross-fitting loop
  param = param[order(idx),]                              # sort nuisance parameters according to original order of observations
  param=cbind(d, 1-d,  t, y, param)
  trimmed=1*( (param[,5]>(1-trim)) | (param[,6]<trim*d*(1-t)) | (param[,7]<trim*(1-d)*t) | (param[,8]<trim*(1-d)*(1-t)))     # identify observations with small propensity scores in subgroups which are dropped
  param=param[trimmed==0,]
  param[,6]=ifelse(param[,6] == 0, .Machine$double.eps, param[,6])
  param[,7]=ifelse(param[,7] == 0, .Machine$double.eps, param[,7])
  param[,8]=ifelse(param[,8] == 0, .Machine$double.eps, param[,8])
  resd1t1=(param[,1]*param[,3])/sum(param[,1]*param[,3])    # normalize weights to add up to one
  resd1t0=(param[,1]*(1-param[,3])*param[,5]/param[,6])/sum(param[,1]*(1-param[,3])*param[,5]/param[,6])
  resd0t1=(param[,2]*param[,3]*param[,5]/param[,7])/sum(param[,2]*param[,3]*param[,5]/param[,7])
  resd0t0=(param[,2]*(1-param[,3])*param[,5]/param[,8])/sum(param[,2]*(1-param[,3])*param[,5]/param[,8])
  reg=resd1t1*(param[,4]-param[,9]-param[,10]+param[,11])
  if (est=="dr"){
    resd1t0=resd1t0*(param[,4]-param[,9])
    resd0t1=resd0t1*(param[,4]-param[,10])
    resd0t0=resd0t0*(param[,4]-param[,11])
    score=sum(1-trimmed)*(reg-resd1t0-resd0t1+resd0t0) # estimate the score function
  }
  if (est=="reg") score=sum(1-trimmed)*(reg)
  if (est=="ipw") score=sum(1-trimmed)*(resd1t1*param[,4]-resd1t0*param[,4]-resd0t1*param[,4]+resd0t0*param[,4])
  ATET=mean(score)        # ATET
  if (is.null(cluster)){
    se=summary(lm(score~1))$coefficients[, "Std. Error"]     # se without clustering
  }else{
    se=sqrt(vcovCL(lm(score~1), cluster = cluster[trimmed==0]))[1,1]}    # se with clustering
  pval= 2*pnorm((-1)*abs(ATET/se))                           # p-value
  list(ATET=ATET, se=se, pval=pval, ntrimmed=sum(trimmed), pscores=param[,5:8], outcomepred=param[,9:11], treat=param[,1], time=param[,3])
}
