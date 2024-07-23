#' Continuous Difference-in-Differences using Double Machine Learning for Repeated Cross-Sections
#' @description This function estimates the average treatment effect on the treated of a continuously distributed treatment in repeated cross-sections based on a Difference-in-Differences (DiD) approach using double machine learning to control for time-varying confounders in a data-driven manner. It supports estimation under various machine learning methods and uses k-fold cross-fitting.
#' @param y Outcome variable. Should not contain missing values.
#' @param d Treatment variable in the treatment period of interest. Should be continuous and not contain missing values.
#' @param t Time variable indicating outcome periods. Should not contain missing values.
#' @param dtreat Value of the treatment under treatment (in the treatment period of interest). This value would be 1 for binary treatments.
#' @param dcontrol Value of the treatment under control (in the treatment period of interest). This value would be 0 for binary treatments.
#' @param t0 Value indicating the pre-treatment outcome period. Default is 0.
#' @param t1 Value indicating the post-treatment outcome period in which the effect is evaluated. Default is 1.
#' @param controls Covariates and/or previous treatment history to be controlled for. Should not contain missing values.
#' @param MLmethod Machine learning method for estimating nuisance parameters using the \code{SuperLearner} package. Must be one of \code{"lasso"} (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"}, \code{"ensemble"}, or \code{"parametric"}.
#' @param psmethod Method for computing generalized propensity scores. Set to 1 for estimating conditional treatment densities using the treatment as dependent variable, or 2 for using the treatment kernel weights as dependent variable. Default is 1.
#' @param trim Trimming threshold (in percentage) for discarding observations with too much influence within any subgroup defined by the treatment group and time. Default is 0.1.
#' @param lognorm Logical indicating if log-normal transformation should be applied when estimating conditional treatment densities using the treatment as dependent variable. Default is FALSE.
#' @param bw Bandwidth for kernel density estimation. Default is NULL, implying that the bandwidth is calculated based on the rule-of-thumb.
#' @param bwfactor Factor by which the bandwidth is multiplied. Default is 0.7 (undersmoothing).
#' @param cluster Optional clustering variable for calculating standard errors.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @details This function estimates the Average Treatment Effect on the Treated (ATET) by Difference-in-Differences in repeated cross-sections while controlling for confounders using double machine learning. The function supports different machine learning methods for estimating nuisance parameters and performs k-fold cross-fitting to improve estimation accuracy. The function also handles binary and continuous outcomes, and provides options for trimming and bandwidth adjustments in kernel density estimation.
#' @return A list with the following components:
#' @return \code{ATET}: Estimate of the Average Treatment Effect on the Treated.
#' @return \code{se}: Standard error of the ATET estimate.
#' @return \code{trimmed}: Number of discarded (trimmed) observations.
#' @return \code{pval}: P-value.
#' @return \code{pscores}: Propensity scores (4 columns): under treatment in period t1, under treatment in period t0, under control in period t1, under control in period t0.
#' @return \code{outcomes}: Conditional outcomes (3 columns): in treatment group in period t0, in control group in period t1, in control group in period t0.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @references Haddad, M., Huber, M., Medina-Reyes, J., Zhang, L. (2024): "Difference-in-Differences under time-varying continuous treatments based on double machine learning"
#' @examples
#' \dontrun{
#' # Example with simulated data
#' n=2000
#' t=rep(c(0, 1), each=n/2)
#' x=0.5*rnorm(n)
#' u=runif(n,0,2)
#' d=x+u+rnorm(n)
#' y=(2*d+x)*t+u+rnorm(n)
#' # true effect is 2
#' results=didcontDML(y=y, d=d, t=t, dtreat=1, dcontrol=0, controls=x, MLmethod="lasso")
#' cat("ATET: ", round(results$ATET, 3), ", Standard error: ", round(results$se, 3))
#' }
#' @importFrom stats rnorm lm predict sd dnorm
#' @importFrom SuperLearner SuperLearner
#' @import np sandwich
#' @export
didcontDML=function(y, d, t, dtreat, dcontrol, t0=0, t1=1, controls, MLmethod="lasso", psmethod=1, trim=0.1, lognorm=FALSE, bw=NULL, bwfactor=0.7, cluster=NULL, k=3) {
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  controls=data.frame(controls)
  y=y[t==t0 | t==t1]; d=d[t==t0 | t==t1]; controls=controls[(t==t0 | t==t1),]; t=t[t==t0 | t==t1]  # drop periods that are not used for estimation
  if(is.null(bw)) bw=sd(d)*2.34/(length(d)^0.25)          # rule-of-thumb bandwidth for treatment based on its standard deviation
  bw=bw*bwfactor                                          # change bandwidth according to bwfactor
  kernwgtdcontrol=npksum(bws=bw, txdat = d, tydat = y, exdat = dcontrol, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  kernwgtdtreat=npksum(bws=bw, txdat = d, tydat = y, exdat = dtreat, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  kernwgtdtreat=kernwgtdtreat/sum(kernwgtdtreat); kernwgtdcontrol=kernwgtdcontrol/sum(kernwgtdcontrol) # normalized kernel weights for treatment comparison
  dd=d
  if(psmethod!=2 & lognorm==TRUE){    #lognormal transformation when estimating joint densities indirectly
    dd[d==0]=0.00001; dd=log(dd)
    if (dcontrol==0) d0=0.00001; if (dtreat==0) d1=0.00001
  }
  if(psmethod==2){                   # estimate joint densities directly
    wtreatt1=kernwgtdtreat*(t==t1)
    wtreatt0=kernwgtdtreat*(t==t0)
    wcontrolt1=kernwgtdcontrol*(t==t1)
    wcontrolt0=kernwgtdcontrol*(t==t0)
  }
  stepsize=ceiling((1/k)*length(d))                               # sample size in folds
  set.seed(1); idx= sample(length(d), replace=FALSE)              # shuffle data
  param=c();
  for (i in 1:k){                                                         # start of cross-fitting loop
    tesample=idx[((i-1)*stepsize+1):(min((i)*stepsize,length(d)))]
    trsample=idx[-tesample]                                                                # cross-fitting loop
    dcontrols=data.frame(d,controls); controls1=data.frame(1,controls)
    mut1=MLfunct(y=y[trsample], x=dcontrols[trsample,], d1=1*(t[trsample]==t1), MLmethod=MLmethod, ybin=ybin) # outcome model in post-treatment period
    mut0=MLfunct(y=y[trsample], x=dcontrols[trsample,], d1=1*(t[trsample]==t0), MLmethod=MLmethod, ybin=ybin) # outcome model in pre-treatment period
    d1controlstest=data.frame(dtreat,controls)[tesample, ]     # test data under treatment dtreat
    colnames(d1controlstest)[1]="d"
    d0controlstest=data.frame(dcontrol,controls)[tesample, ]   # test data under non-treatment dcontrol
    colnames(d0controlstest)[1]="d"
    mucontrolt1=predict(mut1, d0controlstest, onlySL = TRUE)$pred #outcome prediction under non-treatment, post-treatment period
    mutreatt0=predict(mut0, d1controlstest, onlySL = TRUE)$pred  #outcome prediction under treatment, pre-treatment period
    mucontrolt0=predict(mut0, d0controlstest, onlySL = TRUE)$pred   #outcome prediction under non-treatment, pre-treatment period
    if(psmethod!=2){                                                # estimate joint densities indirectly
      ggg1=MLfunct(y=dd[trsample], x=controls1[trsample,], d1=1*(t[trsample]==t1), MLmethod=MLmethod,  ybin=0)
      pred1=predict(ggg1, controls1[tesample,], onlySL = TRUE)$pred
      resid1=dd[tesample]-pred1
      ggg0=MLfunct(y=dd[trsample], x=controls1[trsample,], d1=1*(t[trsample]==t0), MLmethod=MLmethod,  ybin=0)
      pred0=predict(ggg0, controls1[tesample,], onlySL = TRUE)$pred
      resid0=dd[tesample]-pred0
      if(lognorm==TRUE){
        pscontrol1=(dnorm( (log(dcontrol)-pred1)/sqrt(mean(resid1^2)))/dcontrol)  # generalized treatment propensity scores conditional on t1
        pstreat1=(dnorm( (log(dtreat)-pred1)/sqrt(mean(resid1^2)))/dtreat)
        pscontrol0=(dnorm( (log(dcontrol)-pred0)/sqrt(mean(resid0^2)))/dcontrol)  # generalized treatment propensity scores conditional on t0
        pstreat0=(dnorm( (log(dtreat)-pred0)/sqrt(mean(resid0^2)))/dtreat)
      }
      if(lognorm==FALSE){
        pscontrol1=(dnorm( (dcontrol-pred1)/sqrt(mean(resid1^2))))          # generalized treatment propensity scores conditional on t1
        pstreat1=(dnorm( (dtreat-pred1)/sqrt(mean(resid1^2))))
        pscontrol0=(dnorm( (dcontrol-pred0)/sqrt(mean(resid0^2))))          # generalized treatment propensity scores conditional on t0
        pstreat0=(dnorm( (dtreat-pred0)/sqrt(mean(resid0^2))))
      }
      posttreat=1*(t==t1)
      gg=MLfunct(y=posttreat[trsample], x=controls1[trsample,], MLmethod=MLmethod,  ybin=1)
      predpost=predict(gg, controls1[tesample,], onlySL = TRUE)$pred # conditional probability of post-treatment period
      rhotreatt0=pstreat0*(1-predpost)                               # compute the joint densities based on the conditional probabilities/densities of period and treatment
      rhotreatt1=pstreat1*predpost
      rhocontrolt0=pscontrol0*(1-predpost)
      rhocontrolt1=pscontrol1*predpost
    }
    if(psmethod==2){                                                 # estimate joint densities directly
      rhotreatt0=MLfunct(y=wtreatt0[trsample], x=controls1[trsample,], MLmethod=MLmethod,  ybin=0)
      rhotreatt0=predict(rhotreatt0, controls1[tesample,], onlySL = TRUE)$pred
      rhotreatt1=MLfunct(y=wtreatt1[trsample], x=controls1[trsample,], MLmethod=MLmethod,  ybin=0)
      rhotreatt1=predict(rhotreatt1, controls1[tesample,], onlySL = TRUE)$pred
      rhocontrolt0=MLfunct(y=wcontrolt0[trsample], x=controls1[trsample,], MLmethod=MLmethod,  ybin=0)
      rhocontrolt0=predict(rhocontrolt0, controls1[tesample,], onlySL = TRUE)$pred
      rhocontrolt1=MLfunct(y=wcontrolt1[trsample], x=controls1[trsample,], MLmethod=MLmethod,  ybin=0)
      rhocontrolt1=predict(rhocontrolt1, controls1[tesample,], onlySL = TRUE)$pred
    }
    param=rbind(param, cbind(rhotreatt1,rhotreatt0,rhocontrolt1,rhocontrolt0,mutreatt0,mucontrolt1,mucontrolt0))
  }  # close cross-fitting loop
  param = param[sort(idx),]                              # sort nuisance parameters according to original order of observations
  param=cbind(kernwgtdtreat, kernwgtdcontrol,  1*(t==t1), y, param)
  resd1t1=(param[,1]*param[,3])/sum(param[,1]*param[,3])    # normalize weights to add up to one
  resd1t0=(param[,1]*(1-param[,3])*param[,5]/param[,6])/sum(param[,1]*(1-param[,3])*param[,5]/param[,6])
  resd0t1=(param[,2]*param[,3]*param[,5]/param[,7])/sum(param[,2]*param[,3]*param[,5]/param[,7])
  resd0t0=(param[,2]*(1-param[,3])*param[,5]/param[,8])/sum(param[,2]*(1-param[,3])*param[,5]/param[,8])
  trimmed=(resd1t1>trim | resd1t0>trim | resd0t1>trim | resd0t0>trim)     # identify observations with too much influence in one of the subgroups which are to be dropped
  param=param[trimmed==0,]
  reg=resd1t1[trimmed==0]*(param[,4]-param[,9]-param[,10]+param[,11])/sum(resd1t1[trimmed==0])   # renormalize weights to add up to one and multiply
  resd1t0=resd1t0[trimmed==0]*(param[,4]-param[,9])/sum(resd1t0[trimmed==0])
  resd0t1=resd0t1[trimmed==0]*(param[,4]-param[,10])/sum(resd0t1[trimmed==0])
  resd0t0=resd0t0[trimmed==0]*(param[,4]-param[,11])/sum(resd0t0[trimmed==0])
  score=sum(1-trimmed)*(reg-resd1t0-resd0t1+resd0t0) # estimate the score function
  ATET=mean(score)        # ATET
  kfunct=param[,1]/bw*param[,3]
  fd=mean(kfunct)
  scoremod=score-ATET/fd*(kfunct-fd)
  if (is.null(cluster)){
    se=summary(lm(scoremod~1))$coefficients[, "Std. Error"]     # se without clustering
  }else{
    se=sqrt(vcovCL(lm(scoremod~1), cluster = cluster[trimmed==0]))[1,1]}    # se with clustering
  pval= 2*pnorm((-1)*abs(ATET/se))                           # p-value
  list(ATET=ATET, se=se, pval=pval, ntrimmed=sum(trimmed), pscores=param[,5:8], outcomes=param[,9:11])
}
