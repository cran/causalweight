#' Continuous Difference-in-Differences using Double Machine Learning for Panel Data
#' @description This function estimates the average treatment effect on the treated of a continuously distributed treatment in panel data based on a Difference-in-Differences (DiD) approach using double machine learning to control for time-varying confounders in a data-driven manner. It supports estimation under various machine learning methods and uses k-fold cross-fitting.
#' @param ydiff Outcome difference between two pre- and post-treatment periods. Should not contain missing values.
#' @param d Treatment variable in the treatment period of interest. Should be continuous and not contain missing values.
#' @param t Time variable indicating outcome periods. Should not contain missing values.
#' @param dtreat Value of the treatment under treatment (in the treatment period of interest). This value would be 1 for binary treatments.
#' @param dcontrol Value of the treatment under control (in the treatment period of interest). This value would be 0 for binary treatments.
#' @param t1 Value indicating the post-treatment outcome period in which the effect is evaluated, which is the later of the two periods used to generate the outcome difference in \code{ydiff}. For instance, if the pre-treatment outcome is measured in period 0 and the post-treatment outcome is measured in period 1 to generate \code{ydiff}, then \code{t1} is equal to 1. Default is 1.
#' @param controls Covariates and/or previous treatment history to be controlled for. Should not contain missing values.
#' @param MLmethod Machine learning method for estimating nuisance parameters using the \code{SuperLearner} package. Must be one of \code{"lasso"} (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"}, \code{"ensemble"}, or \code{"parametric"}.
#' @param psmethod Method for computing generalized propensity scores. Set to 1 for estimating conditional treatment densities using the treatment as dependent variable, or 2 for using the treatment kernel weights as dependent variable. Default is 1.
#' @param trim Trimming threshold (in percentage) for discarding observations with too much influence within any subgroup defined by the treatment group and time. Default is 0.1.
#' @param lognorm Logical indicating if log-normal transformation should be applied when estimating conditional treatment densities using the treatment as dependent variable. Default is FALSE.
#' @param bw Bandwidth for kernel density estimation. Default is NULL, implying that the bandwidth is calculated based on the rule-of-thumb.
#' @param bwfactor Factor by which the bandwidth is multiplied. Default is 0.7 (undersmoothing).
#' @param cluster Optional clustering variable for calculating standard errors.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @details This function estimates the Average Treatment Effect on the Treated (ATET) by Difference-in-Differences in panel data while controlling for confounders using double machine learning. The function supports different machine learning methods for estimating nuisance parameters and performs k-fold cross-fitting to improve estimation accuracy. The function also handles binary and continuous outcomes, and provides options for trimming and bandwidth adjustments in kernel density estimation.
#' @return A list with the following components:
#' @return \code{ATET}: Estimate of the Average Treatment Effect on the Treated.
#' @return \code{se}: Standard error of the ATET estimate.
#' @return \code{trimmed}: Number of discarded (trimmed) observations.
#' @return \code{pval}: P-value.
#' @return \code{pscores}: Propensity scores of untrimmed observations (2 columns): under treatment, under control.
#' @return \code{outcomepred}: Conditional outcome predictions of untrimmed observations.
#' @return \code{treat}: Treatment status of untrimmed observations.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @references Haddad, M., Huber, M., Medina-Reyes, J., Zhang, L. (2024): "Difference-in-Differences with Time-varying Continuous Treatments using Double/Debiased Machine Learning", working paper, University of Fribourg.
#' @examples
#' \dontrun{
#' # Example with simulated data
#' n=1000
#' x=0.5*rnorm(n)
#' u=runif(n,0,2)
#' d=x+u+rnorm(n)
#' y0=u+rnorm(n)
#' y1=2*d+x+u+rnorm(n)
#' t=rep(1,n)
#' # true effect is 2
#' results=didcontDMLpanel(ydiff=y1-y0, d=d, t=t, dtreat=1, dcontrol=0, controls=x, MLmethod="lasso")
#' cat("ATET: ", round(results$ATET, 3), ", Standard error: ", round(results$se, 3))
#' }
#' @importFrom stats rnorm lm predict sd dnorm
#' @importFrom SuperLearner SuperLearner
#' @import np sandwich
#' @export
didcontDMLpanel<-function(ydiff, d, t, dtreat, dcontrol, t1=1, controls, MLmethod="lasso", psmethod=1, trim=0.1, lognorm=FALSE, bw=NULL, bwfactor=0.7, cluster=NULL, k=3){
ybin=1*(length(unique(ydiff))==2 & min(ydiff)==0 & max(ydiff)==1)  # check if binary outcome
controls=data.frame(controls)
y=ydiff[t==t1]; d=d[t==t1]; controls=controls[t==t1,]   # drop periods that are not used for estimation
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
if(psmethod==2){                   # estimate densities based on kernel weights
  wtreat=kernwgtdtreat
  wcontrol=kernwgtdcontrol
}
stepsize=ceiling((1/k)*length(d))                               # sample size in folds
set.seed(1); idx= sample(length(d), replace=FALSE)              # shuffle data
param=c();
for (i in 1:k){                                                         # start of cross-fitting loop
  tesample=idx[((i-1)*stepsize+1):(min((i)*stepsize,length(d)))]
  trsample=idx[-tesample]                                                                # cross-fitting loop
  dcontrols=data.frame(d,controls); controls1=data.frame(1,controls)
  mut1=MLfunct(y=y[trsample], x=dcontrols[trsample,], MLmethod=MLmethod, ybin=ybin) # outcome model
  d1controlstest=data.frame(dtreat,controls)[tesample, ]     # test data under treatment dtreat
  colnames(d1controlstest)[1]="d"
  d0controlstest=data.frame(dcontrol,controls)[tesample, ]   # test data under non-treatment dcontrol
  colnames(d0controlstest)[1]="d"
  mucontrol=predict(mut1, d0controlstest, onlySL = TRUE)$pred #outcome prediction under non-treatment, post-treatment period
  if(psmethod!=2){
    ggg1=MLfunct(y=dd[trsample], x=controls1[trsample,], MLmethod=MLmethod,  ybin=0)
    pred1=predict(ggg1, controls1[tesample,], onlySL = TRUE)$pred
    resid1=dd[tesample]-pred1
    if(lognorm==TRUE){
      pscontrol=(dnorm( (log(dcontrol)-pred1)/sqrt(mean(resid1^2)))/dcontrol)  # generalized treatment propensity scores
      pstreat=(dnorm( (log(dtreat)-pred1)/sqrt(mean(resid1^2)))/dtreat)
    }
    if(lognorm==FALSE){
      pscontrol=(dnorm( (dcontrol-pred1)/sqrt(mean(resid1^2))))          # generalized treatment propensity scores
      pstreat=(dnorm( (dtreat-pred1)/sqrt(mean(resid1^2))))
    }
  }
  if(psmethod==2){
    pstreat=MLfunct(y=wtreat[trsample], x=controls1[trsample,], MLmethod=MLmethod,  ybin=0)
    pstreat=predict(pstreat, controls1[tesample,], onlySL = TRUE)$pred
    pscontrol=MLfunct(y=wcontrol[trsample], x=controls1[trsample,], MLmethod=MLmethod,  ybin=0)
    pscontrol=predict(pscontrol, controls1[tesample,], onlySL = TRUE)$pred
  }
  param=rbind(param, cbind(pstreat,pscontrol, mucontrol))
}  # close cross-fitting loop
param = param[order(idx),]                              # sort nuisance parameters according to original order of observations
param=cbind(kernwgtdtreat, kernwgtdcontrol,  y, param)
resd1=param[,1]/sum(param[,1])                           # normalize weights
resd0=(param[,2]*param[,4]/param[,5])/sum(param[,2]*param[,4]/param[,5])
trimmed=(resd1>trim | resd0>trim)     # identify observations with too much influence in one of the subgroups which are to be dropped
param=param[trimmed==0,]
reg=resd1[trimmed==0]*(param[,3]-param[,6])/sum(resd1[trimmed==0])    #renormalize weights and multiply
resd0=resd0[trimmed==0]*(param[,3]-param[,6])/sum(resd0[trimmed==0])
score=sum(1-trimmed)*(reg-resd0) # estimate the score function
ATET=mean(score)        # ATET
kfunct=param[,1]/bw
meank=mean(kfunct)
fd=sum(kfunct)/length(kfunct)
scoremod=score-ATET/fd*(kfunct-meank)
if (is.null(cluster)){
  se=summary(lm(scoremod~1))$coefficients[, "Std. Error"]     # se without clustering
}else{
  se=sqrt(vcovCL(lm(scoremod~1), cluster = cluster[trimmed==0]))[1,1]}    # se with clustering
pval= 2*pnorm((-1)*abs(ATET/se))                           # p-value
list(ATET=ATET, se=se, pval=pval, ntrimmed=sum(trimmed), pscores=param[,4:5], outcomepred=param[,6], treat=d[trimmed==0])
}
