#' Continuous Treatment Effect Estimation under Selection on Observables using Double Machine Learning
#' @description This function estimates the average treatment effect (ATE, default) or average treatment effect on the treated (ATET, optional) of a continuously distributed treatment under selection on observables (unconfoundedness), using double machine learning with kernel smoothing around the treatment values of interest.
#' @param y Outcome variable. Should not contain missing values.
#' @param d Treatment variable. Should be continuous and not contain missing values.
#' @param x Covariates. Should not contain missing values.
#' @param dtreat Value of the treatment for which the mean potential outcome is evaluated as the "treatment" level.
#' @param dcontrol Value of the treatment for which the mean potential outcome is evaluated as the "control" level.
#' @param ATET Logical. If \code{FALSE} (default), the average treatment effect (ATE) in the total population is estimated. If \code{TRUE}, the average treatment effect on units with treatment close to \code{dtreat} (ATET) is estimated instead.
#' @param MLmethod Machine learning method for estimating nuisance parameters using the \code{SuperLearner} package. Must be one of \code{"lasso"} (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"}, \code{"ensemble"}, or \code{"parametric"}.
#' @param psmethod Method for computing generalized propensity scores. Set to 1 for estimating conditional treatment densities using the treatment as dependent variable (imposing the parametric assumption of a Gaussian treatment distribution), or 2 for using the treatment kernel weights as dependent variable (nonparametric). Default is 1.
#' @param trim Trimming threshold (in percent) for discarding observations whose (normalized) weight would otherwise dominate the estimate. Default is 0.1.
#' @param lognorm Logical indicating if log-normal transformation should be applied when estimating conditional treatment densities using the treatment as dependent variable (parametric). Default is FALSE.
#' @param bw Bandwidth for kernel density estimation. Default is NULL, implying that the bandwidth is calculated based on the rule-of-thumb.
#' @param bwfactor Factor by which the bandwidth is multiplied. Default is 0.7 (undersmoothing).
#' @param cluster Optional clustering variable for calculating standard errors.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @details Under ATET=FALSE (default), the function estimates the ATE using a doubly robust, kernel-weighted score for continuous treatments, following Kennedy, Ma, McHugh, and Small (2017) and Colangelo and Lee (2020): for each of \code{dtreat} and \code{dcontrol}, every observation contributes both an outcome-regression prediction at that dose and a kernel/generalized-propensity-score-weighted residual correction. Under ATET=TRUE, a different doubly robust function instead makes observations near \code{dcontrol} to match the covariate distribution of units with treatment near \code{dtreat}, targeting the effect on that treated subpopulation rather than the total population.
#' @return A list with the following components:
#' @return \code{effect}: Estimate of the average treatment effect (ATE) or average treatment effect on the treated (ATET).
#' @return \code{se}: Standard error of the effect estimate.
#' @return \code{trimmed}: Number of discarded (trimmed) observations.
#' @return \code{pval}: P-value.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @references Kennedy, E. H., Ma, Z., McHugh, M. D., Small, D. S. (2017): "Non-parametric methods for doubly robust estimation of continuous treatment effects", Journal of the Royal Statistical Society Series B, 79, 1229-1245.
#' @references Colangelo, K., Lee, Y.-Y. (2020): "Double debiased machine learning nonparametric inference with continuous treatments", arXiv preprint 2004.03036.
#' @examples
#' \dontrun{
#' # Example with simulated data
#' n=3000
#' x=0.8*rnorm(n)
#' d=x+rnorm(n)
#' y=2*d+x+rnorm(n)
#' # true effect is 2
#' results=treatcontDML(y=y, d=d, dtreat=1, dcontrol=0, x=x, MLmethod="lasso")
#' cat("ATE: ", round(results$effect, 3), ", Standard error: ", round(results$se, 3))
#' }
#' @importFrom stats rnorm lm predict sd dnorm pnorm
#' @import np sandwich
#' @export
treatcontDML=function(y, d, x, dtreat, dcontrol, ATET=FALSE, MLmethod="lasso", psmethod=1, trim=0.1, lognorm=FALSE, bw=NULL, bwfactor=0.7, cluster=NULL, k=3){
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  x=data.frame(x)
  if(is.null(bw)) bw=sd(d)*2.34/(length(d)^0.25)          # rule-of-thumb bandwidth for treatment based on its standard deviation
  bw=bw*bwfactor                                          # change bandwidth according to bwfactor
  kernwgtdcontrol=npksum(bws=bw, txdat = d, tydat = y, exdat = dcontrol, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  kernwgtdtreat=npksum(bws=bw, txdat = d, tydat = y, exdat = dtreat, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  kernwgtdtreatnorm=kernwgtdtreat/sum(kernwgtdtreat); kernwgtdcontrolnorm=kernwgtdcontrol/sum(kernwgtdcontrol) # normalized kernel weights (sum to 1), used for ATET and trimming
  dd=d
  if(psmethod!=2 & lognorm==TRUE){    #lognormal transformation when estimating joint densities indirectly
    dd[d==0]=0.00001; dd=log(dd)
    if (dcontrol==0) dcontrol=0.00001; if (dtreat==0) dtreat=0.00001  # fix: modify dcontrol/dtreat themselves (not unused d0/d1)
  }
  if(psmethod==2){                   # estimate densities based on kernel weights
    wtreat=kernwgtdtreat
    wcontrol=kernwgtdcontrol
  }
  set.seed(1); idx = sample(length(d), replace=FALSE)   # shuffle data
  folds = split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))  # generate folds
  param=c();
  for (i in 1:k){                                                         # start of cross-fitting loop
    tesample=folds[[i]]
    trsample=idx[!(idx %in% tesample)]                                    # cross-fitting loop
    dcontrols=data.frame(d,x); x1=data.frame(1,x)
    mut1=MLfunct(y=y[trsample], x=dcontrols[trsample,], MLmethod=MLmethod, ybin=ybin) # outcome model E[Y|D,X]
    d1xtest=data.frame(dtreat,x)[tesample, ]                   # test data under treatment dtreat
    colnames(d1xtest)[1]="d"
    d0xtest=data.frame(dcontrol,x)[tesample, ]                 # test data under non-treatment dcontrol
    colnames(d0xtest)[1]="d"
    mutreat=predict(mut1, d1xtest, onlySL = TRUE)$pred          # outcome prediction under treatment (needed for ATE)
    mucontrol=predict(mut1, d0xtest, onlySL = TRUE)$pred        # outcome prediction under non-treatment
    if(psmethod!=2){
      ggg1=MLfunct(y=dd[trsample], x=x1[trsample,], MLmethod=MLmethod,  ybin=0)
      pred1=predict(ggg1, x1[tesample,], onlySL = TRUE)$pred
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
      pstreat=MLfunct(y=wtreat[trsample], x=x1[trsample,], MLmethod=MLmethod,  ybin=0)
      pstreat=predict(pstreat, x1[tesample,], onlySL = TRUE)$pred
      pscontrol=MLfunct(y=wcontrol[trsample], x=x1[trsample,], MLmethod=MLmethod,  ybin=0)
      pscontrol=predict(pscontrol, x1[tesample,], onlySL = TRUE)$pred
    }
    param=rbind(param, cbind(pstreat,pscontrol, mutreat, mucontrol))
  }  # close cross-fitting loop
  param = param[order(idx),]                              # sort nuisance parameters according to original order of observations
  param=cbind(kernwgtdtreat, kernwgtdcontrol, kernwgtdtreatnorm, kernwgtdcontrolnorm, y, param)
  # columns: 1 kernwgtdtreat, 2 kernwgtdcontrol, 3 kernwgtdtreatnorm, 4 kernwgtdcontrolnorm, 5 y, 6 pstreat, 7 pscontrol, 8 mutreat, 9 mucontrol
  
  if (ATET!=TRUE){
    # ATE: doubly robust kernel/GPS-weighted score at each dose, following Kennedy et al (2017) / Colangelo and Lee (2020)
    wtreatraw=param[,1]/param[,6]                                  # K_h(D-dtreat)/pscore(dtreat|X)
    wcontrolraw=param[,2]/param[,7]                                # K_h(D-dcontrol)/pscore(dcontrol|X)
    wtreatsum=wtreatraw/sum(wtreatraw); wcontrolsum=wcontrolraw/sum(wcontrolraw)   # for trimming
    trimmed=(wtreatsum>trim | wcontrolsum>trim)
    wtreat=wtreatraw/mean(wtreatraw[trimmed==0]); wcontrol=wcontrolraw/mean(wcontrolraw[trimmed==0]) # sample-mean normalized (using the untrimmed sample)
    phitreat=param[,8] + wtreat*(param[,5]-param[,8])
    phicontrol=param[,9] + wcontrol*(param[,5]-param[,9])
    phitreat=phitreat[trimmed==0]; phicontrol=phicontrol[trimmed==0]  # use untrimmed observations
    effect=mean(phitreat)-mean(phicontrol)
    scoremod=phitreat-phicontrol-effect
    if (is.null(cluster)){
      se=summary(lm(scoremod~1))$coefficients[, "Std. Error"]
    }else{
      se=sqrt(vcovCL(lm(scoremod~1), cluster = cluster[trimmed==0]))[1,1]}
  }
  
  if (ATET==TRUE){
    # ATET: reweight control-dose observations to represent the covariate distribution of units near dtreat,
    # matching didcontDMLpanel's formula (applied here to outcome levels rather than pre-post differences)
    resd1=param[,3]                                                 # kernwgtdtreatnorm
    resd0=(param[,4]*param[,6]/param[,7])/sum(param[,4]*param[,6]/param[,7])   # kernwgtdcontrolnorm * pstreat/pscontrol, renormalized
    trimmed=(resd1>trim | resd0>trim)
    paramtr=param[trimmed==0,]
    resd1tr=resd1[trimmed==0]; resd0tr=resd0[trimmed==0]
    reg=resd1tr*(paramtr[,5]-paramtr[,9])/sum(resd1tr)              # renormalize weights and multiply by (y-mucontrol)
    resd0tr=resd0tr*(paramtr[,5]-paramtr[,9])/sum(resd0tr)
    score=sum(1-trimmed)*(reg-resd0tr)
    effect=mean(score)
    kfunct=param[,1]/bw
    meank=mean(kfunct)
    fd=sum(kfunct)/length(kfunct)
    scoremod=score-effect/fd*(kfunct[trimmed==0]-meank)
    if (is.null(cluster)){
      se=summary(lm(scoremod~1))$coefficients[, "Std. Error"]
    }else{
      se=sqrt(vcovCL(lm(scoremod~1), cluster = cluster[trimmed==0]))[1,1]}
  }
  pval= 2*pnorm((-1)*abs(effect/se))
  list(effect=effect, se=se, pval=pval, ntrimmed=sum(trimmed))
}