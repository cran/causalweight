#' Difference-in-Differences for Mediation Analysis with Panel Data and Discrete Treatments Using Double Machine Learning
#' @description This function estimates the total effect, natural direct effect, and natural indirect effect for the treated group in the post-treatment period in panel data with discrete treatments. Estimation is based on a difference-in-differences approach for mediation analysis combined with double machine learning to control for (possibly time-varying) confounders in a data-driven manner. The function supports various machine learning methods for estimating nuisance parameters through k-fold cross-fitting.
#' @param y0 Outcome variable in the pre-treatment period. Should not contain missing values.
#' @param y1 Outcome variable in the post-treatment period. Should not contain missing values.
#' @param d Treatment group indicator (discrete). Should not contain missing values.
#' @param m Mediator variable. Should not contain missing values.
#' @param x Covariates to be controlled for. Should not contain missing values.
#' @param dtreat Value of the treatment under treatment (in the treatment period of interest). Default is 1.
#' @param dcontrol Value of the treatment under control (in the treatment period of interest). Default is 0.
#' @param MLmethod Machine learning method for estimating nuisance parameters using the \code{SuperLearner} package. Must be one of \code{"lasso"} (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"}, \code{"ensemble"}, or \code{"parametric"}.
#' @param trim Trimming threshold for discarding observations with too small propensity scores in the control group. Default is 0.05.
#' @param cluster Optional clustering variable for calculating cluster-robust standard errors.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @details This function estimates the total effect, natural direct effect, and natural indirect effect for the treated group in the post-treatment period in panel data with discrete treatments. Estimation is based on the difference-in-differences approach for mediation analysis proposed by Huber and Oberhänsli (2026). Specifically, these effects are computed from normalized sample analogs of the doubly robust expressions in equations (27) and (28). Double machine learning is used to control for confounders in a data-adaptive way. The function supports different machine learning methods for estimating nuisance parameters (conditional mean outcomes and propensity scores) as well as cross-fitting to mitigate overfitting.
#' @return A list with the following components:
#' @return \code{eff}: Estimates of the total effect, natural direct effect, and natural indirect effect for the treated group in the post-treatment period.
#' @return \code{se}: Standard error of the estimates.
#' @return \code{tval}: t-value of the estimates.
#' @return \code{pval}: p-value of the estimates.
#' @return \code{trimmed}: Number of discarded (trimmed) observations.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @references Huber, M., and Oberhänsli, S. J. (2026): "Difference-in-differences for mediation analysis using double machine learning", arXiv preprint 2602.23877.
#' @examples
#' \dontrun{
#' # Example with simulated data
#' n=4000                            # sample size
#' u=rnorm(n)                        # time constant unobservable
#' x=rnorm(n)                        # covariate
#' d=1*(x+0.5*u+rnorm(n)>0)          # treatment
#' m=x+0.5*d+rnorm(n)                # mediator
#' y0=u+rnorm(n)                     # outcome in the pre-treatment period
#' y1=x+1+d+m+u+rnorm(n)             # outcome in the post-treatment period
#' # true NDET is equal to 1; true NIET is equal to 0.5; true ATET is equal to 1.5
#' didmedDMLpanel(y0=y0,y1=y1, d=d, m=m, x=x)
#' }
#' @importFrom stats rnorm lm predict sd dnorm
#' @import clubSandwich
#' @export
didmedDMLpanel<-function(y0,y1, d, m, x, dtreat=1, dcontrol=0, MLmethod="lasso", trim=0.05, cluster=NULL, k=3){
  ydiff=y1-y0
  ybin=1*(length(unique(ydiff))==2 & min(ydiff)==0 & max(ydiff)==1)
  controls=data.frame(x)
  dt=1*(d==dtreat); dc=1*(d==dcontrol)
  set.seed(1); idx = sample(length(d), replace=FALSE)   # shuffle data
  folds = split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))  # generate folds
  param=c();
  for (i in 1:k){                                                         # start of cross-fitting loop
    tesample=folds[[i]]
    trsample=idx[!(idx %in% tesample)]                                                                # cross-fitting loop
    ydifftr=ydiff[trsample]; dttr=dt[trsample]; dctr=dc[trsample]
    Xtr=data.frame(1,controls)[trsample,] # covariates in the training data
    Xte=data.frame(1,controls)[tesample,] # covariates in the test data

    # plug-in nuisance parameters for E[Y1(d',M(d'))|D=d]
    # E(Y1-Y0|D=d',X)
    mudc=MLfunct(y=ydifftr, x=Xtr, d1=1*(dctr==1), MLmethod=MLmethod, ybin=ybin)
    mudc=predict(mudc, Xte, onlySL = TRUE)$pred
    # Pr(D=d|X)
    pidt=MLfunct(y=1*(dttr==1), x=Xtr, MLmethod=MLmethod,  ybin=1)
    pidt=predict(pidt, Xte, onlySL = TRUE)$pred
    # Pr(D=d'|X)
    pidc=MLfunct(y=1*(dctr==1), x=Xtr, MLmethod=MLmethod,  ybin=1)
    pidc=predict(pidc, Xte, onlySL = TRUE)$pred

    Xmtr=data.frame(1,m,controls)[trsample,] # covariates and mediator in the training data
    Xmte=data.frame(1,m,controls)[tesample,] # covariates and mediator in the test data

    # plug-in nuisance parameters for E[Y1(d',M(d))|D=d]
    # E(Y1-Y0|D=d,M,X)
    mudcm=MLfunct(y=ydifftr, x=Xmtr, d1=1*(dctr==1), MLmethod=MLmethod, ybin=ybin)
    mudcm=predict(mudcm, Xmte, onlySL = TRUE)$pred
    # Pr(D=d|M,X)
    pidtm=MLfunct(y=1*(dttr==1), x=Xmtr, MLmethod=MLmethod,  ybin=1)
    pidtm=predict(pidtm, Xmte, onlySL = TRUE)$pred
    # Pr(D=d'|M,X)
    pidcm=MLfunct(y=1*(dctr==1), x=Xmtr, MLmethod=MLmethod,  ybin=1)
    pidcm=predict(pidcm, Xmte, onlySL = TRUE)$pred

    vars=c("pidt","pidc","mudc","pidtm","pidcm","mudcm")
    param <- rbind(param, as.data.frame(mget(vars)))
  } # close cross-fitting loop
  param=param[order(idx),]
  param=cbind(dt, dc, m, y1, y0, ydiff, param)

  # identify observations with small propensity scores in subgroups which are dropped
  trimmed=1*((param$pidc<trim*param$dc) | (param$pidcm<trim*param$dc))
  param=param[trimmed==0,]
  param$pidt=ifelse(param$pidt == 0, .Machine$double.eps, param$pidt)
  param$pidc=ifelse(param$pidc == 0, .Machine$double.eps, param$pidc)
  param$pidtm=ifelse(param$pidtm == 0, .Machine$double.eps, param$pidtm)
  param$pidcm=ifelse(param$pidcm == 0, .Machine$double.eps, param$pidcm)

  # compute E[Y1(d',M(d'))|D=d]
  basecdt=param$dt*(param$y0+param$mudc)/sum(param$dt)
  rescdc=(param$ydiff-param$mudc)*(param$dc*param$pidt/param$pidc)/sum(param$dc*param$pidt/param$pidc)
  ydcmdc=sum(1-trimmed)*(basecdt+rescdc)

  # compute E[Y1(d,M(d0))|D=d]
  basetdt=param$dt*(param$y0+param$mudcm)/sum(param$dt)
  restdc=(param$ydiff-param$mudcm)*(param$dc*param$pidtm/param$pidcm)/sum(param$dc*param$pidtm/param$pidcm)
  ydcmdt=sum(1-trimmed)*(basetdt+restdc)

  # compute E[Y1(d,M(d))|D=d]
  ydtmdt=sum(1-trimmed)*(param$dt*param$y1/sum(param$dt))

  # compute total effect, natural direct effect and natural indirect effect
  score_te=ydtmdt-ydcmdc
  te=mean(score_te)
  score_nde=ydtmdt-ydcmdt
  nde=mean(score_nde)
  score_nie=ydcmdt-ydcmdc
  nie=mean(score_nie)
  if (is.null(cluster)){
    se_te=summary(lm(score_te~1))$coefficients[, "Std. Error"] # se without clustering
    se_nde=summary(lm(score_nde~1))$coefficients[, "Std. Error"]
    se_nie=summary(lm(score_nie~1))$coefficients[, "Std. Error"]
  }else{
    se_te=coef_test((lm(score_te~1)), vcov = "CR2", cluster = cluster[trimmed == 0])$SE # se with clustering
    se_nde=coef_test((lm(score_nde~1)), vcov = "CR2", cluster = cluster[trimmed == 0])$SE
    se_nie=coef_test((lm(score_nie~1)), vcov = "CR2", cluster = cluster[trimmed == 0])$SE}
  #p-values
  pval_te=2*pnorm((-1)*abs(te/se_te))
  pval_nde=2*pnorm((-1)*abs(nde/se_nde))
  pval_nie=2*pnorm((-1)*abs(nie/se_nie))

  eff=rbind(te,nde,nie)
  se=rbind(se_te,se_nde,se_nie)
  tval=rbind(te/se_te, nde/se_nde, nie/se_nie)
  pval=rbind(pval_te,pval_nde,pval_nie)
  trimmed_obs=sum(trimmed)

  output <- cbind(eff, se, tval, pval)
  rownames(output) <- c("Total Effect",
                        "Natural Direct Effect",
                        "Natural Indirect Effect")
  colnames(output) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

  cat("\nDifference-in-Differences for Mediation Analysis:\n\n")
  printCoefmat(output, P.values = TRUE, has.Pvalue = TRUE, digits = 4,
               signif.stars = TRUE, eps.Pvalue = 2e-16)
  cat("Number of trimmed observations:", trimmed_obs, "\n")
  cat("\n",
      paste0("Value of the treatment under treatment: D = ", dtreat),
      "\n",
      paste0("Value of the treatment under control: D =  ", dcontrol),
      sep = "",
      "\n\n")

  invisible(
    list(
      eff = eff,
      se = se,
      tval = tval,
      pval = pval,
      trimmed_obs = trimmed_obs
    )
  )
}
