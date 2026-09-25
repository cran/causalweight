#' Difference-in-Differences for Mediation Analysis with Repeated Cross-Sections and Discrete Treatments Using Double Machine Learning
#' @description This function estimates the total effect, natural direct effect, and natural indirect effect for the treated group in the post-treatment period in repeated cross-sections with discrete treatments. Estimation is based on a difference-in-differences approach for mediation analysis combined with double machine learning to control for (possibly time-varying) confounders in a data-driven manner. The function supports various machine learning methods for estimating nuisance parameters through k-fold cross-fitting.
#' @param y Outcome variable. Should not contain missing values.
#' @param d Treatment group indicator (discrete). Should not contain missing values.
#' @param t Time period indicator (binary). Should be 1 for post-treatment period and 0 for pre-treatment period. Should not contain missing values.
#' @param m Mediator variable. Should not contain missing values.
#' @param x Covariates to be controlled for. Should not contain missing values.
#' @param dtreat Value of the treatment under treatment (in the treatment period of interest). Default is 1.
#' @param dcontrol Value of the treatment under control (in the treatment period of interest). Default is 0.
#' @param MLmethod Machine learning method for estimating nuisance parameters using the \code{SuperLearner} package. Must be one of \code{"lasso"} (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"}, \code{"ensemble"}, or \code{"parametric"}.
#' @param trim Trimming threshold for discarding observations with too small propensity scores in the treated group in the pre-treatment period and in the control group in either period. Default is 0.05.
#' @param cluster Optional clustering variable for calculating cluster-robust standard errors.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @details This function estimates the total effect, natural direct effect, and natural indirect effect for the treated group in the post-treatment period in repeated cross-sections with discrete treatments. Estimation is based on the difference-in-differences approach for mediation analysis proposed by Huber and Oberhänsli (2026). Specifically, these effects are computed from normalized sample analogs of the doubly robust expressions in equations (18) and (20). Double machine learning is used to control for confounders in a data-adaptive way. The function supports different machine learning methods for estimating nuisance parameters (conditional mean outcomes and propensity scores) as well as cross-fitting to mitigate overfitting.
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
#' t=rbinom(n, 1, 0.5)               # time period
#' u=rnorm(n)                        # time constant unobservable
#' x=0.5*t+rnorm(n)                  # time varying covariate
#' d=1*(x+0.5*u+rnorm(n)>0)          # treatment
#' m=x+0.5*d+rnorm(n)                # mediator
#' y=x+(1+d+m)*t+u+rnorm(n)          # outcome
#' # true NDET is equal to 1; true NIET is equal to 0.5; true ATET is equal to 1.5
#' didmedDML(y=y, d=d, t=t, m=m, x=x)
#' }
#' @importFrom stats rnorm lm predict sd dnorm
#' @import clubSandwich
#' @export
didmedDML<-function(y, d, t, m, x, dtreat=1, dcontrol=0, MLmethod="lasso", trim=0.05, cluster=NULL, k=3){
ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
controls=data.frame(x)
dt=1*(d==dtreat); dc=1*(d==dcontrol)
set.seed(1); idx = sample(length(d), replace=FALSE)   # shuffle data
folds = split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))  # generate folds
param=c();
for (i in 1:k){                                                         # start of cross-fitting loop
  tesample=folds[[i]]
  trsample=idx[!(idx %in% tesample)]
  ytr=y[trsample]; dttr=dt[trsample]; dctr=dc[trsample]; ttr=t[trsample]
  Xtr=data.frame(1,controls)[trsample,] # covariates in the training data
  Xte=data.frame(1,controls)[tesample,] # covariates in the test data

  # plug-in nuisance parameters for E[Y_1(d',M(d'))|D=d]
  # E(Y|D=d',T=1,X)
  mdct1=MLfunct(y=ytr, x=Xtr, d1=1*(dctr==1 & ttr==1), MLmethod=MLmethod, ybin=ybin)
  mdct1=predict(mdct1, Xte, onlySL = TRUE)$pred
  # E(Y|D=d,T=0,X)
  mdtt0=MLfunct(y=ytr, x=Xtr, d1=1*(dttr==1 & ttr==0), MLmethod=MLmethod, ybin=ybin)
  mdtt0=predict(mdtt0, Xte, onlySL = TRUE)$pred
  # E(Y|D=d',T=0,X)
  mdct0=MLfunct(y=ytr, x=Xtr, d1=1*(dctr==1 & ttr==0), MLmethod=MLmethod, ybin=ybin)
  mdct0=predict(mdct0, Xte, onlySL = TRUE)$pred
  # Pr(D=d,T=1|X)
  pidtt1=MLfunct(y=1*(dttr==1 & ttr==1), x=Xtr, MLmethod=MLmethod,  ybin=1)
  pidtt1=predict(pidtt1, Xte, onlySL = TRUE)$pred
  # Pr(D=d,T=0|X)
  pidtt0=MLfunct(y=1*(dttr==1 & ttr==0), x=Xtr, MLmethod=MLmethod,  ybin=1)
  pidtt0=predict(pidtt0, Xte, onlySL = TRUE)$pred
  # Pr(D=d',T=1|X)
  pidct1=MLfunct(y=1*(dctr==1 & ttr==1), x=Xtr, MLmethod=MLmethod,  ybin=1)
  pidct1=predict(pidct1, Xte, onlySL = TRUE)$pred
  # Pr(D=d',T=0|X)
  pidct0=MLfunct(y=1*(dctr==1 & ttr==0), x=Xtr, MLmethod=MLmethod,  ybin=1)
  pidct0=predict(pidct0, Xte, onlySL = TRUE)$pred

  Xmtr=data.frame(1,m,controls)[trsample,] # covariates and mediator in the training data
  Xmte=data.frame(1,m,controls)[tesample,] # covariates and mediator in the test data

  # plug-in nuisance parameters for E[Y_1(d',M(d))|D=d]
  # E(Y|D=d',T=1,M,X)
  mudct1=MLfunct(y=ytr, x=Xmtr, d1=1*(dctr==1 & ttr==1), MLmethod=MLmethod, ybin=ybin)
  mudct1=predict(mudct1, Xmte, onlySL = TRUE)$pred
  # E(Y|D=d,T=0,M,X)
  mudtt0=MLfunct(y=ytr, x=Xmtr, d1=1*(dttr==1 & ttr==0), MLmethod=MLmethod, ybin=ybin)
  mudtt0=predict(mudtt0, Xmte, onlySL = TRUE)$pred
  # E(Y|D=d',T=0,M,X)
  mudct0=MLfunct(y=ytr, x=Xmtr, d1=1*(dctr==1 & ttr==0), MLmethod=MLmethod, ybin=ybin)
  mudct0=predict(mudct0, Xmte, onlySL = TRUE)$pred
  # Pr(D=d,T=1|M,X)
  pidtt1m=MLfunct(y=1*(dttr==1 & ttr==1), x=Xmtr, MLmethod=MLmethod,  ybin=1)
  pidtt1m=predict(pidtt1m, Xmte, onlySL = TRUE)$pred
  # Pr(D=d,T=0|M,X)
  pidtt0m=MLfunct(y=1*(dttr==1 & ttr==0), x=Xmtr, MLmethod=MLmethod,  ybin=1)
  pidtt0m=predict(pidtt0m, Xmte, onlySL = TRUE)$pred
  # Pr(D=d',T=1|M,X)
  pidct1m=MLfunct(y=1*(dctr==1 & ttr==1), x=Xmtr, MLmethod=MLmethod,  ybin=1)
  pidct1m=predict(pidct1m, Xmte, onlySL = TRUE)$pred
  # Pr(D=d',T=0|M,X)
  pidct0m=MLfunct(y=1*(dctr==1 & ttr==0), x=Xmtr, MLmethod=MLmethod,  ybin=1)
  pidct0m=predict(pidct0m, Xmte, onlySL = TRUE)$pred

  vars=c("pidtt1","pidtt0","pidct1","pidct0","mdtt0","mdct1","mdct0",
         "pidtt1m","pidtt0m","pidct1m","pidct0m","mudtt0","mudct1","mudct0")
  param <- rbind(param, as.data.frame(mget(vars)))
} # close cross-fitting loop
param=param[order(idx),]
param=cbind(dt, dc, t, m, y, param)
# identify observations with small propensity scores in subgroups which are dropped
trimmed=1*((param$pidtt0<trim*param$dt*(1-param$t)) | (param$pidct1<trim*param$dc*param$t) | (param$pidct0<trim*param$dc*(1-param$t)) |
             (param$pidtt0m<trim*param$dt*(1-param$t)) | (param$pidct1m<trim*param$dc*param$t) | (param$pidct0m<trim*param$dc*(1-param$t)))
param=param[trimmed==0,]
param$pidtt0=ifelse(param$pidtt0 == 0, .Machine$double.eps, param$pidtt0)
param$pidct1=ifelse(param$pidct1 == 0, .Machine$double.eps, param$pidct1)
param$pidct0=ifelse(param$pidct0 == 0, .Machine$double.eps, param$pidct0)
param$pidtt0m=ifelse(param$pidtt0m == 0, .Machine$double.eps, param$pidtt0m)
param$pidct1m=ifelse(param$pidct1m == 0, .Machine$double.eps, param$pidct1m)
param$pidct0m=ifelse(param$pidct0m == 0, .Machine$double.eps, param$pidct0m)

# compute E[Y_1(d',M(d'))|D=d]
basecdtt1=(param$dt*param$t)/sum(param$dt*param$t)*(param$mdtt0+param$mdct1-param$mdct0)
rescdtt0=(param$y-param$mdtt0)*(param$dt*(1-param$t)*param$pidtt1/param$pidtt0)/sum(param$dt*(1-param$t)*param$pidtt1/param$pidtt0)
rescdct1=(param$y-param$mdct1)*(param$dc*param$t*param$pidtt1/param$pidct1)/sum(param$dc*param$t*param$pidtt1/param$pidct1)
rescdct0=(param$y-param$mdct0)*(param$dc*(1-param$t)*param$pidtt1/param$pidct0)/sum(param$dc*(1-param$t)*param$pidtt1/param$pidct0)
ydcmdc=sum(1-trimmed)*(basecdtt1+rescdtt0+rescdct1-rescdct0)

# compute E[Y_1(d',M(d))|D=d]
basetdtt1=(param$dt*param$t)/sum(param$dt*param$t)*(param$mudtt0+param$mudct1-param$mudct0)
restdtt0=(param$y-param$mudtt0)*(param$dt*(1-param$t)*param$pidtt1m/param$pidtt0m)/sum(param$dt*(1-param$t)*param$pidtt1m/param$pidtt0m)
restdct1=(param$y-param$mudct1)*(param$dc*param$t*param$pidtt1m/param$pidct1m)/sum(param$dc*param$t*param$pidtt1m/param$pidct1m)
restdct0=(param$y-param$mudct0)*(param$dc*(1-param$t)*param$pidtt1m/param$pidct0m)/sum(param$dc*(1-param$t)*param$pidtt1m/param$pidct0m)
ydcmdt=sum(1-trimmed)*(basetdtt1+restdtt0+restdct1-restdct0)

# compute E[Y_1(d,M(d))|D=d]
ydtmdt=sum(1-trimmed)*(param$dt*param$t*param$y/sum(param$dt*param$t))

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
    paste0("Value of the treatment under control: D = ", dcontrol),
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
