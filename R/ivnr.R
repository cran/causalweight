#' Instrument-based treatment evaluation under endogeneity and non-response bias
#' @description Non- and semiparaemtric treatment effect estimation under treatment endogeneity and selective non-response in the outcome based on a binary instrument for the treatment and a continous instrument for response.
#' @param y Dependent variable.
#' @param d Treatment, must be binary and must not contain missings.
#' @param r Response, must be a binary indicator for whether the outcome is observed.
#' @param z1 Binary instrument for the treatment, must not contain missings.
#' @param z2 Continuous instrument for response, must not contain missings.
#' @param x A data frame of covariates to be included in the nonparametric estimation, must not contain missings. Factors and ordered varaibles must be appropriately defined as such by \code{factor()} and \code{ordered()}. Default is \code{NULL} (no covariates included). Covariates are only considered if both \code{x} and \code{xpar} are not \code{NULL}.
#' @param xpar Covariates  to be included in the semiparametric estimation, must not contain missings. Default is \code{NULL} (no covariates included). Covariates are only considered if both \code{x} and \code{xpar} are not \code{NULL}.
#' @param ruleofthumb  If 1, bandwidth selection in any kernel function is based on the Silverman (1986) rule of thumb. Otherwise, least squares cross-validation is used. Default is 1.
#' @param wgtfct Weighting function to be used in effect estimation. If set to 1, equation (18) in Fricke et al (2020) is used as weight. If set to 2, equation (19) in Fricke et al (2020) is used as weight. If set to 3, the median of LATEs across values of response probabilities \code{numresprob} is used. Default is 2.
#' @param rtype Regression type used for continuous outcomes in the kernel regressions. Either \code{"ll"} for local linear or \code{"lc"} for local constant regression. Default is \code{"ll"}.
#' @param numresprob number of response probabilities at which the effects are evaluated. An equidistant grid is constructed based on the number provided. Default is 20.
#' @param boot Number of bootstrap replications for estimating standard errors of the effects. Default is 499.
#' @param estlate If set to \code{TRUE} the local average treatment effect on compliers (LATE) is estimated, otherwise the average treatment effect (ATE) is estimated. Default is \code{TRUE}.
#' @param trim Trimming rule for too extreme denominators in the weighting functions or inverses of products of conditional treatment probabilities. Values below \code{trim} are set to \code{trim} to avoid values that are too close to zero in any denominator. Default is 0.01.
#' @references Fricke, H., Frölich, M., Huber, M., Lechner, M. (2020): "Endogeneity and non-response bias in treatment evaluation - nonparametric identification of causal effects by instruments", Journal of Applied Econometrics, forthcoming.
#' @details Non- and semiparametric treatment effect estimation under treatment endogeneity and selective non-response in the outcome based on a binary instrument for the treatment and a continuous instrument for response. The effects are estimated both semi-parametrically (using probit and OLS for the estimation of plug-in parameters like conditional probabilities and outcomes) and fully non-parametrically (based on kernel regression for any conditional probability/mean). Besides the instrument-based estimates, results are also presented under a missing-at-random assumption (MAR) when not using the instrument \code{z2} for response (but only \code{z1} for the treatment). See Fricke et al. (2020) for further details.
#' @return A \code{ivnr} object contains one output component:
#' @return \code{output}: The first row provides the effect estimates under non- and semi-parametric estimation using both instruments, see \code{"nonpara (L)ATE IV"} and \code{"semipara (L)ATE IV"} as well as under a missing-at-random assumption for response when using only the first instrument for the treatment, see \code{"nonpara (L)ATE MAR"} and \code{"semipara (L)ATE MAR"}. The second row provides the standard errors based on bootstrapping the effects. The third row provides the p-values based on the t-statistics.
#' @examples # A little example with simulated data (1000 observations)
#' \dontrun{
#' n=1000          # sample size
#' e<-(rmvnorm(n,rep(0,3), matrix(c(1,0.5,0.5,  0.5,1,0.5,  0.5,0.5,1),3,3)))
#' # correlated error term of treatment, response, and outcome equation
#' x=runif(n,-0.5,0.5)           # observed confounder
#' z1<-(-0.25*x+rnorm(n)>0)*1    # binary instrument for treatment
#' z2<- -0.25*x+rnorm(n)         # continuous instrument for selection
#' d<-(z1-0.25*x+e[,1]>0)*1      # treatment equation
#'  y_star<- -0.25*x+d+e[,2]     # latent outcome
#'  r<-(-0.25*x+z2+d+e[,3]>0)*1  # response equation
#'  y=y_star                     # observed outcome
#'  y[r==0]=0                    # nonobserved outcomes are set to zero
#'  # The true treatment effect is 1
#'  ivnr(y=y,d=d,r=r,z1=z1,z2=z2,x=x,xpar=x,numresprob=4,boot=39)}
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm dnorm quantile
#' @import np mvtnorm
#' @export


ivnr<-function(y,d,r,z1,z2, x=NULL, xpar=NULL,  ruleofthumb=1, wgtfct=2, rtype="ll", numresprob=20, boot=499, estlate=TRUE, trim=0.01){
  y[r==0]=0
  out=latenonrespxxfct(y,d,r,z1,z2, x=x, xpar=xpar, bres1=NULL, bres0=NULL,  bwyz1=NULL, bwdz1=NULL, bwyz0=NULL, bwdz0=NULL, bwps=NULL, bwcox1=NULL, bwcox2=NULL, bw1=NULL, bw2=NULL, bw3=NULL, bw4=NULL, bw5=NULL, bw6=NULL, bw7=NULL, bw8=NULL, bw9=NULL, bw10=NULL, bw11=NULL, bw12=NULL, ruleofthumb=ruleofthumb, wgtfct=wgtfct, rtype=rtype, numresprob=numresprob, estlate=estlate, trim=trim)
  if ((is.null(x)==0) & (is.null(xpar)==0)) results2=bootstrap.late.nr(y,d,r,z1,z2, x=x, xpar=xpar, bres1=out$bres1, bres0=out$bres0,  bwyz1=out$bwyz1, bwdz1=out$bwdz1, bwyz0=out$bwyz0, bwdz0=out$bwdz0, bwps=out$bwps, bwcox1=out$bwcox1, bwcox2=out$bwcox2, bw1=out$bw1, bw2=out$bw2, bw3=out$bw3, bw4=out$bw4, bw5=out$bw5, bw6=out$bw6, bw7=out$bw7, bw8=out$bw8, bw9=out$bw9, bw10=out$bw10, bw11=out$bw11, bw12=out$bw12, ruleofthumb=ruleofthumb, wgtfct=wgtfct, rtype=rtype, numresprob=numresprob, boot=boot, estlate=estlate, trim=trim)
  if (is.null(x) | is.null(xpar)) results2=bootstrap.late.nr(y,d,r,z1,z2, x=NULL, xpar=NULL, bres1=NULL, bres0=NULL,  bwyz1=NULL, bwdz1=NULL, bwyz0=NULL, bwdz0=NULL, bwps=NULL, bwcox1=NULL, bwcox2=NULL, bw1=out$bw1, bw2=out$bw2, bw3=out$bw3, bw4=out$bw4, bw5=out$bw5, bw6=out$bw6, bw7=out$bw7, bw8=out$bw8, bw9=out$bw9, bw10=out$bw10, bw11=out$bw11, bw12=out$bw12, ruleofthumb=ruleofthumb, wgtfct=wgtfct, rtype=rtype, numresprob=numresprob, boot=boot, estlate=estlate, trim=trim)
  results=out$results
  output<-matrix(NA,3,length(results))
  output[1,]<-results
  for (i in 1:length(results)){
    output[2,i]<-sd(results2[,i])
    output[3,i]<-2*pnorm(-abs(output[1,i]/output[2,i]))
    #output[4,i]<-c(2*(min( mean(( results2[,i]-output[1,i])<=output[1,i]) , mean((results2[,i]-output[1,i])>output[1,i]) ) ))
    colnames(output)=c("nonpara (L)ATE IV", "semipara (L)ATE IV", "nonpara (L)ATE MAR", "semipara (L)ATE MAR")
    rownames(output)=c("estimate", "se", "p-val")
  }
  list(output=output)
}






