#' Causal mediation analysis based on inverse probability weighting
#' @description Causal mediation analysis (evaluation of natural direct and indirect effects) based on weighting by the inverse of treatment propensity scores as suggested in Huber (2014).
#' @param y Dependent variable, must not contain missings.
#' @param d Treatment, must be binary (either 1 or 0), must not contain missings.
#' @param m Mediator(s), may be a scalar or a vector, must not contain missings.
#' @param x Pre-treatment confounders of the treatment, mediator, and/or outcome, must not contain missings.
#' @param w Post-treatment confounders of the mediator and the outcome. Default is NULL. Must not contain missings.
#' @param ATET If FALSE, the average treatment effect (ATE) and the corresponding direct and indirect effects are estimated. If TRUE, the average treatment effect on the treated (ATET)  and the corresponding direct and indirect effects are estimated. Default is FALSE.
#' @param trim Trimming rule for discarding observations with extreme propensity scores. In the absence of post-treatment confounders (w=NULL), observations with Pr(D=1|M,X)<\code{trim} or Pr(D=1|M,X)>(1-\code{trim}) are dropped.
#'  In the presence of post-treatment confounders (\code{w} is defined), observations with Pr(D=1|M,W,X)<\code{trim} or Pr(D=1|M,W,X)>(1-\code{trim}) are dropped. Default is 0.05.
#' @param logit If FALSE, probit regression is used for propensity score estimation. If TRUE, logit regression is used. Default is FALSE.
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @details Estimation of causal mechanisms (natural direct and indirect effects) of a binary treatment under a selection on observables assumption assuming that all confounders of the treatment and the mediator, the treatment and the outcome, or the mediator and the outcome are observed. Units are weighted by the inverse of their conditional treatment propensities given the mediator and/or observed confounders, which are estimated by probit or logit regression.
#' The form of weighting depends on whether the observed confounders are exclusively pre-treatment (\code{x}), or also contain post-treatment confounders of the mediator and the outcome (\code{w}). In the latter case, only partial indirect effects (from  \code{d} to \code{m} to \code{y}) can be estimated that exclude any causal paths from \code{d} to \code{w} to \code{m} to \code{y}, see the discussion in Huber (2014). Standard errors are obtained by bootstrapping the effects.
#' @return A medweight object contains two components, \code{results} and \code{ntrimmed}:
#' @return \code{results}: a 3X5 matrix containing the effect estimates in the first row ("effects"), standard errors in the second row ("se"), and p-values in the third row ("p-value").
#' The first column provides the total effect, namely the average treatment effect (ATE) if \code{ATET=FALSE} or the average treatment effect on the treated (ATET) if \code{ATET=TRUE}.
#' The second and third columns provide the direct effects under treatment and control, respectively ("dir.treat", "dir.control"). See equation (6) if \code{w=NULL} (no post-treatment confounders) and equation (13) if \code{w} is defined, respectively, in Huber (2014). If \code{w=NULL}, the fourth and fifth columns provide the indirect effects under treatment and control, respectively ("indir.treat", "indir.control"), see equation (7) in Huber (2014).
#' If w is defined, the fourth and fifth columns provide the partial indirect effects under treatment and control, respectively ("par.in.treat", "par.in.control"), see equation (14) in Huber (2014).
#' @return \code{ntrimmed}: number of discarded (trimmed) observations due to extreme propensity score values.
#' @references Huber, M. (2014): "Identifying causal mechanisms (primarily) based on inverse probability weighting",  Journal of Applied Econometrics, 29, 920-943.
#' @examples # A little example with simulated data (10000 observations)
#' n=10000
#' x=rnorm(n)
#' d=(0.25*x+rnorm(n)>0)*1
#' w=0.2*d+0.25*x+rnorm(n)
#' m=0.5*w+0.5*d+0.25*x+rnorm(n)
#' y=0.5*d+m+w+0.25*x+rnorm(n)
#' # The true direct and partial indirect effects are all equal to 0.5
#' output=medweight(y=y,d=d,m=m,x=x, w=w, trim=0.05, ATET=FALSE, logit=TRUE, boot=19)
#' round(output$results,3)
#' output$ntrimmed
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm
#' @import mvtnorm
#' @export
medweight<-function(y,d,m,x, w=NULL, ATET=FALSE, trim=0.05, logit=FALSE, boot=1999){
  temp=mediation(y=y,d=d,m=m,x=x,w=w,trim=trim, ATET=ATET, logit=logit)
  ntrimmed=temp[length(temp)]
  temp=temp[1:(length(temp)-1)]
  temp2=bootstrap.mediation(y=y,d=d,m=m,x=x,w=w,boot=boot,trim=trim, ATET=ATET, logit=logit)
  se=apply(temp2[,1:(ncol(temp2)-1)], 2, sd)
  temp3=2*pnorm(-abs(temp/se))
  results=rbind(temp, se, temp3)
  if (is.null(w)==TRUE & ATET==FALSE) colnames(results)=c("ATE", "dir.treat", "dir.control", "indir.treat", "indir.control")
  if (is.null(w)==FALSE & ATET==FALSE ) colnames(results)=c("ATE", "dir.treat", "dir.control", "par.in.treat", "par.in.control")
  if (is.null(w)==TRUE & ATET==TRUE) colnames(results)=c("ATET", "dir.treat", "dir.control", "indir.treat", "indir.control")
  if (is.null(w)==FALSE & ATET==TRUE) colnames(results)=c("ATET", "dir.treat", "dir.control", "par.in.treat", "par.in.control")
  rownames(results)=c("effect", "se", "p-value")
  list(results=results, ntrimmed=ntrimmed)
}



