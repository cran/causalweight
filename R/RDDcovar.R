#' Sharp regression discontinuity design conditional on covariates
#' @description Nonparametric (kernel regression-based) sharp regression discontinuity controlling for covariates that are permitted to jointly affect the treatment assignment and the outcome at the threshold of the running variable, see Frölich and Huber (2019).
#' @param y Dependent variable, must not contain missings.
#' @param z Running variable. Must be coded such that the treatment is zero for \code{z} being smaller than zero and one for \code{z} being larger than or equal to zero. Must not contain missings.
#' @param x Covariates, must not contain missings.
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @param bw0 Bandwidth for a kernel regression of \code{y} on \code{z} and \code{x} below the threshold (for treatment equal to zero), using the Epanechnikov kernel. Default is \code{NULL}, implying that the bandwidth is estimated by least squares cross-validation.
#' @param bw1 Bandwidth for a kernel regression of \code{y} on \code{z} and \code{x} above the threshold (for treatment equal to one), using the Epanechnikov kernel. Default is \code{NULL}, implying that the bandwidth is estimated by least squares cross-validation.
#' @param regtype Defines the type of the kernel regressions of \code{y} on \code{z} and \code{x} below and above the threshold.  Must either be set to \code{"ll"} for local linear regression or to \code{"ll"} for local constant regression. Default is \code{"ll"}.
#' @param bwz Bandwidth for the (Epanechnikov) kernel function on \code{z}. Default is \code{NULL}, implying that the bandwidth is estimated by least squares cross-validation.
#' @details Sharp regression discontinuity design conditional on covariates to control for observed confounders jointly affecting the treatment assignment and outcome at the threshold of the running variable as discussed in Frölich and Huber (2019). This is implemented by running kernel regressions of the outcome on the running variable and the covariates separately above and below the threshold and by applying a kernel smoother to the running variable around the threshold. The procedure permits choosing kernel bandwidths by cross-validation, even though this does in general not yield the optimal bandwidths for treatment effect estimation (checking the robustness of the results by varying the bandwidths is therefore highly recommended). Standard errors are based on bootstrapping.
#' @return \code{effect}: Estimated treatment effect at the threshold.
#' @return \code{se}: Bootstrap-based standard error of the effect estimate.
#' @return \code{pvalue}: P-value based on the t-statistic.
#' @return \code{bw0}: Bandwidth for kernel regression of \code{y} on \code{z} and \code{x} below the threshold (for treatment equal to zero).
#' @return \code{bw1}: Bandwidth for kernel regression of \code{y} on \code{z} and \code{x} above the threshold (for treatment equal to one).
#' @return \code{bwz}: Bandwidth for the kernel function on \code{z}.
#' @references  Frölich, M. and Huber, M. (2019): "Including covariates in the regression discontinuity design", Journal of Business & Economic Statistics, 37, 736-748.
#' @examples
#' \dontrun{
#' # load unemployment duration data
#' data(ubduration)
#' # run sharp RDD conditional on covariates with user-defined bandwidths
#' output=RDDcovar(y=ubduration[,1],z=ubduration[,2],x=ubduration[,c(-1,-2)],
#'  bw0=c(0.17, 1, 0.01, 0.05, 0.54, 70000, 0.12, 0.91, 100000),
#'  bw1=c(0.59, 0.65, 0.30, 0.06, 0.81, 0.04, 0.12, 0.76, 1.03),bwz=0.2,boot=19)
#' cat("RDD effect estimate: ",round(c(output$effect),3),", standard error: ",
#'  round(c(output$se),3), ", p-value: ", round(c(output$pvalue),3))}
#' @import np
#' @importFrom stats var
#' @export


RDDcovar=function(y,z,x, boot=1999, bw0=NULL, bw1=NULL, regtype="ll", bwz=NULL){
  d=1*(z>=0)
  xz=data.frame(x,z)
  xzcutoff=data.frame(x,rep(0,length(d)))
  xz0=xz[d==0,]; xz1=xz[d==1,]; d1=d[d==1]; d0=d[d==0]; y1=y[d==1]; y0=y[d==0];
  if (is.null(bw0)==1) regbw0 <- npregbw(ydat=y0, xdat=xz0, bwmethod="cv.ls", ckertype="epanechnikov", regtype=regtype)$bw
  if (is.null(bw0)==0) regbw0 <- bw0
  if (is.null(bw1)==1) regbw1 <- npregbw(ydat=y1, xdat=xz1, bwmethod="cv.ls", ckertype="epanechnikov", regtype=regtype)$bw
  if (is.null(bw1)==0) regbw1 <- bw1
  if (is.null(bwz)==1) regbw <- npregbw(ydat=y, xdat=z, bwmethod="cv.ls", ckertype="epanechnikov", regtype="lc")$bw
  if (is.null(bwz)==0) regbw <- bwz
  est=rdd.x.est(y=y,z=z,x=x, bw0=bw0, bw1=bw1, regtype=regtype, bwz= bwz)
  se=sqrt(var(rdd.x.boot(y=y,z=z,x=x, bw0=bw0, bw1=bw1, regtype=regtype, bwz=bwz, boot=boot)))
  pvalue=2*pnorm(-abs(est/se))
  list(effect=est, se=se, pvalue=pvalue, bw0=bw0, bw1=bw1, bwz=bwz)
}
