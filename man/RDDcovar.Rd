% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RDDcovar.R
\name{RDDcovar}
\alias{RDDcovar}
\title{Sharp regression discontinuity design conditional on covariates}
\usage{
RDDcovar(
  y,
  z,
  x,
  boot = 1999,
  bw0 = NULL,
  bw1 = NULL,
  regtype = "ll",
  bwz = NULL
)
}
\arguments{
\item{y}{Dependent variable, must not contain missings.}

\item{z}{Running variable. Must be coded such that the treatment is zero for \code{z} being smaller than zero and one for \code{z} being larger than or equal to zero. Must not contain missings.}

\item{x}{Covariates, must not contain missings.}

\item{boot}{Number of bootstrap replications for estimating standard errors. Default is 1999.}

\item{bw0}{Bandwidth for a kernel regression of \code{y} on \code{z} and \code{x} below the threshold (for treatment equal to zero), using the Epanechnikov kernel. Default is \code{NULL}, implying that the bandwidth is estimated by least squares cross-validation.}

\item{bw1}{Bandwidth for a kernel regression of \code{y} on \code{z} and \code{x} above the threshold (for treatment equal to one), using the Epanechnikov kernel. Default is \code{NULL}, implying that the bandwidth is estimated by least squares cross-validation.}

\item{regtype}{Defines the type of the kernel regressions of \code{y} on \code{z} and \code{x} below and above the threshold.  Must either be set to \code{"ll"} for local linear regression or to \code{"ll"} for local constant regression. Default is \code{"ll"}.}

\item{bwz}{Bandwidth for the (Epanechnikov) kernel function on \code{z}. Default is \code{NULL}, implying that the bandwidth is estimated by least squares cross-validation.}
}
\value{
\code{effect}: Estimated treatment effect at the threshold.

\code{se}: Bootstrap-based standard error of the effect estimate.

\code{pvalue}: P-value based on the t-statistic.

\code{bw0}: Bandwidth for kernel regression of \code{y} on \code{z} and \code{x} below the threshold (for treatment equal to zero).

\code{bw1}: Bandwidth for kernel regression of \code{y} on \code{z} and \code{x} above the threshold (for treatment equal to one).

\code{bwz}: Bandwidth for the kernel function on \code{z}.
}
\description{
Nonparametric (kernel regression-based) sharp regression discontinuity controlling for covariates that are permitted to jointly affect the treatment assignment and the outcome at the threshold of the running variable, see Frölich and Huber (2019).
}
\details{
Sharp regression discontinuity design conditional on covariates to control for observed confounders jointly affecting the treatment assignment and outcome at the threshold of the running variable as discussed in Frölich and Huber (2019). This is implemented by running kernel regressions of the outcome on the running variable and the covariates separately above and below the threshold and by applying a kernel smoother to the running variable around the threshold. The procedure permits choosing kernel bandwidths by cross-validation, even though this does in general not yield the optimal bandwidths for treatment effect estimation (checking the robustness of the results by varying the bandwidths is therefore highly recommended). Standard errors are based on bootstrapping.
}
\examples{
\dontrun{
# load unemployment duration data
data(ubduration)
# run sharp RDD conditional on covariates with user-defined bandwidths
output=RDDcovar(y=ubduration[,1],z=ubduration[,2],x=ubduration[,c(-1,-2)],
 bw0=c(0.17, 1, 0.01, 0.05, 0.54, 70000, 0.12, 0.91, 100000),
 bw1=c(0.59, 0.65, 0.30, 0.06, 0.81, 0.04, 0.12, 0.76, 1.03),bwz=0.2,boot=19)
cat("RDD effect estimate: ",round(c(output$effect),3),", standard error: ",
 round(c(output$se),3), ", p-value: ", round(c(output$pvalue),3))}
}
\references{
Frölich, M. and Huber, M. (2019): "Including covariates in the regression discontinuity design", Journal of Business & Economic Statistics, 37, 736-748.
}
