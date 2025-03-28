% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATETDML.R
\name{ATETDML}
\alias{ATETDML}
\title{ATET Estimation for Binary Treatments using Double Machine Learning}
\usage{
ATETDML(
  y,
  d,
  x,
  MLmethod = "lasso",
  est = "dr",
  trim = 0.05,
  cluster = NULL,
  k = 3
)
}
\arguments{
\item{y}{Outcome variable. Should not contain missing values.}

\item{d}{Treatment group indicator (binary). Should not contain missing values.}

\item{x}{Covariates to be controlled for. Should not contain missing values.}

\item{MLmethod}{Machine learning method for estimating nuisance parameters using the \code{SuperLearner} package. Must be one of \code{"lasso"} (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"}, \code{"ensemble"}, or \code{"parametric"}.}

\item{est}{Estimation method. Must be one of \code{"dr"} (default) for doubly robust, \code{"ipw"} for inverse probability weighting (not doubly robust!), or \code{"reg"} for regression (not doubly robust!).}

\item{trim}{Trimming threshold (in percentage) for discarding observations with propensity scores too close to one. Default is 0.05, implying that treatment propensity scores larger than 1-0.05=0.95 (such that the probability to be not treated is below 0.05) are trimmed.}

\item{cluster}{Optional clustering variable for calculating cluster-robust standard errors.}

\item{k}{Number of folds in k-fold cross-fitting. Default is 3.}
}
\value{
A list with the following components:

\code{ATET}: Estimate of the Average Treatment Effect on the Treated (ATET) in the post-treatment period.

\code{se}: Standard error of the ATET estimate.

\code{pval}: P-value of the ATET estimate.

\code{trimmed}: Number of discarded (trimmed) observations.

\code{treat}: Treatments of untrimmed observations.

\code{outcome}: Outcomes of untrimmed observations.

\code{pscores}: Treatment propensity scores of untrimmed observations.

\code{outcomepred}: Conditional outcome predictions under nontreatment of untrimmed observations.
}
\description{
This function estimates the average treatment effect on the treated (ATET) for a binary treatment. Combines estimation based on (doubly robust) efficient score functions with double machine learning to control for confounders in a data-driven way.
}
\details{
This function estimates the Average Treatment Effect on the Treated (ATET) under conditional independence, assuming that confounders jointly affecting the treatment and the outcome can be controlled for by observed covariates. Estimation is based on the (doubly robust) efficient score functions for potential outcomes in combination with double machine learning with cross-fitting, see Chernozhukov et al (2018). To this end, one part of the data is used for estimating the model parameters of the treatment and outcome equations based machine learning. The other part of the data is used for predicting the efficient score functions. The roles of the data parts are swapped (using k-fold cross-fitting) and the ATET is estimated based on averaging the predicted efficient score functions in the total sample. Besides double machine learning, the function also provides inverse probability weighting and regression adjustment methods (which are, however, not doubly robust).
}
\examples{
\dontrun{
n=2000                            # sample size
p=50                              # number of covariates
s=2                               # number of covariates that are confounders
x=matrix(rnorm(n*p),ncol=p)       # covariate matrix
beta=c(rep(0.25,s), rep(0,p-s))   # coefficients determining degree of confounding
d=(x\%*\%beta+rnorm(n)>0)*1         # treatment equation
y=x\%*\%beta+0.5*d+rnorm(n)         # outcome equation
# The true ATE is equal to 0.5
output=ATETDML(y,d,x)
cat("ATET: ",round(c(output$ATET),3),", standard error: ",
    round(c(output$se),3), ", p-value: ",round(c(output$pval),3))
output$ntrimmed
}
}
\references{
Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
}
