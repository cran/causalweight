#' paneltestDML: Overidentification test for ATET estimation in panel data 
#' @description This function applies an overidentification test to assess if unconfoundedness of the treatments and conditional common trends as imposed in differences-in-differences jointly hold in panel data when evaluating the average treatment effect on the treated (ATET). 
#' @param y1 Potential outcomes for the treated.
#' @param y0 Potential outcomes for the non-treated.
#' @param d Treatment group indicator (binary). Should not contain missing values.
#' @param x Covariates to be controlled for. Should not contain missing values.
#' @param MLmethod Machine learning method for estimating nuisance parameters using the \code{SuperLearner} package. Must be one of \code{"lasso"} (default), \code{"randomforest"}, \code{"xgboost"}, \code{"svm"}, \code{"ensemble"}, or \code{"parametric"}.
#' @param trim Trimming threshold for discarding observations with too small propensity scores within any subgroup defined by the treatment group and time. Default is 0.05.
#' @param k Number of folds in k-fold cross-fitting. Default is 3.
#' @details The test statistic corresponds to the difference between the ATETs that are based on two distinct doubly robust score functions, namely that under unconfoundedness and that based on difference-in-differences under conditional common trends. Estimation in panel data is based on double machine learning and the function supports different machine learning methods to estimate nuisance parameters (conditional mean outcomes and propensity scores) as well as cross-fitting to mitigate overfitting.
#' ATETselobs and ATETdid equals zero.
#' @return A list with the following components:
#' \item{est}{Test statistic.}
#' \item{se}{Standard error.}
#' \item{pval}{P-value.}
#' \item{ntrimmed}{Number of trimmed or dropped observations due to propensitiy scores below the threshold \code{trim}.}
#' \item{pscore.xy0}{Propensity score under unconfoundedness.}
#' \item{pscore.x}{Propensity score under conditional common trends.}
#' \item{ATETselobs}{ATET based on the selection on observables/unconfoundedness assumption.}
#' \item{seATETselobs}{Standard error of the ATET based on the selection on observables/unconfoundedness assumption.}
#' \item{ATETdid}{ATET based on difference-in-differences invoking the conditional common trends assumption.}
#' \item{seATETdid}{Standard error of the ATET based on difference-in-differences invoking the conditional common trends assumption.}
#' @import SuperLearner
#' @import sandwich
#' @export
#' @references Huber, M., and Oeß, E.-M. (2024): "A joint test of unconfoundedness and common trends", arXiv preprint 2404.16961.
#' @examples
#' \dontrun{
#' n=1000
#' x=matrix(rnorm(n * 5), n, 5)
#' d=1*(x[,1]+2*rnorm(n)>0)
#' t=rbinom(n, 1, 0.5)
#' y0=x[,1]+rnorm(n)
#' y1=y0+rnorm(n)
#' y=ifelse(t == 1, y1, y0)
#' # report p-value (note that unconfoundedness and common trends hold jointly)
#' paneltestDML(y1, y0, d, x)$pval}
paneltestDML=function(y1, y0, d, x, trim = 0.01, MLmethod = "lasso", k = 3){
  
  # Call hdtest to Estimate nuisance parameters
  nuisance = hdtest(y1 = y1, y0 = y0, d = d, x = x, trim = 0.01, MLmethod = MLmethod, k = k)
  
  # Trim observations with extreme propensity scores close to one
  trimmed = nuisance[, 8]            
  nuisance = nuisance[trimmed == 0,]  
  
  # Estimate sample weights
  sumtreated = sum(nuisance[, 1])
  sumweightp = sum((1 - nuisance[, 1]) * nuisance[, 4] / (1 - nuisance[, 4]))
  sumweightpi = sum((1 - nuisance[, 1]) * nuisance[, 5] / (1 - nuisance[, 5]))
  
  
  # Observational level score for the test-statistic
  score = (nuisance[, 1] / sumtreated) * (nuisance[, 6] - nuisance[, 3] - nuisance[, 7]) +
    ((1 - nuisance[, 1]) / sumweightp) * (nuisance[, 4] * (nuisance[, 2] - nuisance[, 6]) / (1 - nuisance[, 4])) -
    ((1 - nuisance[, 1]) / sumweightpi) * (nuisance[, 5] * (nuisance[, 2] - nuisance[, 3] - nuisance[, 7]) / (1 - nuisance[, 5]))
  
  # Sum score over all observations
  theta = sum(score)
  
  # Number of observations
  n = length(nuisance[, 1])
  
  # Estimate standard error
  se = sqrt(mean((n * score - theta)^2) / n)
  
  # Estimate p-value
  pval = 2 * pnorm(-abs(theta / se))
  
  
  ## score ATET using selection on observables/unconfoundedness assumption (Huber 2023)
  scoreselobs = nuisance[, 1] / sumtreated * (nuisance[, 2] - nuisance[, 6]) -
    (1 - nuisance[, 1]) / sumweightp * (nuisance[, 4] * (nuisance[, 2] - nuisance[, 6]) / (1 - nuisance[, 4]))
  
  # ATETselobs point estimate
  ATETselobs = sum(scoreselobs)
  
  # ATETselobs variance
  varATETselobs = mean(scoreselobs^2)
  
  ## score ATET using conditional common trends (Zhao Sant'Anna 2019)
  scoredid = nuisance[, 1] / sumtreated * (nuisance[, 2] - nuisance[, 3] - nuisance[, 7]) -
    (1 - nuisance[, 1]) / sumweightpi * (nuisance[, 5] * (nuisance[, 2] - nuisance[, 3] - nuisance[, 7]) / (1 - nuisance[, 5]))
  
  # ATETdid point estimate
  ATETdid = sum(scoredid)
  
  # ATETdid variance
  varATETdid = mean(scoredid^2)
  
  # Summarise results
  list(
    est = theta, 
    se = se, 
    pval = pval, 
    ntrimmed = sum(trimmed), 
    pscore.xy0 = nuisance[, 4], 
    pscore.x = nuisance[, 5], 
    ATETselobs = ATETselobs, 
    seATETselobs = sqrt(varATETselobs), 
    ATETdid = ATETdid, 
    seATETdid = sqrt(varATETdid)
  )
}