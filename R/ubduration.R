#' Austrian unemployment duration data
#'
#' A dataset containing unemployed females between 46 and 53 years old living in an Austrian region where an extension of the maximum duration of unemployment benefits (from 30 to 209 weeks under particular conditions) for job seekers aged 50 or older was introduced.
#'
#' @format A data frame with 5659 rows and 10 variables:
#' \describe{
#'   \item{y}{Outcome variable: unemployment duration of the jobseeker in weeks (registered at the unemployment office). Variable is numeric.}
#'   \item{z}{Running variable: distance to the age threshold of 50 (implying an extended duration of unemployment benefits), measured in months divided by 12. Variable is numeric.}
#'   \item{marrstatus}{Marital status: 0=other, 1=married, 2=single. Variable is a factor.}
#'   \item{education}{Eductation: 0=low education, 1=medium education, 2=high education. Variable is ordered.}
#'   \item{foreign}{Migrant status: 1=foreigner, 0=Austrian. Variable is a factor.}
#'   \item{rr}{Replacement rate (of previous earnings by unemployment benefits). Variable is numeric.}
#'   \item{lwageljob}{Log wage in last job. Variable is numeric.}
#'   \item{experience}{Ratio of actual to potential work experience. Variable is numeric.}
#'   \item{whitecollar}{1=white collar worker, 0=blue collar worker. Variable is a factor.}
#'   \item{industry}{Industry: 0=other, 1=agriculture, 2=utilities, 3=food, 4=textiles, 5=wood, 6=machines, 7=other manufacturing, 8=construction, 9=tourism, 10=traffic, 11=services. Variable is a factor.}
#' }
#' @docType data
#' @references Lalive, R. (2008): "How Do Extended Benefits Affect Unemployment Duration? A Regression Discontinuity Approach", Journal of Econometrics, 142, 785–806.
#' @references Frölich, M. and Huber, M. (2019): "Including covariates in the regression discontinuity design", Journal of Business & Economic Statistics, 37, 736-748.
#' @examples
#' \dontrun{
#' # load unemployment duration data
#' data(ubduration)
#' # run sharp RDD conditional on covariates with user-defined bandwidths
#' RDDcovar(y=ubduration[,1],z=ubduration[,2],x=ubduration[,c(-1,-2)],
#'  bw0=c(0.17, 1, 0.01, 0.05, 0.54, 70000, 0.12, 0.91, 100000),
#'  bw1=c(0.59, 0.65, 0.30, 0.06, 0.81, 0.04, 0.12, 0.76, 1.03),bwz=0.2,boot=19)
#' cat("RDD effect estimate: ",round(c(output$effect),3),", standard error: ",
#'  round(c(output$se),3), ", p-value: ", round(c(output$pvalue),3))}
"ubduration"
