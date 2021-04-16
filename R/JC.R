#' Job Corps data
#'
#' A dataset from the U.S. Job Corps experimental study with information on the participation of disadvantaged youths in (academic and vocational) training in the first and second year after program assignment.
#'
#' @format A data frame with 9240 rows and 46 variables:
#' \describe{
#'   \item{assignment}{1=randomly assigned to Job Corps, 0=randomized out of Job Corps}
#'   \item{female}{1=female, 0=male}
#'   \item{age}{age in years at assignment}
#'   \item{white}{1=white, 0=non-white}
#'   \item{black}{1=black, 0=non-black}
#'   \item{hispanic}{1=hispanic, 0=non-hispanic}
#'   \item{educ}{years of education at assignment}
#'   \item{educmis}{1=education missing at assignment}
#'   \item{geddegree}{1=has a GED degree at assignment}
#'   \item{hsdegree}{1=has a high school degree at assignment}
#'   \item{english}{1=English mother tongue}
#'   \item{cohabmarried}{1=cohabiting or married at assignment}
#'   \item{haschild}{1=has at least one child, 0=no children at assignment}
#'   \item{everwkd}{1=has ever worked at assignment, 0=has never worked at assignment}
#'   \item{mwearn}{average weekly gross earnings at assignment}
#'   \item{hhsize}{household size at assignment}
#'   \item{hhsizemis}{1=household size missing}
#'   \item{educmum}{mother's years of education at assignment}
#'   \item{educmummis}{1=mother's years of education missing}
#'   \item{educdad}{father's years of education at assignment}
#'   \item{educdadmis}{1=father's years of education missing}
#'   \item{welfarechild}{welfare receipt during childhood in categories from 1 to 4 (measured at assignment)}
#'   \item{welfarechildmis}{1=missing welfare receipt during childhood}
#'   \item{health}{general health at assignment from 1 (excellent) to 4 (poor)}
#'   \item{healthmis}{1=missing health at assignment}
#'   \item{smoke}{extent of smoking at assignment in categories from 1 to 4}
#'   \item{smokemis}{1=extent of smoking missing}
#'   \item{alcohol}{extent of alcohol consumption at assignment in categories from 1 to 4}
#'   \item{alcoholmis}{1=extent of alcohol consumption missing}
#'   \item{everwkdy1}{1=has ever worked one year after assignment, 0=has never worked one year after assignment}
#'   \item{earnq4}{weekly earnings in fourth quarter after assignment}
#'   \item{earnq4mis}{1=missing weekly earnings in fourth quarter after assignment}
#'   \item{pworky1}{proportion of weeks employed in first year after assignment}
#'   \item{pworky1mis}{1=missing proportion of weeks employed in first year after assignment}
#'   \item{health12}{general health 12 months after assignment from 1 (excellent) to 4 (poor)}
#'   \item{health12mis}{1=missing general health 12 months after assignment}
#'   \item{trainy1}{1=enrolled in education and/or vocational training in the first year after assignment, 0=no education or training in the first year after assignment}
#'   \item{trainy2}{1=enrolled in education and/or vocational training in the second year after assignment, 0=no education or training in the second year after assignment}
#'   \item{pworky2}{proportion of weeks employed in second year after assignment}
#'   \item{pworky3}{proportion of weeks employed in third year after assignment}
#'   \item{pworky4}{proportion of weeks employed in fourth year after assignment}
#'   \item{earny2}{weekly earnings in second year after assignment}
#'   \item{earny3}{weekly earnings in third year after assignment}
#'   \item{earny4}{weekly earnings in fourth year after assignment}
#'   \item{health30}{general health 30 months after assignment from 1 (excellent) to 4 (poor)}
#'   \item{health48}{general health 48 months after assignment from 1 (excellent) to 4 (poor)}
#' }
#' @docType data
#' @references Schochet, P. Z., Burghardt, J., Glazerman, S. (2001): "National Job Corps study: The impacts of Job Corps on participants' employment and related outcomes", Mathematica Policy Research, Washington, DC.
#' @examples
#' \dontrun{
#' data(JC)
#' # Dynamic treatment effect evaluation of training in 1st and 2nd year
#' # define covariates at assignment (x0) and after one year (x1)
#' x0=JC[,2:29]; x1=JC[,30:36]
#' # define treatment (training) in first year (d1) and second year (d2)
#' d1=JC[,37]; d2=JC[,38]
#' # define outcome (weekly earnings in fourth year after assignment)
#' y2=JC[,44]
#' # assess dynamic treatment effects (training in 1st+2nd year vs. no training)
#' output=dyntreatDML(y2=y2, d1=d1, d2=d2, x0=x0, x1=x1)
#' cat("dynamic ATE: ",round(c(output$effect),3),", standard error: ",
#'     round(c(output$se),3), ", p-value: ",round(c(output$pval),3))}
"JC"

