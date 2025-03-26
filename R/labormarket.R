#' Temporary Work Agency (TWA) Assignments and Permanent Employment in Sicily
#' 
#' A dataset containing 1120 individuals from Sicily, an Italian region in which a labor market  
#' intervention through a temporary work agency (TWA) was offered. 
#' It includes 229 individuals assigned to TWA jobs in the first 6 months of 2001 and 891 controls 
#' aged 18-40, who were in the labor force but not in stable employment on 1 January 2001, 
#' and who were not on a TWA assignment during the first half of that year.
#' 
#' @format A data frame with 1120 rows and 32 variables:
#' \describe{
#'   \item{age}{age at time 2}
#'   \item{blu1}{blue-collar at time 1}
#'   \item{child}{number of children at time 1}
#'   \item{ct}{dummy Catania}
#'   \item{dist}{distance from nearest agency}
#'   \item{emp1}{employed at time 1}
#'   \item{hour1}{weekly hours of work at time 1}
#'   \item{hour3}{weekly hours of work at time 3}
#'   \item{id}{identification code}
#'   \item{italy}{italian nationality}
#'   \item{male}{male gender}
#'   \item{me}{dummy Messina}
#'   \item{nofl1}{out of labor force at time 1}
#'   \item{nyu1}{fraction of school-to-work without employment}
#'   \item{out3}{employed at time 3}
#'   \item{pa}{dummy Palermo}
#'   \item{single}{non married}
#'   \item{tp}{dummy Trapani}
#'   \item{train1}{received professional training before treatment}
#'   \item{treat}{temp at time 2 - treatment status}
#'   \item{unemp1}{unemployed at time 1}
#'   \item{wage1}{monthly wage at time 1}
#'   \item{wage2}{monthly wage at time 2}
#'   \item{wage3}{monthly wage at time 3}
#'   \item{yu1}{years without employment at time 1}
#'   \item{white1}{white collar at time 1}
#'   \item{unemp3}{unemployed at time 3}
#'   \item{nofl3}{out of labor force at time 3}
#'   \item{edu0}{low education}
#'   \item{edu1}{medium education}
#'   \item{edu2}{high education}
#'   \item{edu}{education level}
#' }
#' @docType data
#' @references Ichino, A., Mealli, F., and Nannicini, T. (2008): "From Temporary Help Jobs to Permanent Employment: What Can We Learn from Matching Estimators and their Sensitivity?", Journal of Applied Econometrics, 23(3), 305-327.
"labormarket"