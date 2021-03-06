#' Wage expectations of students in Switzerland
#'
#' A dataset containing information on wage expectations of 804 students
#' at the University of Fribourg and the University of Applied Sciences in Bern in the year 2017.
#'
#' @format A data frame with 804 rows and 39 variables:
#' \describe{
#'   \item{wexpect1}{wage expectations after finishing studies: 0=less than 3500 CHF gross per month; 1=3500-4000 CHF; 2=4000-4500 CHF;...; 15=10500-11000 CHF; 16=more than 11000 CHF}
#'   \item{wexpect2}{wage expectations 3 years after studying: 0=less than 3500 CHF gross per month; 1=3500-4000 CHF; 2=4000-4500 CHF;...; 15=10500-11000 CHF; 16=more than 11000 CHF}
#'   \item{wexpect1othersex}{expected wage of other sex after finishing studies in percent of own expected wage}
#'   \item{wexpect2othersex}{expected wage of other sex 3 years after studying in percent of own expected wage}
#'   \item{male}{1=male; 0=female}
#'   \item{business}{1=BA in business}
#'   \item{econ}{1=BA in economics}
#'   \item{communi}{1=BA in communication}
#'   \item{businform}{1=BA in business informatics}
#'   \item{plansfull}{1=plans working fulltime after studies}
#'   \item{planseduc}{1=plans obtaining further education (e.g. MA) after studies}
#'   \item{sectorcons}{1=planned sector: construction}
#'   \item{sectortradesales}{1=planned sector: trade and sales}
#'   \item{sectortransware}{1=planned sector: transport and warehousing}
#'   \item{sectorhosprest}{1=planned sector: hospitality and restaurant}
#'   \item{sectorinfocom}{1=planned sector: information and communication}
#'   \item{sectorfininsur}{1=planned sector: finance and insurance}
#'   \item{sectorconsult}{1=planned sector: consulting}
#'   \item{sectoreduscience}{1=planned sector: education and science}
#'   \item{sectorhealthsocial}{1=planned sector: health and social services}
#'   \item{typegenstratman}{1=planned job type: general or strategic management}
#'   \item{typemarketing}{1=planned job type: marketing}
#'   \item{typecontrol}{1=planned job type: controlling}
#'   \item{typefinance}{1=planned job type: finance}
#'   \item{typesales}{1=planned job type: sales}
#'   \item{typetechengin}{1=planned job type: technical/engineering}
#'   \item{typehumanres}{1=planned job type: human resources}
#'   \item{posmanager}{1=planned position: manager}
#'   \item{age}{age in years}
#'   \item{swiss}{1=Swiss nationality}
#'   \item{hassiblings}{1=has one or more siblings}
#'   \item{motherhighedu}{1=mother has higher education}
#'   \item{fatherhighedu}{1=father has higher education}
#'   \item{motherworkedfull}{1=mother worked fulltime at respondent's age 4-6}
#'   \item{motherworkedpart}{1=mother worked parttime at respondent's age 4-6}
#'   \item{matwellbeing}{self-assessed material wellbeing compared to average Swiss: 1=much worse; 2=worse; 3=as average Swiss; 4=better; 5=much better}
#'   \item{homeowner}{1=home ownership}
#'   \item{treatmentinformation}{1=if information on median wages in Switzerland was provided (randomized treatment)}
#'   \item{treatmentorder}{1=if order of questions on professional plans and personal information in survey has been reversed (randomized treatment), meaning that personal questions are asked first and professional ones later}
#' }
#' @docType data
#' @references Fernandes, A., Huber, M., and Vaccaro, G. (2020): "Gender Differences in Wage Expectations", arXiv preprint arXiv:2003.11496.
#' @examples
#' data(wexpect)
#' attach(wexpect)
#' # effect of randomized wage information (treatment) on wage expectations 3 years after
#' # studying (outcome)
#' treatweight(y=wexpect2,d=treatmentinformation,x=cbind(male,business,econ,communi,
#' businform,age,swiss,motherhighedu,fatherhighedu),boot=199)
#' # direct effect of gender (treatment) and indirect effect through choice of field of
#' # studies (mediator) on wage expectations (outcome)
#' medweight(y=wexpect2,d=male,m=cbind(business,econ,communi,businform),
#' x=cbind(treatmentinformation,age,swiss,motherhighedu,fatherhighedu),boot=199)
"wexpect"

