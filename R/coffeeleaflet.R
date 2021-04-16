#' Information leaflet on coffee production and environmental awareness of high school / university students in Bulgaria
#'
#' A dataset on the impact of an information leaflet about coffee production on students' awareness about environmental issues
#' collected at Bulgarian highschools and universities in the year 2015.
#'
#' @format A data frame with 522 rows and 48 variables:
#' \describe{
#'   \item{grade}{school grade}
#'   \item{sex}{1=male, 0=female}
#'   \item{age}{age in years}
#'   \item{mob}{month of birth}
#'   \item{bulgnationality}{dummy for Bulgarian nationality}
#'   \item{langbulg}{dummy for Bulgarian mother tongue}
#'   \item{mumage}{mother's age in years}
#'   \item{mumedu}{mother's education (1=lower secondary or less, 2=upper secondary, 3=higher)}
#'   \item{mumprof}{mother's profession (1=manager, 2=specialist, 3=worker, 4=self-employed,5=not working,6=retired,7=other)}
#'   \item{dadage}{father's age in years}
#'   \item{dadedu}{father's education (1=lower secondary or less, 2=upper secondary, 3=higher)}
#'   \item{dadprof}{father's profession (1=manager, 2=specialist, 3=worker, 4=self-employed,5=not working,6=retired,7=other)}
#'   \item{material}{material situation of the family (1=very bad,...,5=very good)}
#'   \item{withbothpar}{dummy for living with both parents}
#'   \item{withmum}{dummy for living with mother only}
#'   \item{withdad}{dummy for living with father only}
#'   \item{withneither}{dummy for living with neither mother nor father}
#'   \item{oldsiblings}{number of older siblings}
#'   \item{youngsiblings}{numer of younger siblings}
#'   \item{schoolmaths}{school dummy (for highschool with maths specialization)}
#'   \item{schoolrakdelsve}{school dummy}
#'   \item{schoolvazov}{school dummy}
#'   \item{schoolfinance}{school dummy}
#'   \item{schoolvarna}{school dummy (for highschool in city of Varna)}
#'   \item{schoolspanish}{school dummy (for Spanish highschool)}
#'   \item{schooltechuni}{school dummy (for technical university)}
#'   \item{schoolvidin}{school dummy (for highschool in city of Vidin)}
#'   \item{schooluni}{school dummy (for university)}
#'   \item{citysofia}{dummy for the capital city of Sofia}
#'   \item{cityvarna}{dummy for the city of Varna}
#'   \item{cityvidin}{dummy for the city of Vidin}
#'   \item{treatment}{treatment (1=leaflet on environmental impact of coffee growing, 0=control group)}
#'   \item{drinkcoffee}{drinks coffee (1=never, 2=not more than 1 time per week,3=several times per week, 4=1 time per day, 5=several times per day)}
#'   \item{cupsest}{outcome: guess how many cups of coffee per capita are consumed in Bulgaria per year}
#'   \item{devi_cupsest}{outcome: deviation of guess from true coffee consumption per capita and year in Bulgaria}
#'   \item{impworldecon}{outcome: assess the importance of coffee for world economy (1=not at all important,..., 5=very important)}
#'   \item{impincome}{assess the importance of coffee as a source of income for people in Africa and Latin America (1=not at all important,..., 5=very important)}
#'   \item{awarewaste}{outcome: awareness of waste production due to coffee production (1=not aware,..., 5=fully aware)}
#'   \item{awarepesticide}{outcome: awareness of pesticide use due to coffee production (1=not aware,..., 5=fully aware)}
#'   \item{awaredeforestation}{outcome: awareness of deforestation due to coffee production (1=not aware,..., 5=fully aware)}
#'   \item{awarewastewater}{outcome: awareness of waste water due to coffee production (1=not aware,..., 5=fully aware)}
#'   \item{awarebiodiversityloss}{outcome: awareness of biodiversity loss due to coffee production (1=not aware,..., 5=fully aware)}
#'   \item{awareunfairworking}{outcome: awareness of unfair working conditions due to coffee production (1=not aware,..., 5=fully aware)}
#'   \item{reusepurposeful}{outcome: can coffee waste be reused purposefully (1=no, 2=maybe, 3=yes)}
#'   \item{reusesoil}{outcome: can coffee waste be reused as soil (1=no, 2=maybe, 3=yes)}
#'   \item{choiceprice}{importance of price when buying coffee (1=not important at all,..., 5=very important, 6=I don't drink coffee)}
#'   \item{choicetastepleasure}{importance of pleasure or taste when buying coffee (1=not important at all,..., 5=very important, 6=I don't drink coffee)}
#'   \item{choiceenvironsocial}{importance of environmental or social impact when buying coffee (1=not important at all,..., 5=very important, 6=I don't drink coffee)}
#' }
#' @docType data
#' @references Faldzhiyskiy, S. (Ecosystem Europe, Bulgaria) and Huber, M. (University of Fribourg): "The impact of an information leaflet about coffee production on students' awareness about environmental issues".
#' @examples
#' \dontrun{
#' data(coffeeleaflet)
#' attach(coffeeleaflet)
#' data=na.omit(cbind(awarewaste,treatment,grade,sex,age))
#' # effect of information leaflet (treatment) on awareness of waste production
#' treatweight(y=data[,1],d=data[,2],x=data[,3:5],boot=199)}
"coffeeleaflet"

