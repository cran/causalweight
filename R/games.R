#' Sales of video games
#'
#' A dataset containing information on 3956 video games, including sales as well as expert and user ratings.
#' @format A data frame with 3956 rows and 9 variables:
#' \describe{
#'   \item{name}{factor variable providing the name of the video game}
#'   \item{genre}{factor variable indicating the genre of the game (e.g. Action, Sports...)}
#'   \item{platform}{factor variable indicating the hardware platform of the game (e.g. PC,...)}
#'   \item{esrbrating}{factor variable indicating the age recommendation for the game(E is age 6+, T is 13+, M is 17+)}
#'   \item{publisher}{factor variable indicating the publisher of the game}
#'   \item{year}{numeric variable indicating the year the video game was released}
#'   \item{metascore}{numeric variable providing a weighted average rating of the game by professional critics}
#'   \item{userscore}{numeric variable providing the average user rating of the game}
#'   \item{sales}{numeric variable indicating the total global sales (in millions) of the game up to the year 2018}
#' }
#' @import fastDummies
#' @docType data
#' @references Wittwer, J. (2020): "Der Erfolg von Videospielen - eine empirische Untersuchung moeglicher Erfolgsfaktoren", BA thesis, University of Fribourg.
#' @examples
#' \dontrun{
#' #load data
#' data(games)
#' #generate missing indicator
#' mis=1*(apply(is.na(games),1,sum)>0)
#' #select non-missing observations
#' games_nomis=games[mis==0,]
#' #turn year into a factor variable
#' games_nomis$year=factor(games_nomis$year)
#' #attach data
#' attach(games_nomis)
#' #load library for generating dummies
#' library(fastDummies)
#' #generate dummies for genre
#' dummies=dummy_cols(genre, remove_most_frequent_dummy = TRUE)
#' #drop original variable
#' genredummies=dummies[,2:ncol(dummies)]
#' #make dummies numeric
#' genredummies=apply(genredummies, 2, function(genredummies) as.numeric(genredummies))
#' #generate dummies for year
#' dummies=dummy_cols(year, remove_most_frequent_dummy = TRUE)
#' #drop original variable
#' yeardummies=dummies[,2:ncol(dummies)]
#' #make dummies numeric
#' yeardummies=apply(yeardummies, 2, function(yeardummies) as.numeric(yeardummies))
#' # mediation analysis with metascore as treatment, userscore as mediator, sales as outcome
#' x=cbind(genredummies,yeardummies)
#' output=medweightcont(y=sales,d=metascore, d0=60, d1=80, m=userscore, x=x, boot=199)
#' round(output$results,3)
#' output$ntrimmed}
"games"

