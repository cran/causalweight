#' Data on daily spending and coupon receipt 
#' A dataset containing information on the purchasing behavior of 1582 retail store customers across 32 coupon campaigns.
#'
#' @format A data frame with 50,624 rows and 27 variables:
#' \describe{
#'   \item{customer_id}{customer identifier}
#'   \item{period}{period of observation: 1 = 1st period to 32 = last period}
#'   \item{age_range}{age of customer: 1 = 18-25; 2 = 26-35; 3 = 36-45; 4 = 46-55; 5 = 56-70; 6 = 71 plus}
#'   \item{married}{marital status: 1 = married; 0 = unmarried}
#'   \item{rented}{dwelling type: 1 = rented; 0 = owned}
#'   \item{family_size}{number of family members: 1 = 1; 2 = 2; 3 = 3; 4 = 4; 5 = 5 plus}
#'   \item{income_bracket}{income group: 1 = lowest to 12 = highest}
#'   \item{dailyspending_preperiod}{customer's daily spending at the retailer in previous period}
#'   \item{purchase_ReadyEatFood_preperiod}{purchases of ready-to-eat food in previous period: 1 = yes, 0 = no}
#'   \item{purchase_MeatSeafood_preperiod}{purchases of meat and seafood products in previous period: 1 = yes, 0 = no}
#'   \item{purchase_OtherFood_preperiod}{purchases of other food products in previous period: 1 = yes, 0 = no}
#'   \item{purchase_Drugstore_preperiod}{purchases of drugstore products in previous period: 1 = yes, 0 = no}
#'   \item{purchase_OtherNonfood_preperiod}{purchases of other non-food products in previous period: 1 = yes, 0 = no}
#'   \item{coupons_Any_preperiod}{coupon reception in previous period: 1 = customer received at least one coupon; 0 = customer did not receive any coupon}
#'   \item{coupons_ReadyEatFood_preperiod}{coupon reception in previous period: 1 = customer received at least one ready-to-eat food coupon; 0 = customer did not receive any ready-to-eat food coupon}
#'   \item{coupons_MeatSeafood_preperiod}{coupon reception in previous period: 1 = customer received at least one meat/seafood coupon; 0 = customer did not receive any meat/seafood coupon}
#'   \item{coupons_OtherFood_preperiod}{coupon reception in previous period: 1 = customer received at least one coupon applicable to other food items; 0 = customer did not receive any coupon applicable to other food items}
#'   \item{coupons_Drugstore_preperiod}{coupon reception in previous period: 1 = customer received at least one drugstore coupon; 0 = customer did not receive any drugstore coupon}
#'   \item{coupons_OtherNonfood_preperiod}{coupon reception in previous period: 1 = customer received at least one coupon applicable to other non-food items; 0 = customer did not receive any coupon applicable to other non-food items}
#'   \item{coupons_Any_redeemed_preperiod}{coupon redemption in previous period: 1 = customer redeemed at least one coupon; 0 = customer did not redeem any coupon}
#'   \item{coupons_Any}{treatment: 1 = customer received at least one coupon in current period; 0 = customer did not receive any coupon}
#'   \item{coupons_ReadyEatFood}{treatment: 1 = customer received at least one ready-to-eat food coupon; 0 = customer did not receive any ready-to-eat food coupon}
#'   \item{coupons_MeatSeafood}{treatment: 1 = customer received at least one meat/seafood coupon; 0 = customer did not receive any meat/seafood coupon}
#'   \item{coupons_OtherFood}{treatment: 1 = customer received at least one coupon applicable to other food items; 0 = customer did not receive any coupon applicable to other food items}
#'   \item{coupons_Drugstore}{treatment: 1 = customer received at least one drugstore coupon; 0 = customer did not receive any drugstore coupon}
#'   \item{coupons_OtherNonfood}{treatment: 1 = customer received at least one coupon applicable to other non-food items; 0 = customer did not receive any coupon applicable to other non-food items}
#'   \item{dailyspending}{outcome: customer's daily spending at the retailer in current period}
#'    }
#' @docType data
#' @references Langen, Henrika, and Huber, Martin (2023): "How causal machine learning can leverage marketing strategies: Assessing and improving the performance of a coupon campaign." PLoS ONE, 18 (1): e0278937. 
"couponsretailer"

