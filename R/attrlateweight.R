#' Local average treatment effect estimation in multiple follow-up periods with outcome attrition based on inverse probability weighting
#' @description Instrumental variable-based evaluation of local average treatment effects using weighting by the inverse of the instrument propensity score.
#' @param y1 Outcome variable in the first outcome period.
#' @param y2 Outcome variable in the second outcome period.
#' @param s1 Selection indicator for first outcome period. Must be one if \code{y1} is observed (non-missing) and zero if \code{y1} is not observed (missing).
#' @param s2 Selection indicator for second outcome period. Must be one if \code{y1} is observed (non-missing) and zero if \code{y1} is not observed (missing).
#' @param d Treatment, must be binary (either 1 or 0), must not contain missings.
#' @param z Instrument for the endogenous treatment, must be binary (either 1 or 0), must not contain missings.
#' @param x0 Baseline (pre-instrument) confounders of the instrument and outcome, must not contain missings.
#' @param x1 Confounders in outcome period 1 (may include outcomes of period 1 \code{y1})
#' @param weightmax Trimming rule based on the maximum relative weight a single observation may obtain in estimation - observations with higher weights are discarded. Default is 0.1 (no observation can be assigned more than 10 percent of weights)
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @param cluster A cluster ID for block or cluster bootstrapping when units are clustered rather than iid. Must be numerical. Default is NULL (standard bootstrap without clustering).
#' @details  Estimation of local average treatment effects of a binary endogenous treatment on outcomes in two follow up periods that are prone to attrition. Treatment endogeneity is tackled by a binary instrument that is assumed to be conditionally valid given observed baseline confounders \code{x0}. Outcome attrition is tackled by either assuming that it is missing at random (MAR), i.e. selection w.r.t. observed variables  \code{d}, \code{z}, \code{x0}, \code{x1} (in the case of \code{y2}), and \code{s1} (in the case of \code{y2}); or by assuming latent ignorability (LI), i.e. selection w.r.t. the treatment compliance type as well as \code{z}, \code{x0}, \code{x1} (in the case of \code{y2}), and \code{s1} (in the case of \code{y2}). Units are weighted by the inverse of their conditional instrument and selection propensities, which are estimated by probit regression. Standard errors are obtained by bootstrapping the effect.
#' @return An attrlateweight object contains one component \code{results}:
#' @return \code{results}: a 4X4 matrix containing the effect estimates in the first row ("effects"), standard errors in the second row ("se"), p-values in the third row ("p-value"), and the number of trimmed observations due to too large weights in the fourth row ("trimmed obs"). The first column provides the local average treatment effect (LATE) on \code{y1} among compliers under missingness at random (MAR). The second column provides the local average treatment effect (LATE) on \code{y2}  under missingness at random (MAR). The third column provides the local average treatment effect (LATE) on \code{y1}  under latent ignorability (LI). The forth column provides the local average treatment effect (LATE) on \code{y2} under latent ignorability (LI).
#' @references Frölich, M., Huber, M. (2014): "Treatment Evaluation With Multiple Outcome Periods Under Endogeneity and Attrition", Journal of the American Statistical Association, 109, 1697-1711.
#' @examples # A little example with simulated data (4000 observations)
#' n=4000
#' e=(rmvnorm(n,rep(0,3), matrix(c(1,0.3,0.3,  0.3,1,0.3,  0.3,0.3,1),3,3) ))
#' x0=runif(n,0,1)
#' z=(0.25*x0+rnorm(n)>0)*1
#' d=(1.2*z-0.25*x0+e[,1]>0.5)*1
#' y1_star=0.5*x0+0.5*d+e[,2]
#' s1=(0.25*x0+0.25*d+rnorm(n)>-0.5)*1
#' y1=s1*y1_star
#' x1=(0.5*x0+0.5*rnorm(n))
#' y2_star=0.5*x0+x1+d+e[,3]
#' s2=s1*((0.25*x0+0.25*x1+0.25*d+rnorm(n)>-0.5)*1)
#' y2=s2*y2_star
#' # The true LATEs on y1 and y2 are equal to 0.5 and 1, respectively.
#' output=attrlateweight(y1=y1,y2=y2,s1=s1,s2=s2,d=d,z=z,x0=x0,x1=x1, boot=19)
#' round(output$results,3)
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm quantile predict
#' @import mvtnorm
#' @importFrom LARF Generate.Powers
#' @importFrom hdm rlasso rlassologit
#' @export
attrlateweight<-function(y1,y2,s1,s2,d,z,x0,x1, weightmax=0.1,  boot=1999, cluster=NULL){
  y1[s1==0]=0;y2[s2==0]=0
  temp=attrlate(y1=y1,y2=y2,r1=s1,r2=s2,d=d,z=z,x0=x0,x1=x1, weightmax=weightmax)
  ntrimmed=temp[(length(temp)-3):length(temp)]
  temp=temp[1:(length(temp)-4)]
  temp2=bootstrap.attrlate(y1=y1,y2=y2,r1=s1,r2=s2,d=d,z=z,x0=x0,x1=x1,weightmax=weightmax, boot=boot, cluster=cluster)
  se=apply(temp2[,1:(ncol(temp2)-4)], 2, sd)
  temp3=2*pnorm(-abs(temp/se))
  results=rbind(temp, se, temp3,ntrimmed)
  colnames(results)=c("mar.y1", "mar.y2", "li.y1", "li.y2")
  rownames(results)=c("effect", "se", "p-value", "trimmed obs")
  list(results=results)
}

