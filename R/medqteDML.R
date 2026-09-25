#' Natural direct and indirect quantile treatment effects with double machine learning
#' @description Estimation of natural direct and indirect quantile treatment effects (QTEs) of a binary treatment operating through a scalar mediator, based on double machine learning. The method inverts post-lasso-based double/debiased machine learning estimates of the cumulative distribution functions (cdfs) of the four potential outcomes \code{Y(1,M(1))}, \code{Y(1,M(0))}, \code{Y(0,M(1))}, \code{Y(0,M(0))}, computed from the efficient influence functions of Hsu, Huber, and Yen (2026), and provides pointwise and uniform confidence bands via a multiplier bootstrap.
#' @param y Outcome/dependent variable, must be a numeric scalar variable, must not contain missings.
#' @param d Treatment, must be binary (either 1 or 0), must not contain missings.
#' @param m Mediator, must be a numeric scalar variable, must not contain missings.
#' @param x (Potential) pre-treatment confounders of the treatment, mediator, and/or outcome, must not contain missings.
#' @param tau Vector of quantile ranks at which the QTEs are estimated, each strictly between 0 and 1. Default is \code{seq(0.2,0.9,0.1)}.
#' @param a Grid of values at which the cdfs of the potential outcomes are estimated before inversion into quantiles. Default is \code{NULL}, in which case a grid of at least 50 points is built from empirical quantiles of \code{y}, spanning somewhat beyond the range of \code{tau}: the potential-outcome distributions are shifted relative to the pooled distribution of \code{y}, so a grid tied only to \code{tau} itself is usually too narrow and forces the quantile inversion to extrapolate. A custom \code{a} can be supplied for more control, but should be finer and wider than \code{tau}; \code{medqteDML} warns if the estimated cdfs do not bracket the requested \code{tau}.
#' @param kfold Number of folds in the cross-fitting procedure. Default is 5.
#' @param trim Trimming threshold applied to the estimated propensity scores \code{P(D=1|X)} and \code{P(D=1|M,X)}, which are winsorized to \code{[trim, 1-trim]}. Default is 0.05.
#' @param nboot Number of multiplier bootstrap replications used for the pointwise confidence intervals and uniform confidence bands. Default is 1000. The bootstrap reuses the estimated efficient influence functions rather than refitting the nuisance parameters, so increasing \code{nboot} is inexpensive relative to increasing \code{length(a)} or \code{kfold}.
#' @param alpha Significance level for the confidence intervals and bands. Default is 0.05.
#' @param q Lower quantile used to compute the rescaled quantile spread underlying the uniform confidence bands, see Chernozhukov, Fernandez-Val, and Melly (2013). Default is 0.1, a common choice in that literature.
#' @details The nuisance parameters (propensity scores \code{P(D=1|X)} and \code{P(D=1|M,X)}, and the conditional cdf \code{P(Y<=a|D,M,X)} estimated by distribution regression) are estimated by post-lasso, with covariates selected via a double/triple selection step (lasso selection from the \code{D|X}, \code{D|M,X}, and \code{1{Y<=a}|D,M,X} models, unioned) at every grid point \code{a} and every fold, following the implementation underlying the empirical application in Hsu, Huber, and Yen (2026). The estimator itself is the K-fold cross-fitting estimator of the cdf of each potential outcome based on the (triply robust) efficient influence function, see Section 2 of Hsu, Huber, and Yen (2026) and, for the underlying regression imputation approach, Algorithm 2 in Farbmacher et al. (2022). Quantiles are obtained by inverting the estimated cdfs on the grid \code{a}, after a rearrangement step (Chernozhukov, Fernandez-Val, and Galichon 2010) that enforces monotonicity. Statistical inference is based on a multiplier bootstrap that reuses the estimated efficient influence functions, following Section 2.5 of Hsu, Huber, and Yen (2026); both pointwise confidence intervals (percentile/reflection method) and uniform confidence bands (based on the bootstrap Kolmogorov-Smirnov max-t-statistic and the rescaled quantile spread) are provided.
#' @return A \code{medqteDML} object contains the following components:
#' @return \code{effects}: a list of five data frames, one per estimand: \code{total} (the total effect, \code{Y(1,M(1))} minus \code{Y(0,M(0))}); \code{dir.treat} and \code{dir.control} (the direct effects with the mediator fixed at its value under treatment, \code{M(1)}, or under non-treatment, \code{M(0)}, respectively -- \code{NDQTE'} and \code{NDQTE} in Hsu, Huber, and Yen 2026); and \code{indir.treat} and \code{indir.control} (the indirect effects with the treatment fixed at \code{D=1} or \code{D=0}, respectively -- \code{NIQTE} and \code{NIQTE'} in Hsu, Huber, and Yen 2026). Each data frame has columns \code{tau} (the quantile rank), \code{effect} (the point estimate), \code{se} (bootstrap standard deviation), \code{pw.lower}/\code{pw.upper} (pointwise confidence interval), and \code{unif.lower}/\code{unif.upper} (uniform confidence band).
#' @return \code{cdf}: a data frame with the estimated cdfs of the four potential outcomes across the grid \code{a}.
#' @return \code{quantiles}: a data frame with the estimated quantiles \code{Q(1,M(1))}, \code{Q(1,M(0))}, \code{Q(0,M(1))}, \code{Q(0,M(0))} at each \code{tau}.
#' @return \code{boot}: a list with the raw multiplier bootstrap draws of the five QTEs (rows are \code{tau}, columns are bootstrap replications), for users who want to construct alternative confidence measures.
#' @return \code{ntrimmed}: number of observations whose estimated propensity score \code{P(D=1|X)} or \code{P(D=1|M,X)} was winsorized (capped to \code{trim} or \code{1-trim}) in at least one fold or grid point.
#' @references Chernozhukov, V., Fernandez-Val, I., and Galichon, A. (2010): "Quantile and Probability Curves Without Crossing", Econometrica, 78, 1093-1125.
#' @references Chernozhukov, V., Fernandez-Val, I., and Melly, B. (2013): "Inference on Counterfactual Distributions", Econometrica, 81, 2205-2268.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., Robins, J. (2018): "Double/debiased machine learning for treatment and structural parameters", The Econometrics Journal, 21, C1-C68.
#' @references Farbmacher, H., Huber, M., Laffers, L., Langen, H., and Spindler, M. (2022): "Causal mediation analysis with double machine learning", The Econometrics Journal, 25, 277-300.
#' @references Hsu, Y.-C., Huber, M., and Yen, Y.-M. (2026): "Estimation of Direct and Indirect Quantile Treatment Effects with Double Machine Learning", Journal of Business & Economic Statistics, DOI: 10.1080/07350015.2026.2654889.
#' @examples # A little example with simulated data
#' \dontrun{
#' n=5000                                               # sample size
#' p=10                                                 # number of covariates
#' x=matrix(rnorm(n*p),ncol=p)                          # covariate matrix
#' d=rbinom(n,1,plogis(0.5*x[,1]-0.3*x[,2]))            # treatment equation
#' m=rbinom(n,1,plogis(-0.3+0.8*d+0.4*x[,1]-0.3*x[,2])) # mediator equation
#' y=2*d+1.5*m-0.4*x[,1]+0.4*x[,2]+rnorm(n)             # outcome equation
#' # Direct effect is 2 at every quantile; indirect effect is around 0.3
#' output=medqteDML(y=y,d=d,m=m,x=x,tau=c(0.3,0.5,0.7))
#' output$effects$dir.control
#' output$effects$indir.treat
#' output$effects$dir.treat
#' output$effects$indir.control
#' output$ntrimmed
#' }
#' @importFrom stats binomial glm lm predict quantile approx aggregate rnorm sd qnorm
#' @import hdm
#' @export
medqteDML=function(y, d, m, x, tau=seq(0.2,0.9,0.1), a=NULL, kfold=5, trim=0.05, nboot=1000, alpha=0.05, q=0.1){

  if (is.null(a)){
    ngrid=max(50, 5*length(tau))
    pgrid=seq(max(0.001,min(tau)-0.2), min(0.999,max(tau)+0.2), length.out=ngrid)
    a=sort(unique(as.numeric(quantile(y, probs=pgrid, na.rm=TRUE, type=1))))
  }

  fit=medqteplassofit(y=y, d=d, m=m, x=x, a=a, kfold=kfold, trim=trim)

  ## diagnostic: warn if any potential-outcome cdf fails to bracket a requested
  ## tau on the grid a, since findqx then has to extrapolate (approx(...,rule=2)),
  ## which silently clamps to the grid boundary rather than a real quantile
  cdfs=list(F11=fit$px11, F10=fit$px10, F00=fit$px00, F01=fit$px01)
  for (nm in names(cdfs)){
    rng=range(cdfs[[nm]])
    if (min(tau) < rng[1] | max(tau) > rng[2]) warning(sprintf(
      "requested tau (range %.3f-%.3f) falls outside the estimated cdf range for %s (%.3f-%.3f); the corresponding quantile(s) were extrapolated and clamped to a grid boundary rather than genuinely estimated. Widen a (or its default padding) to fix.",
      min(tau), max(tau), nm, rng[1], rng[2]))
  }

  q11=findqx(px=fit$px11, ax=a, taux=tau); q10=findqx(px=fit$px10, ax=a, taux=tau)
  q00=findqx(px=fit$px00, ax=a, taux=tau); q01=findqx(px=fit$px01, ax=a, taux=tau)

  TQTE=q11-q00; NDQTE=q10-q00; NIQTE=q11-q10; NDQTE1=q11-q01; NIQTE1=q01-q00

  n=length(y)
  boot.TQTE=matrix(0,length(tau),nboot); boot.NDQTE=matrix(0,length(tau),nboot); boot.NIQTE=matrix(0,length(tau),nboot)
  boot.NDQTE1=matrix(0,length(tau),nboot); boot.NIQTE1=matrix(0,length(tau),nboot)

  for (b in 1:nboot){
    xi=rnorm(n)     # E[xi]=0, Var(xi)=1, E[exp(|xi|)]<Inf, as required for the multiplier bootstrap
    pb11=mbootp(p=fit$px11, psi=fit$if11, xi=xi); pb10=mbootp(p=fit$px10, psi=fit$if10, xi=xi)
    pb00=mbootp(p=fit$px00, psi=fit$if00, xi=xi); pb01=mbootp(p=fit$px01, psi=fit$if01, xi=xi)

    qb11=findqx(px=pb11, ax=a, taux=tau); qb10=findqx(px=pb10, ax=a, taux=tau)
    qb00=findqx(px=pb00, ax=a, taux=tau); qb01=findqx(px=pb01, ax=a, taux=tau)

    boot.TQTE[,b]=qb11-qb00; boot.NDQTE[,b]=qb10-qb00; boot.NIQTE[,b]=qb11-qb10
    boot.NDQTE1[,b]=qb11-qb01; boot.NIQTE1[,b]=qb01-qb00
  }

  mkeffects=function(est, datab){
    pw=pciper(datab=datab, est=est, alpha=alpha)
    unif=sciqs(datab=datab, est=est, q=q, alpha=alpha)
    data.frame(tau=tau, effect=est, se=apply(datab,1,sd), pw.lower=pw[,1], pw.upper=pw[,2], unif.lower=unif[,1], unif.upper=unif[,2])
  }

  effects=list(total=mkeffects(TQTE,boot.TQTE), dir.treat=mkeffects(NDQTE1,boot.NDQTE1), dir.control=mkeffects(NDQTE,boot.NDQTE),
               indir.treat=mkeffects(NIQTE,boot.NIQTE), indir.control=mkeffects(NIQTE1,boot.NIQTE1))

  cdf=data.frame(a=a, F11=fit$px11, F10=fit$px10, F00=fit$px00, F01=fit$px01)
  quantiles=data.frame(tau=tau, Q11=q11, Q10=q10, Q00=q00, Q01=q01)
  boot=list(total=boot.TQTE, dir.treat=boot.NDQTE1, dir.control=boot.NDQTE, indir.treat=boot.NIQTE, indir.control=boot.NIQTE1)
  ntrimmed=sum(apply(fit$trimmed,1,any))   # units winsorized in at least one grid point/fold, see @return

  list(effects=effects, cdf=cdf, quantiles=quantiles, boot=boot, ntrimmed=ntrimmed)
}
