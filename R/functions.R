hdtreat_newscore = function(y, d, x, MLmethod = "lasso", k = 3, zeta, seed){
  ybin <- 1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  x <- data.frame(x)
  stepsize <- ceiling((1/k)*length(d))
  set.seed(seed)
  idx <- sample(length(d), replace=FALSE)
  score <- c()
  # cross-fitting procedure that splits sample in training and testing data
  for (i in 1:k){
    tesample <- idx[((i-1)*stepsize+1):(min((i)*stepsize,length(d)))]
    trsample <- idx[-tesample]
    eydx <- MLfunct(y = y[trsample], x = x[trsample,], d1 = d[trsample], MLmethod = MLmethod, ybin = ybin)
    eydxte <- predict(eydx, x[tesample,], onlySL = TRUE)$pred  #predict conditional outcome in test data
    score <- rbind(score, cbind(d[tesample],y[tesample],eydxte,zeta[tesample]))
  }
  score <- score[order(idx),]
  score
}

treatDML_newscore = function(y, d, x, dtreat = 1, dcontrol = 0,
  MLmethod = "lasso", k = 3, zeta_sigma = 0.5, seed = 123){
  n <- dim(x)[1]
  zeta <- rnorm(n,0,sd = zeta_sigma) # sample from a normal distribution to avoid degenerated distribution
  dtre <- 1*(d==dtreat)
  dcon <- 1*(d==dcontrol)
  scorestreat <- hdtreat_newscore(y = y,d = dtre, x = x, MLmethod = MLmethod, k = k, zeta = zeta, seed = seed)
  scorescontrol <- hdtreat_newscore(y = y, d = dcon, x = x, MLmethod = MLmethod, k = k ,zeta = zeta, seed = seed)
  tscores <- scorestreat[,3]
  cscores <- scorescontrol[,3]
  meantreat <- mean(tscores)
  meancontrol <- mean(cscores)
  effect <- mean((tscores - cscores)^2+scorestreat[,4])
  se <- sqrt((mean(((tscores-cscores)^2-effect)^2)+var(scorestreat[,4]))/length(tscores))
  pval <- 2*pnorm((-1)*abs(effect/se))
  list(effect = effect, se = se, pval = pval, meantreat = meantreat,
    meancontrol = meancontrol)
}

computescores = function(Y,X,D){
  pscoreforest <- regression_forest(X = X,Y = c(D), num.trees = 200)
  pscore <- predict(pscoreforest,X)$predictions
  yforest1 <- regression_forest(X = X[D==1,], Y = Y[D==1], num.trees = 200)
  condy1 <- predict(yforest1,X)$predictions
  yforest0 <- regression_forest(X = X[D==0,], Y = Y[D==0], num.trees = 200)
  condy0 <- predict(yforest0,X)$predictions
  n <- length(D)
  weightsum1 <- sum(D/pscore)
  weightsum0 <- sum((1-D)/(1-pscore))
  scores1 <- (n*D*(Y-condy1)/(pscore))/weightsum1+(condy1)
  scores0 <- (n*(1-D)*(Y-condy0)/(1-pscore))/weightsum0+(condy0)
  c(scores1-scores0)
}

# uniform inference using the doubly robust score based on bootstrap:

treatDML_bootstrap=function(y,d,x, dtreat = 1, dcontrol = 0, seed = 123, s = NULL, normalized = TRUE, trim = 0.01,
  MLmethod = "lasso", k = 3, B = 2000, importance = 0.95, alpha = 0.1, share = 0.5){
  # add sample split
  idx <- sample(length(d),length(d)*share,replace=FALSE)
  d1 <- d[idx]
  d2 <- d[-idx]
  y1 <- y[idx]
  y2 <- y[-idx]
  x1 <- x[idx,]
  x2 <- x[-idx,]
  # determine subsets with large heterogeneity
  scorestr <- computescores(Y = y1, X = x1, D = d1)
  data1 <- data.frame(scorestr,x1)
  randomf <- ranger(scorestr~., data = data1, num.trees = 200)
  quantilevariable <- which(ranger::importance(randomf)>=quantile(ranger::importance(randomf),importance))
  quantilevariable <-  quantilevariable - 1 # because of the intercept!
  nm <- length(quantilevariable)
  index <- matrix(NA,length(y2),2*nm)
  for (j in 1:nm){
    index[,2*j-1] <- x2[,quantilevariable[j]]<quantile(x2[,quantilevariable[j]],0.5) # define subsets
    index[,2*j] <- x2[,quantilevariable[j]]>=quantile(x2[,quantilevariable[j]],0.5)
  }
  # bootstrap
  M <- dim(index)[2]
  n2 <- length(y2)
  bootscore <- matrix(NA,n2,M)
  results_effect <- rep(NA,M)
  results_pval <- rep(NA,M)
  results_est <- rep(NA,M)
  samplesize_subgroup <- rep(NA,M)
  # estimation using DR score for each subset (using other part of data)
  for (m in 1:M){
    dtre <- 1*(d2[index[,m]]==dtreat)
    dcon <- 1*(d2[index[,m]]==dcontrol)
    scorestreat <- hdtreat(y = y2[index[,m]], d = dtre, x = x2[index[,m],], s = s, trim = trim, MLmethod = MLmethod, k = k)
    scorescontrol <- hdtreat(y = y2[index[,m]], d = dcon, x = x2[index[,m],], s = s, trim = trim, MLmethod = MLmethod,k = k)
    trimmed <- 1*(scorescontrol[,7]+scorestreat[,7]>0)  #number of trimmed observations
    scorestreat <- scorestreat[trimmed==0,]
    scorescontrol <- scorescontrol[trimmed==0,]
    if (normalized==FALSE){
      tscores <- (scorestreat[,1]*scorestreat[,2]*(scorestreat[,3]-scorestreat[,4])/(scorestreat[,5])+scorestreat[,6]*scorestreat[,4])/mean(scorestreat[,6])
      cscores <- (scorescontrol[,1]*scorescontrol[,2]*(scorescontrol[,3]-scorescontrol[,4])/(scorescontrol[,5])+scorescontrol[,6]*scorescontrol[,4])/mean(scorescontrol[,6])
    }
    if (normalized!=FALSE){
      ntreat <- nrow(scorestreat)
      weightsumtreat <- sum(scorestreat[,1]*scorestreat[,2]/(scorestreat[,5]))
      tscores <- (ntreat*scorestreat[,1]*scorestreat[,2]*(scorestreat[,3]-scorestreat[,4])/(scorestreat[,5]))/weightsumtreat+(scorestreat[,6]*scorestreat[,4])/mean(scorestreat[,6])
      ncontrol <- nrow(scorescontrol)
      weightsumcontrol <- sum(scorescontrol[,1]*scorescontrol[,2]/(scorescontrol[,5]))
      cscores <- (ncontrol*scorescontrol[,1]*scorescontrol[,2]*(scorescontrol[,3]-scorescontrol[,4])/(scorescontrol[,5]))/weightsumcontrol+(scorescontrol[,6]*scorescontrol[,4])/mean(scorescontrol[,6])
    }
    meantreat <- mean(tscores)
    meancontrol <- mean(cscores)
    effect <- meantreat - meancontrol
    se <- sqrt(mean((tscores-cscores-effect)^2))
    pval <- 2*pnorm((-1)*abs(sqrt(length(tscores))*effect/se))
    results_est[m] <- sqrt(length(tscores))*se^{-1}*abs(effect)
    results_effect[m] <- effect
    results_pval[m] <- pval
    bootscore[1:length(tscores),m] <- se^(-1)*((tscores-cscores)-effect)
    samplesize_subgroup[m] <- length(tscores)
  }
  process <- matrix(NA,B,M)
  sup_bootstrap_max <- rep(NA,B)
  sup_bootstrap_L2 <- rep(NA,B)
  for (b in 1:B){
    epsilon <- rnorm(n2)
    for (m in 1:M){
      process[b,m] <- apply(matrix(bootscore[,m]),2,function(x) samplesize_subgroup[m]^(-1/2)*sum(epsilon*x,na.rm=T))
    }
    sup_bootstrap_max[b] <- max(abs(process[b,]))
    sup_bootstrap_L2[b] <- sqrt(sum(process[b,]^2))
  }
  reject_max <- max(results_est)>quantile(sup_bootstrap_max,1-alpha)
  reject_L2 <- sqrt(sum(results_est^2))>quantile(sup_bootstrap_L2,1-alpha)
  reject_standard <- sum(results_est>qnorm(1-(alpha/(2*M))))>=1
  list(test_reject = reject_max, effect = results_effect, pval = results_pval, test_standard = reject_standard, test_L2 = reject_L2)
}

mediation<-function(y,d,m,x,w=NULL,s=NULL,z=NULL, selpop=FALSE, trim=0.05, ATET=FALSE, logit=FALSE){
  if (is.null(w)==TRUE){
    if (logit==FALSE){
      if (is.null(z)==TRUE){
        pscore.mx=glm(d~cbind(m,x),family=binomial(probit))$fitted
        pscore.x=glm(d~x,family=binomial(probit))$fitted
      }
      if (is.null(s)==FALSE & is.null(z)==TRUE) pscore.s=glm(s~cbind(d,m,x),family=binomial(probit))$fitted
      if (is.null(s)==FALSE & is.null(z)==FALSE) {
        pscore.s=glm(s~cbind(d,m,x,z),family=binomial(probit))$fitted
        pscore.mx=glm(d~cbind(m,x,pscore.s),family=binomial(probit))$fitted
        pscore.x=glm(d~cbind(x,pscore.s),family=binomial(probit))$fitted
      }
    }
    if (logit==TRUE){
      if (is.null(z)==TRUE){
        pscore.mx=glm(d~cbind(m,x),family=binomial(logit))$fitted
        pscore.x=glm(d~x,family=binomial(logit))$fitted
      }
      if (is.null(s)==FALSE & is.null(z)==TRUE) pscore.s=glm(s~cbind(d,m,x),family=binomial(logit))$fitted
      if (is.null(s)==FALSE & is.null(z)==FALSE) {
        pscore.s=glm(s~cbind(d,m,x,z),family=binomial(logit))$fitted
        pscore.mx=glm(d~cbind(m,x,pscore.s),family=binomial(logit))$fitted
        pscore.x=glm(d~cbind(x,pscore.s),family=binomial(logit))$fitted
      }
    }
    if (is.null(s)==TRUE | selpop==TRUE){
      ind=((pscore.mx<trim) | (pscore.mx>(1-trim)) )
      y=y[ind==0]; d=d[ind==0]; pscore.mx=pscore.mx[ind==0]; pscore.x=pscore.x[ind==0]
    }
    if (is.null(s)==FALSE & selpop==FALSE) {
      ind=((pscore.mx<trim) | (pscore.mx>(1-trim)) | pscore.s<trim)
      y=y[ind==0]; d=d[ind==0]; s=s[ind==0]; pscore.mx=pscore.mx[ind==0]; pscore.x=pscore.x[ind==0]; pscore.s=pscore.s[ind==0]
    }
    if ((is.null(s)==TRUE) | (is.null(s)==FALSE & selpop==TRUE)){
      if (is.null(s)==FALSE & selpop==TRUE) {y=y[s==1]; d=d[s==1]; pscore.mx=pscore.mx[s==1]; pscore.x=pscore.x[s==1]}
      if (ATET==FALSE){
        y1m1=sum(y*d/pscore.x)/sum(d/pscore.x)
        y1m0=sum(y*d*(1-pscore.mx)/(pscore.mx*(1-pscore.x)))/sum(d*(1-pscore.mx)/(pscore.mx*(1-pscore.x)))
        y0m0=sum(y*(1-d)/(1-pscore.x))/sum((1-d)/(1-pscore.x))
        y0m1=sum(y*(1-d)*pscore.mx/((1-pscore.mx)*pscore.x))/sum((1-d)*pscore.mx/((1-pscore.mx)*pscore.x))
      }
      if (ATET==TRUE){
        y1m1=sum(y*d)/sum(d)
        y1m0=sum(y*d*(1-pscore.mx)*pscore.x/(pscore.mx*(1-pscore.x)))/sum(d*(1-pscore.mx)*pscore.x/(pscore.mx*(1-pscore.x)))
        y0m0=sum(y*(1-d)*pscore.x/(1-pscore.x))/sum((1-d)*pscore.x/(1-pscore.x))
        y0m1=sum(y*(1-d)*pscore.mx/((1-pscore.mx)))/sum((1-d)*pscore.mx/((1-pscore.mx)))
      }
    }

   if (is.null(s)==FALSE & selpop==FALSE){
      if (ATET==FALSE){
        y1m1=sum(y*d*s/(pscore.s*pscore.x))/sum(d*s/(pscore.s*pscore.x))
        y1m0=sum(y*d*s*(1-pscore.mx)/(pscore.s*pscore.mx*(1-pscore.x)))/sum(d*s*(1-pscore.mx)/(pscore.s*pscore.mx*(1-pscore.x)))
        y0m0=sum(y*(1-d)*s/(pscore.s*(1-pscore.x)))/sum((1-d)*s/(pscore.s*(1-pscore.x)))
        y0m1=sum(y*(1-d)*s*pscore.mx/(pscore.s*(1-pscore.mx)*pscore.x))/sum((1-d)*s*pscore.mx/(pscore.s*(1-pscore.mx)*pscore.x))
      }
      if (ATET==TRUE){
        y1m1=sum(y*d*s/pscore.s)/sum(d*s/pscore.s)
        y1m0=sum(y*d*s*(1-pscore.mx)*pscore.x/(pscore.s*pscore.mx*(1-pscore.x)))/sum(d*s*(1-pscore.mx)*pscore.x/(pscore.s*pscore.mx*(1-pscore.x)))
        y0m0=sum(y*(1-d)*s*pscore.x/(pscore.s*(1-pscore.x)))/sum((1-d)*s*pscore.x/(pscore.s*(1-pscore.x)))
        y0m1=sum(y*(1-d)*s*pscore.mx/(pscore.s*(1-pscore.mx)))/sum((1-d)*s*pscore.mx/(pscore.s*(1-pscore.mx)))
      }
    }
    results=c(y1m1-y0m0, y1m1 - y0m1, y1m0-y0m0, y1m1 - y1m0, y0m1 - y0m0, sum(ind))
  }

  if (is.null(w)==FALSE){
    if (logit==FALSE){
      pscore.mwx=glm(d~cbind(m,w,x),family=binomial(probit))$fitted
      pscore.x=glm(d~x,family=binomial(probit))$fitted
      pscore.wx=glm(d~cbind(w,x),family=binomial(probit))$fitted
      pscore.mx=glm(d~cbind(m,x),family=binomial(probit))$fitted
    }
    if (logit==TRUE){
      pscore.mwx=glm(d~cbind(m,w,x),family=binomial(logit))$fitted
      pscore.x=glm(d~x,family=binomial(logit))$fitted
      pscore.wx=glm(d~cbind(w,x),family=binomial(logit))$fitted
      pscore.mx=glm(d~cbind(m,x),family=binomial(logit))$fitted
    }
    ind=((pscore.mwx<trim) | (pscore.mwx>(1-trim)) )
    y=y[ind==0]; d=d[ind==0]; pscore.mx=pscore.mx[ind==0]; pscore.x=pscore.x[ind==0]; pscore.mwx=pscore.mwx[ind==0]; pscore.wx=pscore.wx[ind==0]

    if (ATET==FALSE){
      y1m1=sum(y*d/pscore.x)/sum(d/pscore.x)
      y1m0<-(sum(y*d*(1-pscore.mwx)/((1-pscore.x)*pscore.mwx))/sum(d*(1-pscore.mwx)/((1-pscore.x)*pscore.mwx)))
      y0m0=sum(y*(1-d)/(1-pscore.x))/sum((1-d)/(1-pscore.x))
      y0m1<-(sum(y*(1-d)* pscore.mwx/(pscore.x*(1-pscore.mwx)))/sum((1-d)* pscore.mwx/(pscore.x*(1-pscore.mwx))))
      y1m0p=(sum( y*d/pscore.mwx * (1-pscore.mwx)/(1-pscore.wx)* (pscore.wx)/pscore.x )/sum(d/pscore.mwx * (1-pscore.mwx)/(1-pscore.wx)* (pscore.wx)/pscore.x ))
      y0m1p=sum( y*(1-d)/(1-pscore.mwx) * (pscore.mwx)/(pscore.wx)* (1-pscore.wx)/(1-pscore.x) )/sum((1-d)/(1-pscore.mwx) * (pscore.mwx)/(pscore.wx)* (1-pscore.wx)/(1-pscore.x) )
    }
    if (ATET==TRUE){
      y1m1=sum(y*d)/sum(d)
      y1m0<-(sum(y*d*pscore.x*(1-pscore.mwx)/((1-pscore.x)*pscore.mwx))/sum(d*pscore.x*(1-pscore.mwx)/((1-pscore.x)*pscore.mwx)))
      y0m0=(sum(y*(1-d)*pscore.x/(1-pscore.x))/sum((1-d)*pscore.x/(1-pscore.x)))
      y0m1<-(sum(y*(1-d)* pscore.mwx/((1-pscore.mwx)))/sum((1-d)* pscore.mwx/((1-pscore.mwx))))
      y1m0p=(sum( y*d*pscore.x/pscore.mwx * (1-pscore.mwx)/(1-pscore.wx)* (pscore.wx)/pscore.x )/sum(d*pscore.x/pscore.mwx * (1-pscore.mwx)/(1-pscore.wx)* (pscore.wx)/pscore.x ))
      y0m1p=sum( y*(1-d)*pscore.x/(1-pscore.mwx) * (pscore.mwx)/(pscore.wx)* (1-pscore.wx)/(1-pscore.x) )/sum((1-d)*pscore.x/(1-pscore.mwx) * (pscore.mwx)/(pscore.wx)* (1-pscore.wx)/(1-pscore.x) )
    }
    results=c(y1m1-y0m0, y1m1 - y0m1, y1m0-y0m0, y1m1 - y1m0p, y0m1p - y0m0, sum(ind))

  }
  results
}

bootstrap.mediation<-function(y,d,m,x,w=NULL,s=NULL,z=NULL,boot=1999, selpop=FALSE, trim=0.05, ATET=FALSE, logit=FALSE, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(y)
    bsamples=matrix(NA,boot,6)
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE)
      yb=y[sboot]
      db=d[sboot]
      if (is.null(s)==FALSE) sb=s[sboot]
      if (is.null(s)==TRUE) sb=NULL
      if (is.null(ncol(m))) mb<-m[sboot]
      if (is.null(ncol(m))==0) mb<-m[sboot,]
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      if ( (is.null(w)==FALSE) & (length(w)==length(y))) wb=w[sboot]
      if ( (is.null(w)==FALSE) & (length(w)!=length(y))) wb=w[sboot,]
      if (is.null(w)==TRUE) wb=NULL
      if ( (is.null(z)==FALSE) & (length(z)==length(y))) zb=z[sboot]
      if ( (is.null(z)==FALSE) & (length(z)!=length(y))) zb=z[sboot,]
      if (is.null(z)==TRUE) zb=NULL

      bsamples[i,]=c(mediation(y=yb,d=db,m=mb,x=xb,w=wb, s=sb, z=zb, selpop=selpop, trim=trim, ATET=ATET, logit=logit))
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; mb=c(); wb=c(); sb=c(); zb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]])
        if (is.null(s)==FALSE) sb<-c(sb,s[key==sboot[k]])
        if (is.null(ncol(m))) mb<-c(mb,m[key==sboot[k]])
        if (is.null(ncol(m))==0) mb=rbind(mb,m[key==sboot[k],])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
        if ((is.null(w)==FALSE) & is.null(ncol(w))) wb<-c(wb,w[key==sboot[k]])
        if ((is.null(w)==FALSE) & is.null(ncol(w))==0) wb=rbind(wb,w[key==sboot[k],])
        if ((is.null(z)==FALSE) & is.null(ncol(z))) zb<-c(zb,z[key==sboot[k]])
        if ((is.null(z)==FALSE) & is.null(ncol(z))==0) zb=rbind(zb,z[key==sboot[k],])
      }
      if (is.null(w)==TRUE) wb=NULL
      if (is.null(s)==TRUE) sb=NULL
      if (is.null(z)==TRUE) zb=NULL
      est=c(mediation(y=yb,d=db,m=mb,x=xb,w=wb, s=sb, z=zb, trim=trim, ATET=ATET, logit=logit))
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0) cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  bsamples
}


ipw<-function(y,d,x,s,z, selpop=FALSE, trim=0.05, ATET=FALSE, logit=FALSE){
   if (logit==FALSE) {
      if ( (is.null(s))  | (is.null(s)==0 & is.null(z)) )  pscore.x=glm(d~x,family=binomial(probit))$fitted
      if (is.null(s)==0 & is.null(z)) selscore=glm(s~cbind(d,x),family=binomial(probit))$fitted
      if (is.null(s)==0 & is.null(z)==0) {
        selscore=glm(s~cbind(d,x,z),family=binomial(probit))$fitted
        pscore.x=glm(d~cbind(x,selscore),family=binomial(probit))$fitted
      }
    }
    if (logit==TRUE)  {
      if ( (is.null(s))  | (is.null(s)==0 & is.null(z)) )  pscore.x=glm(d~x,family=binomial(logit))$fitted
      if (is.null(s)==0 & is.null(z)) selscore=glm(s~cbind(d,x),family=binomial(logit))$fitted
      if (is.null(s)==0 & is.null(z)==0) {
        selscore=glm(s~cbind(d,x,z),family=binomial(logit))$fitted
        pscore.x=glm(d~cbind(x,selscore),family=binomial(logit))$fitted
      }
    }
    if (ATET==FALSE){
      if (is.null(s)){
        ind=((pscore.x<trim) | (pscore.x>(1-trim)) )
        y=y[ind==0]; d=d[ind==0];  pscore.x=pscore.x[ind==0]
        y1=sum(y*d/pscore.x)/sum(d/pscore.x)
        y0=(sum(y*(1-d)/(1-pscore.x))/sum((1-d)/(1-pscore.x)))
      }
      if ((is.null(s)==0 & is.null(z)) | (is.null(s)==0 & is.null(z)==0 & selpop==FALSE))  {
        ind=((pscore.x<trim) | (pscore.x>(1-trim)) | selscore<trim )
        y=y[ind==0]; d=d[ind==0]; s=s[ind==0];  pscore.x=pscore.x[ind==0]; selscore=selscore[ind==0]
        y1=sum(y*d*s/(pscore.x*selscore))/sum(d*s/(pscore.x*selscore))
        y0=(sum(y*(1-d)*s/((1-pscore.x)*selscore))/sum((1-d)*s/((1-pscore.x)*selscore)))
      }
      if  (is.null(s)==0 & is.null(z)==0 & selpop==TRUE)  {
        ind=((pscore.x<trim) | (pscore.x>(1-trim)) )
        y=y[ind==0 & s==1]; d=d[ind==0 & s==1];  pscore.x=pscore.x[ind==0 & s==1]
        y1=sum(y*d/pscore.x)/sum(d/pscore.x)
        y0=(sum(y*(1-d)/(1-pscore.x))/sum((1-d)/(1-pscore.x)))
      }
    }
    if (ATET==TRUE){
      if (is.null(s)){
        ind= (pscore.x>(1-trim))
        y=y[ind==0]; d=d[ind==0];  pscore.x=pscore.x[ind==0]
        y1=sum(y*d)/sum(d)
        y0=(sum(y*(1-d)*pscore.x/(1-pscore.x))/sum((1-d)*pscore.x/(1-pscore.x)))
      }
      if ((is.null(s)==0 & is.null(z)) | (is.null(s)==0 & is.null(z)==0 & selpop==FALSE))  {
        ind= (pscore.x>(1-trim) | selscore<trim )
        y=y[ind==0]; d=d[ind==0]; s=s[ind==0]; pscore.x=pscore.x[ind==0]; selscore=selscore[ind==0]
        y1=sum(y*d*s/selscore)/sum(d*s/selscore)
        y0=(sum(y*(1-d)*s*pscore.x/((1-pscore.x)*selscore))/sum((1-d)*s*pscore.x/((1-pscore.x)*selscore)))
      }
      if  (is.null(s)==0 & is.null(z)==0 & selpop==TRUE)  {
        ind= (pscore.x>(1-trim))
        y=y[ind==0 & s==1]; d=d[ind==0 & s==1];  pscore.x=pscore.x[ind==0 & s==1]
        y1=sum(y*d)/sum(d)
        y0=(sum(y*(1-d)*pscore.x/(1-pscore.x))/sum((1-d)*pscore.x/(1-pscore.x)))
      }
    }
  results=c(y1-y0, y1, y0, sum(ind))
  results
}

bootstrap.ipw<-function(y,d,x,s=NULL,z=NULL, selpop=FALSE, boot=1999,trim=0.05, ATET=FALSE, logit=FALSE, cluster=NULL){
  if (is.null(cluster)){
   obs<-length(y)
    bsamples=matrix(NA,boot,4)
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE)
      yb=y[sboot]; db<-d[sboot]
      if (is.null(s)==0) sb=s[sboot]
      if (is.null(s)) sb=NULL
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      if (is.null(z)==0){
        if (is.null(ncol(z))) zb<-z[sboot]
        if (is.null(ncol(z))==0) zb<-z[sboot,]
      }
      if (is.null(z)) zb=NULL
      bsamples[i,]=c(ipw(y=yb,d=db,x=xb, s=sb, z=zb, selpop=selpop, trim=trim, ATET=ATET, logit=logit))
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; sb=c(); zb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]])
        if (is.null(s)==0) sb<-c(sb,s[key==sboot[k]])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
        if ((is.null(z)==FALSE) & is.null(ncol(z))) zb<-c(zb,z[key==sboot[k]])
        if ((is.null(z)==FALSE) & is.null(ncol(z))==0) zb=rbind(zb,z[key==sboot[k],])
      }
      if (is.null(s)) sb=NULL
      if (is.null(z)) zb=NULL
      est=c(ipw(y=yb,d=db,x=xb, s=sb, z=zb, selpop=selpop, trim=trim, ATET=ATET, logit=logit))
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0){
    cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  }
  bsamples
}



late<-function(y,d,z, x,trim=0.05, LATT=FALSE, logit=FALSE){
  if (logit==FALSE) pscore.x=glm(z~x,family=binomial(probit))$fitted
  if (logit==TRUE)  pscore.x=glm(z~x,family=binomial(logit))$fitted
  if (LATT==FALSE){
    ind=((pscore.x<trim) | (pscore.x>(1-trim)) )
    y=y[ind==0]; d=d[ind==0]; z=z[ind==0]; pscore.x=pscore.x[ind==0]
    firststage=sum(d*z/pscore.x)/(sum(z/pscore.x))-sum(d*(1-z)/(1-pscore.x))/(sum((1-z)/(1-pscore.x)))
    ITT=sum(y*z/pscore.x)/(sum(z/pscore.x))-sum(y*(1-z)/(1-pscore.x))/(sum((1-z)/(1-pscore.x)))
  }
  if (LATT==TRUE){
    ind= (pscore.x>(1-trim))
    y=y[ind==0]; d=d[ind==0]; z=z[ind==0]; pscore.x=pscore.x[ind==0]
    firststage=sum(d*z)/(sum(z))-sum(d*(1-z)*pscore.x/(1-pscore.x))/(sum((1-z)*pscore.x/(1-pscore.x)))
    ITT=sum(y*z)/(sum(z))-sum(y*(1-z)*pscore.x/(1-pscore.x))/(sum((1-z)*pscore.x/(1-pscore.x)))
  }
  results=c(ITT/firststage,  firststage, ITT, sum(ind))
  results
}


bootstrap.late<-function(y,d,z,x,boot=1999,trim=0.05, LATT=FALSE, logit=FALSE, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(y)
    bsamples=matrix(NA,boot,4)
    for(i in 1:boot){
      sboot=sample(1:obs,obs,TRUE)
      yb=y[sboot]; db=d[sboot]; zb=z[sboot];
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      bsamples[i,]=c(late(y=yb,d=db, z=zb, x=xb, trim=trim, LATT=LATT, logit=logit))
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; zb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]]); zb<-c(zb,z[key==sboot[k]])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
      }
      est=c(late(y=yb,d=db, z=zb, x=xb, trim=trim, LATT=LATT, logit=logit))
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0){
    cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  }
  bsamples
}


effects.late.x<-function(y,d,m,zd,  x, zm, trim=0.05, csquared=FALSE, bwreg=bwreg, bwm=bwm, cminobs=40, logit=FALSE){
  if (is.null(bwreg) | is.null(bwm)) temp<-npcdensbw(ydat=m, xdat=data.frame(zm,x), ckertype="gaussian", bwmethod="normal-reference")
  if (is.null(bwreg))  bwreg<-temp$xbw
  if (is.null(bwm))  bwm<-temp$ybw
  m.dist=npcdist(bws=c(bwm,bwreg), tydat=m, txdat=data.frame(zm,x), ckertype="gaussian")$condist
  x<-as.matrix(cbind(x))
  zm<-as.matrix(cbind(zm))
  nobs=length(y)
  dzd=d*zd
  one_dzd=(1-d)*(1-zd)
  if (logit==FALSE) {
    pscore2=glm(zd~x,family=binomial(probit))$fitted
    pred.d<-fitted.values(glm(d~zd+x,family=binomial(probit)))
  }
  if (logit==TRUE) {
    pscore2=glm(zd~x,family=binomial(logit))$fitted
    pred.d<-fitted.values(glm(d~zd+x,family=binomial(logit)))
  }
  pred.m<-fitted.values(lm(m~cbind(zm,pred.d,x)))
  templm<-lm(y~cbind(pred.d,pred.m,x))
  templm2<-lm(m~cbind(pred.d,x))$coef[2]

  c.d=d*(zd-pscore2)
  c.d_1=(d-1)*(zd-pscore2)
  c.den.d=lm(c.d~zm+x)$fitted

    mm<-sort(m)
     c.num<-c()
    for (i in 1:nrow(x)){
      ind= ((m<=m[i]) | (m<=mm[cminobs]))
      c.num.d=lm(c.d[ind==1]~zm[ind==1,]+x[ind==1,])
      c.num.d_1=lm(c.d_1[ind==1]~zm[ind==1,]+x[ind==1,])
      c.num<-c(c.num, (d[i]*(c(1,zm[i,],x[i,])%*%c.num.d$coef)+(1-d[i])*(c(1,zm[i,],x[i,])%*%c.num.d_1$coef)))
    }
  c=c.num/c.den.d*m.dist

  ind= (is.infinite(c)==0) & (is.na(c)==0)
  y=y[ind==1]; d=d[ind==1]; m=m[ind==1]; x=x[ind==1,]; zm=zm[ind==1,]; zd=zd[ind==1]; c=c[ind==1]; dzd=dzd[ind==1]; one_dzd=one_dzd[ind==1]; pscore2=pscore2[ind==1]


  regs=cbind(c,m,x)
  if (csquared==TRUE)  regs=cbind(regs, c^2)
    if (logit==FALSE){
      pscore1<-glm(zd~regs,family=binomial(probit))$fitted.values
      pscored<-glm(d~regs,family=binomial(probit))$fitted.values
      pscored[pscored<0]=0; pscored[pscored>1]=1
      pscoredz<-glm(dzd~regs,family=binomial(probit))$fitted.values
      pscoreone_dz<-glm(one_dzd~regs,family=binomial(probit))$fitted.values
      pscoreone_dz[pscoreone_dz<0]=0; pscoreone_dz[pscoreone_dz>1]=1
    }
    if (logit==TRUE){
      pscore1<-glm(zd~regs,family=binomial(logit))$fitted.values
      pscored<-glm(d~regs,family=binomial(logit))$fitted.values
      pscored[pscored<0]=0; pscored[pscored>1]=1
      pscoredz<-glm(dzd~regs,family=binomial(logit))$fitted.values
      pscoreone_dz<-glm(one_dzd~regs,family=binomial(logit))$fitted.values
      pscoreone_dz[pscoreone_dz<0]=0; pscoreone_dz[pscoreone_dz>1]=1
    }


  wgt=zd/pscore2/sum(zd/pscore2)-(1-zd)/(1-pscore2)/sum((1-zd)/(1-pscore2))
  firststage=(d*wgt)
  omega=1-(pscore1-pscore2)/(pscoredz-pscored*pscore2)
  one_omega=1/omega
  ind= (is.infinite(firststage)==0) & (is.infinite(wgt)==0) & (is.infinite(one_omega)==0) & (is.infinite(omega)==0)  &  (is.na(firststage)==0) &  (is.na(wgt)==0) & (is.na(one_omega)==0) & (is.na(omega)==0) & (d*omega*wgt/sum(d*omega*wgt)<=trim) & (d*wgt/sum(d*wgt)<=trim) & ((d-1)*one_omega*wgt/sum((d-1)*one_omega*wgt)<=trim) & ((d-1)*wgt/sum((d-1)*wgt)<=trim)
  y=y[ind==1]; d=d[ind==1]; omega=omega[ind==1]; wgt=wgt[ind==1]; firststage=firststage[ind==1]; one_omega=one_omega[ind==1]
  firststage=sum(firststage)
  y1m0=sum(y*d*omega*wgt)/firststage
  y1m1=sum( y*d*wgt )/firststage
  y0m1=sum(y*(d-1)*one_omega*wgt )/firststage
  y0m0= sum(y*(d-1)*wgt)/firststage
  results=c(y1m1-y0m0, y1m1-y0m1, y1m0-y0m0,  y1m1-y1m0, y0m1-y0m0, templm$coef[2], templm$coef[3]*templm2, nobs-length(y))
  results
}

bootstrap.mediation.late.x<-function(y,d,m,zd,zm,x, boot=1999,trim=0.05, csquared=FALSE, bwreg=bwreg, bwm= bwm, cminobs=40, logit=FALSE, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(y)
    bsamples=matrix(NA,boot,8)
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE)
      yb=y[sboot]; db<-d[sboot]; zdb<-zd[sboot]; mb=m[sboot]
      if (is.null(ncol(zm))) zmb<-zm[sboot]
      if (is.null(ncol(zm))==0) zmb<-zm[sboot,]
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      bsamples[i,]=effects.late.x(y=yb,d=db,m=mb,zd=zdb, zm=zmb, x=xb, trim=trim, csquared=csquared, bwreg=bwreg, bwm=bwm, cminobs=cminobs, logit=logit)
    }
  }

  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; zdb=c(); zmb=c(); mb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]]);
        zdb<-c(zdb,zd[key==sboot[k]]); mb<-c(mb,m[key==sboot[k]])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
        if (is.null(ncol(zm))) zmb<-c(zmb,zm[key==sboot[k]])
        if (is.null(ncol(zm))==0) zmb=rbind(zmb,zm[key==sboot[k],])
      }
      est=effects.late.x(y=yb,d=db,m=mb,zd=zdb, zm=zmb, x=xb, trim=trim, csquared=csquared, bwreg=bwreg, bwm=bwm, cminobs=cminobs, logit=logit)
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0){
    cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  }
  bsamples
}


mediation.cont<-function(y,d,m,x, d0, d1, ATET=FALSE, trim=0.05, lognorm=FALSE, bw){
  if(lognorm==TRUE){
    dd=d;dd[d==0]=0.00001
    ggg=glm(log(dd)~x)
    if (d0==0) d0=0.00001; if (d1==0) d1=0.00001
    pscore1d0=(dnorm( (log(d0)-cbind(1,x)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/d0)
    pscore1d1=(dnorm( (log(d1)-cbind(1,x)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/d1)
    ggg=glm(log(dd)~cbind(x,m))
    pscore2d0=(dnorm( (log(d0)-cbind(1,x,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/d0)
    pscore2d1=(dnorm( (log(d1)-cbind(1,x,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/d1)
  }
  if(lognorm==FALSE){
    ggg=glm(d~x)
    pscore1d0=(dnorm( (d0-cbind(1,x)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
    pscore1d1=(dnorm( (d1-cbind(1,x)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
    ggg=glm(d~cbind(x,m))
    pscore2d0=(dnorm( (d0-cbind(1,x,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
    pscore2d1=(dnorm( (d1-cbind(1,x,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
  }
  kernwgtd0=npksum(bws=bw, txdat = d, tydat = y, exdat = d0, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  kernwgtd1=npksum(bws=bw, txdat = d, tydat = y, exdat = d1, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  if (ATET==FALSE) ind= ((kernwgtd1/pscore1d1)/sum(kernwgtd1/pscore1d1)<=trim) & ((kernwgtd1*pscore2d0/(pscore2d1*pscore1d0))/sum(kernwgtd1*pscore2d0/(pscore2d1*pscore1d0))<=trim) & ((kernwgtd0/pscore1d0)/sum(kernwgtd0/pscore1d0)<=trim) & ((kernwgtd0*pscore2d1/(pscore2d0*pscore1d1))/sum(kernwgtd0*pscore2d1/(pscore2d0*pscore1d1))<=trim)
  if (ATET==TRUE)  ind= ((kernwgtd1)/sum(kernwgtd1)<=trim) & ((kernwgtd1*pscore2d0*pscore1d1/(pscore2d1*pscore1d0))/sum(kernwgtd1*pscore2d0*pscore1d1/(pscore2d1*pscore1d0))<=trim) & ((kernwgtd0*pscore1d1/pscore1d0)/sum(kernwgtd0*pscore1d1/pscore1d0)<=trim) & ((kernwgtd0*pscore2d1/pscore2d0)/sum(kernwgtd0*pscore2d1/pscore2d0)<=trim)
  y=y[ind];  pscore1d0=pscore1d0[ind];pscore1d1=pscore1d1[ind];pscore2d0=pscore2d0[ind]; pscore2d1=pscore2d1[ind];
  kernwgtd0=kernwgtd0[ind]; kernwgtd1= kernwgtd1[ind]
  if (ATET==FALSE){
    yd1m1=sum(y*kernwgtd1/pscore1d1)/sum(kernwgtd1/pscore1d1)
    yd1m0=sum(y*kernwgtd1*pscore2d0/(pscore2d1*pscore1d0))/sum(kernwgtd1*pscore2d0/(pscore2d1*pscore1d0))
    yd0m0=sum(y*kernwgtd0/pscore1d0)/sum(kernwgtd0/pscore1d0)
    yd0m1=sum(y*kernwgtd0*pscore2d1/(pscore2d0*pscore1d1))/sum(kernwgtd0*pscore2d1/(pscore2d0*pscore1d1))
  }
  if (ATET==TRUE){
    yd1m1=sum(y*kernwgtd1)/sum(kernwgtd1)
    yd1m0=sum(y*kernwgtd1*pscore2d0*pscore1d1/(pscore2d1*pscore1d0))/sum(kernwgtd1*pscore2d0*pscore1d1/(pscore2d1*pscore1d0))
    yd0m0=sum(y*kernwgtd0*pscore1d1/pscore1d0)/sum(kernwgtd0*pscore1d1/pscore1d0)
    yd0m1=sum(y*kernwgtd0*pscore2d1/pscore2d0)/sum(kernwgtd0*pscore2d1/pscore2d0)
  }
  results=c(yd1m1-yd0m0, yd1m1 - yd0m1, yd1m0-yd0m0, yd1m1 - yd1m0, yd0m1 - yd0m0, sum(1-ind))
}

bootstrap.mediation.cont<-function(y,d,m,x,d0,d1,ATET=FALSE, trim=0.05, lognorm=FALSE, bw, boot=1999, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(y)
    bsamples=matrix(NA,boot,6)
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE)
      yb=y[sboot]
      db<-d[sboot]
      if (is.null(ncol(m))) mb<-m[sboot]
      if (is.null(ncol(m))==0) mb<-m[sboot,]
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      bsamples[i,]=c(mediation.cont(y=yb,d=db,m=mb,x=xb, d0=d0, d1=d1,  ATET=ATET, trim=trim, lognorm=lognorm, bw=bw))
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; mb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]])
        if (is.null(ncol(m))) mb<-c(mb,m[key==sboot[k]])
        if (is.null(ncol(m))==0) mb=rbind(mb,m[key==sboot[k],])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
      }
      est=c(mediation.cont(y=yb,d=db,m=mb,x=xb, d0=d0, d1=d1,  ATET=ATET, trim=trim, lognorm=lognorm, bw=bw))
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0) cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  bsamples
}

attrlate<-function(y1,y2,r1,r2,d,z,x0,x1, weightmax=0.1){
  pz<-fitted.values(glm(z~x0,family=binomial(probit)))
  temp1<-glm(r1~cbind(x0,d),family=binomial(probit))
  pr1<-pnorm(cbind(1,x0,1)%*%(temp1$coef))
  temp2<-glm(r2~cbind(x0,x1,r1,d),family=binomial(probit))
  pr2<-pnorm(cbind(1,x0,x1,1,1)%*%(temp2$coef))
  temp3<-glm(z~cbind(x0,d),family=binomial(probit))
  pz0<-pnorm(cbind(1,x0,1)%*%(temp3$coef))
  temp4<-glm(z~cbind(x0,r1,d),family=binomial(probit))
  pz01<-pnorm(cbind(1,x0,1,1)%*%(temp4$coef))
  temp5<-glm(z~cbind(x0,x1,r1,d),family=binomial(probit))
  pz1<-pnorm(cbind(1,x0,x1,1,1)%*%(temp5$coef))
  temp6<-glm(z~cbind(x0,x1,r2,d),family=binomial(probit))
  pz12<-pnorm(cbind(1,x0,x1,1,1)%*%(temp6$coef))

  temp7<-glm(r1~cbind(x0,z,d),family=binomial(probit))
  pr1z1<-pnorm(cbind(1,x0,1,1)%*%(temp7$coef))
  pr1z0<-pnorm(cbind(1,x0,0,1)%*%(temp7$coef))

  temp8<-glm(r2~cbind(x0,x1,r1,z,d),family=binomial(probit))
  pr2z1<-pnorm(cbind(1,x0,x1,1,1,1)%*%(temp8$coef))
  pr2z0<-pnorm(cbind(1,x0,x1,1,0,1)%*%(temp8$coef))

  ggg=(mean(d/pz*(z-pz)/(1-pz)))
  weightz1<- (r1*d*z*1/(pz*pr1*((pz01-pz)/(pz0-pz))*ggg))
  relweightz1<-weightz1/sum(weightz1)
  indz1a<-1-(relweightz1>weightmax)

  weightz0<-(r1*d*(1-z)*1/((1-pz)*pr1*((pz01-pz)/(pz0-pz))*ggg))
  relweightz0<-weightz0/sum(weightz0)
  relweightz0[is.na(relweightz0)]=0
  indz0a<-1-(relweightz0>weightmax)

  y11<-(y1*r1*(d/pz)*(z-pz)/(1-pz)*(1/pr1)*((pz0-pz)/(pz01-pz)) )*(1/mean(d/pz*(z-pz)/(1-pz)))
  y11t<-mean(y11[(indz1a*indz0a)==1])

  indexli1=sum((indz1a*indz0a)==0)

  weightz1<- (r2*r1*d*z*1/(pz*pr1*((pz01-pz)/(pz0-pz))*pr2*((pz12-pz)/(pz1-pz))*ggg))
  relweightz1<-weightz1/sum(weightz1)
  indz1a<-1-(relweightz1>weightmax)

  weightz0<-(r2*r1*d*(1-z)*1/((1-pz)*pr1*((pz01-pz)/(pz0-pz))*pr2*((pz12-pz)/(pz1-pz))*ggg))
  relweightz0<-weightz0/sum(weightz0)
  relweightz0[is.na(relweightz0)]=0
  indz0a<-1-(relweightz0>weightmax)

  y12<-(y2*r1*r2*(d/pz)*(z-pz)/(1-pz)*(1/pr1)*((pz0-pz)/(pz01-pz)) *(1/pr2)*((pz1-pz)/(pz12-pz)) )*(1/mean(d/pz*(z-pz)/(1-pz)))
  y12t<-mean(y12[(indz1a*indz0a)==1])

  indexli2=sum((indz1a*indz0a)==0)

  weightz1<-((r1*d*z/pz*(1/pr1z1))/mean(d/pz*(z-pz)/(1-pz) ))
  relweightz1<-weightz1/sum(weightz1)
  indz1a<-1-(relweightz1>weightmax)

  weightz0<-((r1*d*(1-z)/(1-pz)*(1/pr1z0))/mean(d/pz*(z-pz)/(1-pz) ))
  relweightz0<-weightz0/sum(weightz0)
  relweightz0[is.na(relweightz0)]=0
  indz0a<-1-(relweightz0>weightmax)

  y11mar<-(y1*r1*d*z/pz*(1/pr1z1))/mean(d/pz*(z-pz)/(1-pz) )-(y1*r1*d*(1-z)/(1-pz)*(1/pr1z0))/mean(d/pz*(z-pz)/(1-pz))
  y11mart<-mean(y11mar[(indz1a*indz0a)==1])

  indexmar1=sum((indz1a*indz0a)==0)

  weightz1<-((r1*r2*d*z/pz*(1/pr1z1)*(1/pr2z1))/mean(d/pz*(z-pz)/(1-pz)))
  relweightz1<-weightz1/sum(weightz1)
  indz1a<-1-(relweightz1>weightmax)

  weightz0<-((r1*r2*d*((1-z)/(1-pz)*(1/pr1z0)*(1/pr2z0)))/mean(d/pz*(z-pz)/(1-pz)))
  relweightz0<-weightz0/sum(weightz0)
  relweightz0[is.na(relweightz0)]=0
  indz0a<-1-(relweightz0>weightmax)

  y12mar<-(y2*r1*r2*d*(z/pz*(1/pr1z1)*(1/pr2z1)-(1-z)/(1-pz)*(1/pr1z0)*(1/pr2z0)))/mean(d/pz*(z-pz)/(1-pz))
  y12mart<-mean(y12mar[(indz1a*indz0a)==1])

  indexmar2=sum((indz1a*indz0a)==0)

  pz<-fitted.values(glm(z~x0,family=binomial(probit)))
  pr1<-pnorm(cbind(1,x0,0)%*%(temp1$coef))
  pr2<-pnorm(cbind(1,x0,x1,1,0)%*%(temp2$coef))
  pz0<-pnorm(cbind(1,x0,0)%*%(temp3$coef))
  pz01<-pnorm(cbind(1,x0,1,0)%*%(temp4$coef))
  pz1<-pnorm(cbind(1,x0,x1,1,0)%*%(temp5$coef))
  pz12<-pnorm(cbind(1,x0,x1,1,0)%*%(temp6$coef))
  pr1z1<-pnorm(cbind(1,x0,1,0)%*%(temp7$coef))
  pr1z0<-pnorm(cbind(1,x0,0,0)%*%(temp7$coef))
  pr2z1<-pnorm(cbind(1,x0,x1,1,1,0)%*%(temp8$coef))
  pr2z0<-pnorm(cbind(1,x0,x1,1,0,0)%*%(temp8$coef))
  ggg=(mean((1-d)/pz*(z-pz)/(1-pz)))
  weightz1<- (-r1*(1-d)*z*1/(pz*pr1*((pz01-pz)/(pz0-pz))*ggg))
  relweightz1<-weightz1/sum(weightz1)
  indz1a<-1-(relweightz1>weightmax)

  weightz0<-(-r1*(1-d)*(1-z)*1/((1-pz)*pr1*((pz01-pz)/(pz0-pz))*ggg))
  relweightz0<-weightz0/sum(weightz0)
  relweightz0[is.na(relweightz0)]=0
  indz0a<-1-(relweightz0>weightmax)

  y01<-(y1*r1*((1-d)/pz)*(z-pz)/(1-pz)*(1/pr1)*((pz0-pz)/(pz01-pz)) )*(1/mean((1-d)/pz*(z-pz)/(1-pz)))
  y01t<-mean(y01[(indz1a*indz0a)==1])

  indexli1=indexli1+sum((indz1a*indz0a)==0)

  weightz1<-(-r2*r1*(1-d)*z*1/(pz*pr1*((pz01-pz)/(pz0-pz))*pr2*((pz12-pz)/(pz1-pz))*ggg))
  relweightz1<-weightz1/sum(weightz1)
  indz1a<-1-(relweightz1>weightmax)

  weightz0<-(-r2*r1*(1-d)*(1-z)*1/((1-pz)*pr1*((pz01-pz)/(pz0-pz))*pr2*((pz12-pz)/(pz1-pz))*ggg))
  relweightz0<-weightz0/sum(weightz0)
  relweightz0[is.na(relweightz0)]=0
  indz0a<-1-(relweightz0>weightmax)

  y02<-(y2*r1*r2*((1-d)/pz)*(z-pz)/(1-pz)*(1/pr1)*((pz0-pz)/(pz01-pz)) *(1/pr2)*((pz1-pz)/(pz12-pz)) )*(1/mean((1-d)/pz*(z-pz)/(1-pz)))
  y02t<-mean(y02[(indz1a*indz0a)==1])

  indexli2=indexli2+sum((indz1a*indz0a)==0)

  weightz1<-((-r1*(1-d)*(z/pz*(1/pr1z1)))/mean((1-d)/pz*(z-pz)/(1-pz)))
  relweightz1<-weightz1/sum(weightz1)
  indz1a<-1-(relweightz1>weightmax)

  weightz0<-((-r1*(1-d)*((1-z)/(1-pz)*(1/pr1z0)))/mean((1-d)/pz*(z-pz)/(1-pz)))
  relweightz0<-weightz0/sum(weightz0)
  relweightz0[is.na(relweightz0)]=0
  indz0a<-1-(relweightz0>weightmax)

  y01mar<-(y1*r1*(1-d)*(z/pz*(1/pr1z1)-(1-z)/(1-pz)*(1/pr1z0)))/mean((1-d)/pz*(z-pz)/(1-pz))
  y01mart<-mean(y01mar[(indz1a*indz0a)==1])

  indexmar1=indexmar1+sum((indz1a*indz0a)==0)

  weightz1<-((-r1*r2*(1-d)*(z/pz*(1/pr1z1)*(1/pr2z1)))/mean((1-d)/pz*(z-pz)/(1-pz)))
  relweightz1<-weightz1/sum(weightz1)
  indz1a<-1-(relweightz1>weightmax)

  weightz0<-((-r1*r2*(1-d)*((1-z)/(1-pz)*(1/pr1z0)*(1/pr2z0)))/mean((1-d)/pz*(z-pz)/(1-pz)))
  relweightz0<-weightz0/sum(weightz0)
  relweightz0[is.na(relweightz0)]=0
  indz0a<-1-(relweightz0>weightmax)

  y02mar<-(y2*r1*r2*(1-d)*(z/pz*(1/pr1z1)*(1/pr2z1)-(1-z)/(1-pz)*(1/pr1z0)*(1/pr2z0)))/mean((1-d)/pz*(z-pz)/(1-pz))
  y02mart<-mean(y02mar[(indz1a*indz0a)==1])

  indexmar2=indexmar2+sum((indz1a*indz0a)==0)

  late1t<-y11t-y01t
  late2t<-y12t-y02t

  latemar1t<-y11mart-y01mart
  latemar2t<-y12mart-y02mart

  results=c(latemar1t, latemar2t, late1t, late2t, indexmar1, indexmar2, indexli1, indexli2 )
}


bootstrap.attrlate<-function(y1,y2,r1,r2,d,z,x0,x1,weightmax=0.1, boot=1999, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(d)
    bsamples=matrix(NA,boot,8)
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE)
      y1b=y1[sboot]; y2b=y2[sboot]; r1b=r1[sboot]; r2b=r2[sboot]; db<-d[sboot]; zb<-z[sboot]
      if (is.null(ncol(x0))) x0b<-x0[sboot]
      if (is.null(ncol(x0))==0) x0b<-x0[sboot,]
      if (is.null(ncol(x1))) x1b<-x1[sboot]
      if (is.null(ncol(x1))==0) x1b<-x1[sboot,]
      bsamples[i,]=c(attrlate(y1=y1b,y2=y2b,r1=r1b,r2=r2b,d=db,z=zb,x0=x0b,x1=x1b, weightmax=weightmax))
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      y1b=c(); y2b=c(); r1b=c(); r2b=c(); db=c(); zb=c(); x0b=c(); x1b=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); y1b<-c(y1b,y1[key==sboot[k]]);zb<-c(zb,z[key==sboot[k]]); y2b<-c(y2b,y2[key==sboot[k]])
        r1b<-c(r1b,r1[key==sboot[k]]); r2b<-c(r2b,r2[key==sboot[k]])
        if (is.null(ncol(x0))) x0b<-c(x0b,x0[key==sboot[k]])
        if (is.null(ncol(x0))==0) x0b=rbind(x0b,x0[key==sboot[k],])
        if (is.null(ncol(x1))) x1b<-c(x1b,x1[key==sboot[k]])
        if (is.null(ncol(x1))==0) x1b=rbind(x1b,x1[key==sboot[k],])
      }
      est=c(attrlate(y1=y1b,y2=y2b,r1=r1b,r2=r2b,d=db,z=zb,x0=x0b,x1=x1b, weightmax=weightmax))
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0){
    cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  }
  bsamples
}


# DiD functions
ipwcount<-function(y,d,x=NULL,trim=0.02){
  if (is.null(x)==FALSE) pscore1=glm(d~x,family=binomial(logit))$fitted;
  if (is.null(x)==TRUE)  pscore1=rep(mean(d),length(d))
  ind=(pscore1>(1-trim))
  reweight=(sum(y[ind==0]*(1-d[ind==0])* pscore1[ind==0]/(1-pscore1[ind==0]))/sum((1-d[ind==0])* pscore1[ind==0]/(1-pscore1[ind==0])))
  list(est=reweight, dropped=sum(ind))
}

ipw.did<-function(y,d,t, x=NULL, trim=0.02){
  index=1*cbind(d, t, (d*t)+(1-d)*(1-t))
  treat=1*(d*t)
  means=matrix(NA,3,1)
  totaldropped=0
  for (j in 1:3){xxx=NULL; yyy=y[index[,j]==1]; ttreat=treat[index[,j]==1];
  if (is.null(x)==FALSE & is.null(ncol(x)) ) xxx=x[index[,j]==1]
  if (is.null(x)==FALSE &  is.null(ncol(x))==0 ) xxx=x[index[,j]==1,]
  temp=ipwcount(y=yyy,d=ttreat,x=xxx,trim=trim)
  means[j,1]=temp$est
  totaldropped=totaldropped+temp$dropped
  }
  c(mean(y[d==1 & t==1])-means[1,1]-(means[2,1]-means[3,1]), totaldropped)
}

bootstrap.did<-function(y,d,t, x=NULL,boot=1999,trim=0.05, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(y)
    bsamples=matrix(NA,boot,1)
    for(i in 1:boot){
      sboot=sample(1:obs,obs,TRUE)
      yb=y[sboot]; db=d[sboot]; tb=t[sboot];
      if (is.null(x)==0){
        if (is.null(ncol(x))) xb<-x[sboot]
        if (is.null(ncol(x))==0) xb<-x[sboot,]
      }
      if (is.null(x)==1) xb=NULL
      bsamples[i,1]=c(ipw.did(y=yb,d=db,t=tb, x=xb, trim=trim)[1])
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; tb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]]); tb<-c(tb,t[key==sboot[k]])
        if (is.null(x)==0){
          if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
          if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
        }
      }
      if (is.null(x)==1) xb=NULL
      est=c(ipw.did(y=yb,d=db,t=tb, x=xb, trim=trim)[1])
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0){
    cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  }
  bsamples
}



hdmed=function(y,d,m,x,k=3, trim=0.05, order=1, normalized=TRUE){
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)
  if (order>1) {x=Generate.Powers(cbind(x),lambda=order); x=as.matrix(x,nrow(x),ncol(x))}
  stepsize=ceiling((1/k)*length(d))
  set.seed(1); idx= sample(length(d), replace=FALSE)
  score=c(); selall=c()
  # crossfitting procedure that splits sample in training an testing data
  for (i in 1:k){
    tesample=idx[((i-1)*stepsize+1):(min((i)*stepsize,length(d)))]
    dtr=d[-tesample]; dte=d[tesample]; ytr=y[-tesample]; yte=y[tesample]; ytr1= ytr[dtr==1]; ytr0= ytr[dtr==0];
    mtr=m[-tesample]; mte=m[tesample]; mtr1=mtr[dtr==1]; mtr0=mtr[dtr==0];
    if (is.null(ncol(x)) | ncol(x)==1) {
      xtr=x[-tesample]; xte=x[tesample]; xtr1=xtr[dtr==1]; xtr0=xtr[dtr==0]; xtr11=xtr[dtr==1 & mtr==1]; xtr10=xtr[dtr==1 & mtr==0]; xtr01=xtr[dtr==0 & mtr==1]; xtr00=xtr[dtr==0 & mtr==0]
    }
    if (is.null(ncol(x))==0 & ncol(x)>1) {
      xtr=x[-tesample,]; xte=x[tesample,]; xtr1=xtr[dtr==1,]; xtr0=xtr[dtr==0,]; xtr11=xtr[dtr==1 & mtr==1,]; xtr10=xtr[dtr==1 & mtr==0,]; xtr01=xtr[dtr==0 & mtr==1,]; xtr00=xtr[dtr==0 & mtr==0,]
    }
      ytr11=ytr[dtr==1 & mtr==1]; ytr10=ytr[dtr==1 & mtr==0]; ytr01=ytr[dtr==0 & mtr==1]; ytr00=ytr[dtr==0 & mtr==0];
      # tr stands for first training data, te for test data, "1" and "0" for subsamples with treated and nontreated
      tr=data.frame(ytr,dtr,xtr,mtr);
      tr1=data.frame(ytr1,xtr1,mtr1); tr0=data.frame(ytr0,xtr0,mtr0);
      te=data.frame(yte,xte,mte);
      # predict Pr(M=1|D=1,X) in test data
      pm1=rlassologit(mtr1~xtr1)
      pm1te=predict(pm1, xte, type="response")
      # predict Pr(M=1|D=0,X) in test data
      pm0=rlassologit(mtr0~xtr0)
      pm0te=predict(pm0, xte, type="response")
      # predict Pr(D=1|X) in test data
      pd=rlassologit(dtr~xtr)
      pdte=predict(pd, xte, type="response")
      if (ybin!=1) {
        # predict E(Y| D=1, M=1, X) in test data
        eymx11=rlasso(ytr11~xtr11)
        eymx11te=predict(eymx11, xte)
        # predict E(Y| D=0, M=1, X) in test data
        eymx01=rlasso(ytr01~xtr01)
        eymx01te=predict(eymx01, xte)
        # predict E(Y| D=1, M=0, X) in test data
        eymx10=rlasso(ytr10~xtr10)
        eymx10te=predict(eymx10, xte)
        # predict E(Y| D=0, M=0, X) in test data
        eymx00=rlasso(ytr00~xtr00)
        eymx00te=predict(eymx00, xte)
        #  predict E(Y|D=1, X) in test data
        eyx1=rlasso(ytr1~xtr1)
        eyx1te=predict(eyx1, xte)
        #  predict E(Y|D=0, X) in test data
        eyx0=rlasso(ytr0~xtr0)
        eyx0te=predict(eyx0, xte)
      }
      if (ybin==1) {
        eymx11=rlassologit(ytr11~xtr11)
        eymx11te=predict(eymx11, xte, type="response")
        # predict E(Y| D=0, M=1, X) in test data
        eymx01=rlassologit(ytr01~xtr01)
        eymx01te=predict(eymx01, xte, type="response")
        # predict E(Y| D=1, M=0, X) in test data
        eymx10=rlassologit(ytr10~xtr10)
        eymx10te=predict(eymx10, xte, type="response")
        # predict E(Y| D=0, M=0, X) in test data
        eymx00=rlassologit(ytr00~xtr00)
        eymx00te=predict(eymx00, xte, type="response")
        #  predict E(Y|D=1, X) in test data
        eyx1=rlassologit(ytr1~xtr1)
        eyx1te=predict(eyx1, xte, type="response")
        #  predict E(Y|D=0, X) in test data
        eyx0=rlassologit(ytr0~xtr0)
        eyx0te=predict(eyx0, xte, type="response")
      }
      # predict E(Y| D=0, M, X) in test data
      eymx0te=mte*eymx01te+(1-mte)*eymx00te
      # predict E(Y| D=1, M, X) in test data
      eymx1te=mte*eymx11te+(1-mte)*eymx10te

    # predict score functions for E(Y(1,M(0))) in the test data
    eta10=(eymx11te*pm0te+eymx10te*(1-pm0te))
    eta01=(eymx01te*pm1te+eymx00te*(1-pm1te))
    sel= 1*(((pdte*pm1te)>=trim) & ((1-pdte)>=trim)  & (pdte>=trim) &  (((1-pdte)*pm0te)>=trim)   )
    score=rbind(score, cbind(dte, pm0te, pdte, pm1te, yte, eymx1te, eta10, eyx1te, eymx0te, eta01, eyx0te))[sel==1,]
    selall=c(selall,sel)
  }
  # compute scores
  if (normalized==FALSE){
    y1m0=score[,1]*score[,2]/(score[,3]*score[,4])*(score[,5]-score[,6])+(1-score[,1])/(1-score[,3])*(score[,6]-score[,7] )+score[,7]
    y1m1=score[,8] + score[,1]*(score[,5]-score[,8])/score[,3]
    y0m1=(1-score[,1])*score[,4]/((1-score[,3])*score[,2])*(score[,5]-score[,9])+score[,1]/score[,3]*(score[,9]-score[,10])+score[,10]
    y0m0=score[,11] + (1-score[,1])*(score[,5]-score[,11])/(1-score[,3])
  }
  if (normalized!=FALSE){
    nobs=nrow(score)
    sumscores1=sum(score[,1]*score[,2]/(score[,3]*score[,4]))
    sumscores2=sum((1-score[,1])/(1-score[,3]))
    sumscores3=sum(score[,1]/score[,3])
    sumscores4=sum((1-score[,1])*score[,4]/((1-score[,3])*score[,2]))
    y1m0=(nobs*score[,1]*score[,2]/(score[,3]*score[,4])*(score[,5]-score[,6]))/sumscores1+(nobs*(1-score[,1])/(1-score[,3])*(score[,6]-score[,7]))/sumscores2+score[,7]
    y1m1=score[,8] + (nobs*score[,1]*(score[,5]-score[,8])/score[,3])/sumscores3
    y0m1=(nobs*(1-score[,1])*score[,4]/((1-score[,3])*score[,2])*(score[,5]-score[,9]))/sumscores4+(nobs*score[,1]/score[,3]*(score[,9]-score[,10]))/sumscores3+score[,10]
    y0m0=score[,11] + (nobs*(1-score[,1])*(score[,5]-score[,11])/(1-score[,3]))/sumscores2
  }

  # compute mean potential outcomes
  my1m1=mean(y1m1); my0m1=mean(y0m1); my1m0=mean(y1m0); my0m0=mean(y0m0)
  # compute effects
  tot=my1m1-my0m0; dir1=my1m1-my0m1; dir0=my1m0-my0m0; indir1=my1m1-my1m0; indir0=my0m1-my0m0;
  #compute variances
  vtot=mean((y1m1-y0m0-tot)^2); vdir1=mean((y1m1-y0m1-dir1)^2); vdir0=mean((y1m0-y0m0-dir0)^2);
  vindir1=mean((y1m1-y1m0-indir1)^2); vindir0=mean((y0m1-y0m0-indir0)^2); vcontrol=mean((y0m0-my0m0)^2)
  c(tot, dir1, dir0, indir1, indir0, my0m0, vtot, vdir1, vdir0, vindir1, vindir0, vcontrol, sum(selall))
}

# function for mediation with high dimensional covariates based on Bayes rule
hdmedalt=function(y,d,m,x, trim=0.05, order=1, fewsplits=FALSE, normalized=TRUE){
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)
  #generate higher order terms for lasso
    xm=cbind(x,m)
    if (order>1) {x=Generate.Powers(cbind(x),lambda=order); x=as.matrix(x,nrow(x),ncol(x))}
    if (order>1) xm=Generate.Powers(xm,lambda=order); xm=as.matrix(xm,nrow(xm),ncol(xm))
    stepsize=ceiling((1/3)*length(d))
    nobs = min(3*stepsize,length(d)); set.seed(1); idx = sample(nobs);
    sample1 = idx[1:stepsize]; sample2 = idx[(stepsize+1):(2*stepsize)];
    sample3 = idx[(2*stepsize+1):nobs];
    #nobs = min(4*stepsize,length(d)); set.seed(1); idx = sample(nobs);
    #sample1 = idx[1:stepsize]; sample2 = idx[(stepsize+1):(2*stepsize)];
    #sample3 = idx[(2*stepsize+1):(3*stepsize)]; sample4 = idx[(3*stepsize+1):nobs]
    score=c(); selall=c()
  # crossfitting procedure that splits sample in training an testing data
  for (i in 1:3){
    if (i==1) {tesample=sample1; musample=sample2; deltasample=sample3}
    if (i==2) {tesample=sample3; musample=sample1; deltasample=sample2}
    if (i==3) {tesample=sample2; musample=sample3; deltasample=sample1}
    trsample=c(musample,deltasample); dte=d[tesample]; yte=y[tesample]
    # in case that fewsplits is one, psample and musample are merged
    if (fewsplits==1){musample=c(musample,deltasample);deltasample=musample}
      x=as.matrix(x,nrow(x),ncol(x)); xm=as.matrix(xm,nrow(xm),ncol(xm))
      # fit Pr(D=1|M,X) in total of training data
      pmx=rlassologit(d[trsample]~xm[trsample,])
      # predict Pr(D=1|M,X) in test data
      pmxte=predict(pmx, xm[tesample,], type="response")
      # fit Pr(D=1|X) in total of training data
      px=rlassologit(d[trsample]~x[trsample,])
      # predict Pr(D=1|X) in test data
      pxte=predict(px, x[tesample,], type="response")
      # fit E(Y|M,X,D=1) in first training data
      if (ybin!=1){
        eymx1=rlasso(y[musample[d[musample]==1]]~xm[musample[d[musample]==1],])
        # predict E(Y|M,X,D=1) in test data
        eymx1te=predict(eymx1, xm[tesample,])
        # predict E(Y|M,X,D=1) in delta sample
        eymx1trte=predict(eymx1, xm[deltasample,])
      }
      if (ybin==1){
       eymx1=rlassologit(y[musample[d[musample]==1]]~xm[musample[d[musample]==1],])
        # predict E(Y|M,X,D=1) in test data
        eymx1te=predict(eymx1, xm[tesample,], type="response")
        # predict E(Y|M,X,D=1) in delta sample
        eymx1trte=predict(eymx1, xm[deltasample,], type="response")
      }
      # fit E[E(Y|M,X,D=1)|D=0,X] in delta sample
      dtrte=d[deltasample]; xtrte=x[deltasample,]
      regweymx1=rlasso(eymx1trte[dtrte==0]~xtrte[dtrte==0,])
      # predict E[E(Y|M,X,D=1)|D=0,X] in the test data
      regweymx1te=predict(regweymx1, x[tesample,])
      #  fit E(Y|X,D=1) in total of training data with D=1 by running Y~X
      if (ybin!=1){
        temp=rlasso(y[trsample[d[trsample]==1]]~x[trsample[d[trsample]==1],])
        # predict E(Y|X,D=1) in the test data
        eyx1te=predict(temp, x[tesample,])
        # fit E(Y|M,X,D=0) in first training data
        eymx0=rlasso(y[musample[d[musample]==0]]~xm[musample[d[musample]==0],])
        # predict E(Y|M,X,D=0) in test data
        eymx0te=predict(eymx0, xm[tesample,])
        # predict E(Y|M,X,D=0) in delta sample
        eymx0trte=predict(eymx0, xm[deltasample,])
      }
      if (ybin==1){
        temp=rlassologit(y[trsample[d[trsample]==1]]~x[trsample[d[trsample]==1],])
        # predict E(Y|X,D=1) in the test data
        eyx1te=predict(temp, x[tesample,], type="response")
        # fit E(Y|M,X,D=0) in first training data
        eymx0=rlassologit(y[musample[d[musample]==0]]~xm[musample[d[musample]==0],])
        # predict E(Y|M,X,D=0) in test data
        eymx0te=predict(eymx0, xm[tesample,], type="response")
        # predict E(Y|M,X,D=0) in delta sample
        eymx0trte=predict(eymx0, xm[deltasample,], type="response")
      }
      # fit E[E(Y|M,X,D=0)|D=1,X] in delta sample
      regweymx0=rlasso(eymx0trte[dtrte==1]~xtrte[dtrte==1,])
      # predict E[E(Y|M,X,D=0)|D=1, X] in the test data
      regweymx0te=predict(regweymx0, x[tesample,])
      if (ybin!=1){
        #  fit E(Y|X,D=0) in total of training data with D=0 by running Y~X
        temp=rlasso(y[trsample[d[trsample]==0]]~x[trsample[d[trsample]==0],])
        # predict E(Y|X,D=0) in the test data
        eyx0te=predict(temp, x[tesample,])
      }
      if (ybin==1){
        #  fit E(Y|X,D=0) in total of training data with D=0 by running Y~X
        temp=rlassologit(y[trsample[d[trsample]==0]]~x[trsample[d[trsample]==0],])
        # predict E(Y|X,D=0) in the test data
        eyx0te=predict(temp, x[tesample,], type="response")
      }
    # select observations satisfying trimming restriction
    sel= 1*((((1-pmxte)*pxte)>=trim) & ((1-pxte)>=trim)  & (pxte>=trim) &  (((pmxte*(1-pxte)))>=trim)   )

    #select elements of the score functions
    score=rbind(score, cbind(dte, pmxte, pxte, yte, eymx0te, regweymx0te, eyx0te, eymx1te, regweymx1te, eyx1te)[sel==1,])

    # collect selection dummies
    selall=c(selall,sel)
  }

  # compute scores for potential outcomes
  if (normalized==FALSE) {
    y0m1=((1-score[,1])*score[,2]/((1-score[,2])*score[,3])*(score[,4]-score[,5])+score[,1]/score[,3]*(score[,5]-score[,6])+score[,6])
    y0m0=(score[,7] + (1-score[,1])*(score[,4]-score[,7])/(1-score[,3]))
    y1m0=(score[,1]*(1-score[,2])/(score[,2]*(1-score[,3]))*(score[,4]-score[,8])+(1-score[,1])/(1-score[,3])*(score[,8]-score[,9])+score[,9])
    y1m1=(score[,10] + score[,1]*(score[,4]-score[,10])/score[,3])
  }
  if (normalized!=FALSE) {
    nobs=nrow(score)
    sumscore1=sum((1-score[,1])*score[,2]/((1-score[,2])*score[,3]))
    sumscore2=sum(score[,1]/score[,3])
    sumscore3=sum((1-score[,1])/(1-score[,3]))
    sumscore4=sum(score[,1]*(1-score[,2])/(score[,2]*(1-score[,3])))
    y0m1=(nobs*(1-score[,1])*score[,2]/((1-score[,2])*score[,3])*(score[,4]-score[,5]))/sumscore1+(nobs*score[,1]/score[,3]*(score[,5]-score[,6]))/sumscore2+score[,6]
    y0m0=score[,7] + (nobs*(1-score[,1])*(score[,4]-score[,7])/(1-score[,3]))/sumscore3
    y1m0=(nobs*score[,1]*(1-score[,2])/(score[,2]*(1-score[,3]))*(score[,4]-score[,8]))/sumscore4+(nobs*(1-score[,1])/(1-score[,3])*(score[,8]-score[,9]))/sumscore3+score[,9]
    y1m1=score[,10] + (nobs*score[,1]*(score[,4]-score[,10])/score[,3])/sumscore2
  }
  # compute mean potential outcomes
  my1m1=mean(y1m1); my0m1=mean(y0m1); my1m0=mean(y1m0); my0m0=mean(y0m0)
  # compute effects
  tot=my1m1-my0m0; dir1=my1m1-my0m1; dir0=my1m0-my0m0; indir1=my1m1-my1m0; indir0=my0m1-my0m0;
  #compute variances
  vtot=mean((y1m1-y0m0-tot)^2); vdir1=mean((y1m1-y0m1-dir1)^2); vdir0=mean((y1m0-y0m0-dir0)^2);
  vindir1=mean((y1m1-y1m0-indir1)^2); vindir0=mean((y0m1-y0m0-indir0)^2); vcontrol=mean((y0m0-my0m0)^2)
  # report effects, mean of Y(0,M(0)), variances, number of non-trimmed observations
  c(tot, dir1, dir0, indir1, indir0, my0m0, vtot, vdir1, vdir0, vindir1, vindir0, vcontrol, sum(selall))
}


# DYNAMIC TREATMENT EFFECTS WITH DOUBLE MACHINE LEARNING

MLfunct=function(y, x, d1=NULL, d2=NULL, MLmethod="lasso",  ybin=0){
  if (is.null(d1)==0 & is.null(d2)==0) { y=y[d1==1 & d2==1]; x=x[d1==1 & d2==1,]}
  if (is.null(d1)==0 & is.null(d2)==1) { y=y[d1==1]; x=x[d1==1,]}
  if  (MLmethod=="lasso"){
    if (ybin==1) model=SuperLearner(Y = y, X = x, family = binomial(), SL.library = "SL.glmnet")
    if (ybin!=1) model=SuperLearner(Y = y, X = x, family = gaussian(), SL.library = "SL.glmnet")
  }
  if  (MLmethod=="randomforest"){
    if (ybin==1) model=SuperLearner(Y = y, X = x, family = binomial(), SL.library = "SL.ranger")
    if (ybin!=1) model=SuperLearner(Y = y, X = x, family = gaussian(), SL.library = "SL.ranger")
  }
  if  (MLmethod=="xgboost"){
    if (ybin==1) model=SuperLearner(Y = y, X = x, family = binomial(), SL.library = "SL.xgboost")
    if (ybin!=1) model=SuperLearner(Y = y, X = x, family = gaussian(), SL.library = "SL.xgboost")
  }
 if  (MLmethod=="svm"){
   if (ybin==1) model=SuperLearner(Y = y, X = x, family = binomial(), SL.library = "SL.svm")
   if (ybin!=1) model=SuperLearner(Y = y, X = x, family = gaussian(), SL.library = "SL.svm")
 }
  if  (MLmethod=="ensemble"){
    if (ybin==1) model=SuperLearner(Y = y, X = x, family = binomial(), SL.library = c("SL.glmnet", "SL.xgboost", "SL.svm", "SL.ranger"))
    if (ybin!=1) model=SuperLearner(Y = y, X = x, family = gaussian(), SL.library = c("SL.glmnet", "SL.xgboost", "SL.svm", "SL.ranger"))
  }
  if  (MLmethod=="parametric"){
    if (ybin==1) model=SuperLearner(Y = y, X = x, family = binomial(), SL.library = "SL.glm")
    if (ybin!=1) model=SuperLearner(Y = y, X = x, family = gaussian(), SL.library = "SL.lm")
  }
  model
}

hddyntreat=function(y2,d1,d2,x0,x1, s=NULL, trim=0.05, MLmethod="lasso", fewsplits=fewsplits){
  ybin=1*(length(unique(y2))==2 & min(y2)==0 & max(y2)==1)  # check if binary outcome
  x0=data.frame(x0); x0x1=data.frame(x0,x1); d1x0x1=data.frame(d1,x0x1);
  # crossfitting procedure that splits sample in training an testing data
  stepsize=ceiling((1/3)*length(y2))
  nobs = min(3*stepsize,length(y2)); set.seed(1); idx = sample(nobs);
  sample1 = idx[1:stepsize]; sample2 = idx[(stepsize+1):(2*stepsize)]; sample3 = idx[(2*stepsize+1):nobs];
  score=c(); sel=c(); trimmed=c()
  for (i in 1:3){
    if (i==3) {trsample1=sample1; trsample2=sample2; tesample=sample3}
    if (i==1) {trsample1=sample2; trsample2=sample3; tesample=sample1}
    if (i==2) {trsample1=sample3; trsample2=sample1; tesample=sample2}
    # total training sample
    trsample=c(trsample1,trsample2)
    # in case that fewsplits is one, both training data are merged
    if (fewsplits==1){trsample1=c(trsample1,trsample2);trsample2=trsample1}
    if (is.null(s)) {gte=rep(1,length(tesample)); ste=gte} #check if weighted estimation should be performed
    if (is.null(s)==0) {
      g=MLfunct(y=s[trsample], x=x0[trsample,], MLmethod=MLmethod,  ybin=1)
      gte=predict(g, x0[tesample,], onlySL = TRUE)$pred     #predict weighting function in test data
      ste=s[tesample]
    }
    p1=MLfunct(y=d1[trsample], x=x0[trsample,], MLmethod=MLmethod,  ybin=1)
    p1te=predict(p1, x0[tesample,], onlySL = TRUE)$pred     #predict ps1 in test data

    p2=MLfunct(y=d2[trsample], x=d1x0x1[trsample,], MLmethod=MLmethod, ybin=1)
    p2te=predict(p2, d1x0x1[tesample,], onlySL = TRUE)$pred  #predict ps2 in test data

    y2d1d2=MLfunct(y=y2[trsample1], x=x0x1[trsample1,], d1=d1[trsample1], d2=d2[trsample1], MLmethod=MLmethod, ybin=ybin)
    y2d1d2te=predict(y2d1d2,  x0x1[tesample,], onlySL = TRUE)$pred  #predict E[Y2|D1,D2,X0,X1] in test data
    y2d1d2tr2=predict(y2d1d2, x0x1[trsample2,], onlySL = TRUE)$pred  #predict E[Y2|D1,D2,X0,X1] in second training data

    y1d1=MLfunct(y=y2d1d2tr2, x=x0[trsample2,], d1=d1[trsample2], MLmethod=MLmethod, ybin=0)
    y1d1te=predict(y1d1, x0[tesample,], onlySL = TRUE)$pred  #predict E[Y2|D1,D2,X0,X1] in test data

    # observations not satisfying trimming restriction
    trimmed=1*((p1te*p2te)<trim)
    score=rbind(score, cbind(gte,d1[tesample],d2[tesample],y2[tesample],y2d1d2te,p1te,p2te,y1d1te,ste,trimmed))
  }
  score = score[order(idx),]
  score
}


# ATE ESTIMATION BASED ON DML
hdtreat=function(y,d,x,s=NULL, trim=0.01, MLmethod="lasso", k=3){
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  x=data.frame(x)
  stepsize=ceiling((1/k)*length(d))
  set.seed(1); idx= sample(length(d), replace=FALSE)
  score=c();
  # crossfitting procedure that splits sample in training an testing data
  for (i in 1:k){
    tesample=idx[((i-1)*stepsize+1):(min((i)*stepsize,length(d)))]
    trsample=idx[-tesample]
    if (is.null(s)) {gte=rep(1,length(tesample)); ste=gte} #check if weighted estimation should be performed
    if (is.null(s)==0) {
      g=MLfunct(y=s[trsample], x=x[trsample,], MLmethod=MLmethod,  ybin=1)
      gte=predict(g, x[tesample,], onlySL = TRUE)$pred     #predict weighting function in test data
      ste=s[tesample]
    }
    ps=MLfunct(y=d[trsample], x=x[trsample,], MLmethod=MLmethod,  ybin=1)
    pste=predict(ps, x[tesample,], onlySL = TRUE)$pred     #predict propensity score in test data
    eydx=MLfunct(y=y[trsample], x=x[trsample,], d1=d[trsample], MLmethod=MLmethod, ybin=ybin)
    eydxte=predict(eydx, x[tesample,], onlySL = TRUE)$pred  #predict conditional outcome in test data
    # observations not satisfying trimming restriction
    trimmed=1*((pste)<trim)
    score=rbind(score, cbind(gte,d[tesample],y[tesample],eydxte,pste, ste,trimmed))
  }
  score = score[order(idx),]
  score
}


# ATE ESTIMATION WITH SAMPLE SELECTION BASED ON DML
hdseltreat=function(y,d,x,s, z, trim=0.01, MLmethod="lasso", k=3, selected=0){
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  x=data.frame(x)
  dx=data.frame(d,x)
  score=c();
  # crossfitting procedure that splits sample in training an testing data
  stepsize=ceiling((1/3)*length(y))
    nobs = min(3*stepsize,length(y)); set.seed(1); idx = sample(nobs);
    sample1 = idx[1:stepsize]; sample2 = idx[(stepsize+1):(2*stepsize)]; sample3 = idx[(2*stepsize+1):nobs];
    for (i in 1:3){
      if (i==3) {trsample1=sample1; trsample2=sample2; tesample=sample3}
      if (i==1) {trsample1=sample2; trsample2=sample3; tesample=sample1}
      if (i==2) {trsample1=sample3; trsample2=sample1; tesample=sample2}
      # total training sample
      trsample=c(trsample1,trsample2)
      ytr=y[trsample]; xtr=x[trsample,]; dtr=d[trsample]; str=s[trsample]; xte=x[tesample,]
    if (is.null(z)) {
      g=MLfunct(y=str, x=dx[trsample,], MLmethod=MLmethod,  ybin=1)
      gte=predict(g, dx[tesample,], onlySL = TRUE)$pred  # predict selection model under MAR
      ps=MLfunct(y=dtr, x=xtr, MLmethod=MLmethod,  ybin=1)
      pste=predict(ps, xte, onlySL = TRUE)$pred     #predict propensity score in test data
      eydx=MLfunct(y=ytr[str==1], x=xtr[str==1,], d1=dtr[str==1], MLmethod=MLmethod, ybin=ybin)
      eydxte=predict(eydx, xte, onlySL = TRUE)$pred  #predict conditional outcome in test data
     }
    if (is.null(z)==0) {
      dxz=data.frame(dx,z)
      g=MLfunct(y=s[trsample1], x=dxz[trsample1,], MLmethod=MLmethod,  ybin=1)
      gtotal=predict(g, dxz, onlySL = TRUE)$pred # predict selection model based on instrument
      gte=gtotal[tesample]
      xg=data.frame(x,gtotal)
      xgtr=xg[trsample2,]; xgte=xg[tesample,]
      ytr2=y[trsample2]; dtr2=d[trsample2]; str2=s[trsample2];
      ps=MLfunct(y=dtr2, x=xgtr, MLmethod=MLmethod,  ybin=1)
      pste=predict(ps, xgte, onlySL = TRUE)$pred     #predict propensity score in test data
      eydx=MLfunct(y=ytr2[str2==1], x=xgtr[str2==1,], d1=dtr2[str2==1], MLmethod=MLmethod, ybin=ybin)
      eydxte=predict(eydx, xgte, onlySL = TRUE)$pred  #predict conditional outcome in test data
      }
    # observations not satisfying trimming restriction
    if (selected!=1) trimmed=1*((pste*gte)<trim)
    if (selected==1) trimmed=1*(pste<trim)
    score=rbind(score, cbind(gte,d[tesample],y[tesample],eydxte,pste, s[tesample],trimmed))
  }
  score = score[order(idx),]
  score
}



# LATE WITH ATTRITION


latenonresp=function(y,d,r,z1,z2, bw1=NULL, bw2=NULL, bw3=NULL, bw4=NULL, bw5=NULL, bw6=NULL, bw7=NULL, bw8=NULL, bw9=NULL, bw10=NULL, bw11=NULL, bw12=NULL, ruleofthumb=1,  wgtfct=2, rtype="ll", numresprob=100, trim=0.01){

  Pz1=mean(z1); yphi=r*(z1-Pz1); yphi2=(z1-Pz1); Pco=mean(d[z1==1])-mean(d[z1==0]); n=length(r)

  if (is.null(bw1)) {
    if (ruleofthumb!=1) bw1<-npregbw(ydat=yphi[d==1], xdat=z2[d==1], regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw1<-npudensbw(dat=z2[d==1], ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  phi1a<-(npreg(bws=bw1, tydat=yphi[d==1], txdat=z2[d==1], exdat=z2, regtype=rtype, ckertype="gaussian")$mean)
  phi1apar<-cbind(1,z2)%*%coef(lm(yphi[d==1]~z2[d==1]))

  if (is.null(bw2)) {
    if (ruleofthumb!=1) bw2<-npregbw(ydat=yphi2[d==1], xdat=z2[d==1], regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw2<-bw1
  }
  phi1b<-(npreg(bws=bw2, tydat=yphi2[d==1], txdat=z2[d==1], exdat=z2, regtype=rtype, ckertype="gaussian")$mean)
  phi1bpar<-cbind(1,z2)%*%coef(lm(yphi2[d==1]~z2[d==1]))
  phi1=phi1a/phi1b; phi1par=phi1apar/phi1bpar

  if (is.null(bw3)) {
    if (ruleofthumb!=1) bw3<-npregbw(ydat=yphi[d==0], xdat=z2[d==0], regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw3<-npudensbw(dat=z2[d==0], ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  phi0a<-(npreg(bws=bw3, tydat=yphi[d==0], txdat=z2[d==0], exdat=z2, regtype=rtype, ckertype="gaussian")$mean)
  phi0apar<-cbind(1,z2)%*%coef(lm(yphi[d==0]~z2[d==0]))

  if (is.null(bw4)) {
    if (ruleofthumb!=1) bw4<-npregbw(ydat=yphi2[d==0], xdat=z2[d==0], regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw4<-bw3
  }
  phi0b<-(npreg(bws=bw4, tydat=yphi2[d==0], txdat=z2[d==0], exdat=z2, regtype=rtype, ckertype="gaussian")$mean)
  phi0bpar<-cbind(1,z2)%*%coef(lm(yphi2[d==0]~z2[d==0]))

  phi0=phi0a/phi0b; phi0par=phi0apar/phi0bpar

  yrd=y*r*d*(z1-Pz1);yr1_d=y*r*(1-d)*(z1-Pz1)
  minselprob=max(min(phi1),min(phi0),0.01); maxselprob=min(max(phi1),max(phi0),1)

  if (is.null(bw5)) {
    if (ruleofthumb!=1) bw5<-npregbw(ydat=yrd, xdat=phi1, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw5<-npudensbw(dat=phi1, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  predval=seq(from=minselprob,to=maxselprob,length.out=numresprob)
  eyrd<-(npreg(bws=bw5, tydat=yrd, txdat=phi1, exdat=predval, regtype=rtype, ckertype="gaussian")$mean)
  eyrdpar<-cbind(1,predval)%*%lm(yrd~phi1par)$coef

  if (is.null(bw6)) {
    if (ruleofthumb!=1) bw6<-npregbw(ydat=yr1_d, xdat=phi0, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw6<-npudensbw(dat=phi0, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  eyr1_d<-(npreg(bws=bw6, tydat=yr1_d, txdat=phi0, exdat=predval, regtype=rtype, ckertype="gaussian")$mean)
  eyr1_dpar<-cbind(1,predval)%*%lm(yr1_d~phi0par)$coef

  late=(eyrd+eyr1_d)*1/(Pz1*(1-Pz1)*predval*Pco)
  latepar=(eyrdpar+eyr1_dpar)*1/(Pz1*(1-Pz1)*predval*Pco)


  if (is.null(bw7)) {
    if (ruleofthumb!=1) bw7<-npudensbw(dat=phi1[d==1], ckertype="gaussian", bwmethod="cv.ls")
    if (ruleofthumb==1) bw7<-npudensbw(dat=phi1[d==1],ckertype="gaussian", bwmethod="normal-reference")
  }


  if (is.null(bw8)) {
    if (ruleofthumb!=1) bw8<-npudensbw(dat=phi1par[d==1], ckertype="gaussian", bwmethod="cv.ls")
    if (ruleofthumb==1) bw8<-npudensbw(dat=phi1par[d==1], ckertype="gaussian", bwmethod="normal-reference")
  }

  if (is.null(bw9)) {
    if (ruleofthumb!=1) bw9<-npudensbw(dat=phi0[d==0], ckertype="gaussian", bwmethod="cv.ls")
    if (ruleofthumb==1) bw9<-npudensbw(dat=phi0[d==0], ckertype="gaussian", bwmethod="normal-reference")
  }

  if (is.null(bw10)) {
    if (ruleofthumb!=1) bw10<-npudensbw(dat=phi0par[d==0], ckertype="gaussian", bwmethod="cv.ls")
    if (ruleofthumb==1) bw10<-npudensbw(dat=phi0par[d==0], ckertype="gaussian", bwmethod="normal-reference")
  }

  grid=c(predval); lgrid=length(grid)
  fphi1para=npudens(bws=bw6, tdat=phi1par[d==1], edat=grid)$dens
  fphi0para=npudens(bws=bw8, tdat=phi0par[d==0], edat=grid)$dens
  fphi1a=npudens(bws=bw5, tdat=phi1[d==1], edat=grid)$dens
  fphi0a=npudens(bws=bw7, tdat=phi0[d==0], edat=grid)$dens


  if (wgtfct==1){
    y1sq=(y^2)*r*d*((z1/Pz1^2)+(1-z1)/((1-Pz1)^2))
    y0sq=(y^2)*r*(1-d)*((z1/Pz1^2)+(1-z1)/((1-Pz1)^2))
    if (is.null(bw11)) {
      if (ruleofthumb!=1) bw11=npregbw(ydat=y1sq, xdat=phi1, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
      if (ruleofthumb==1) bw11=npudensbw(dat=phi1, ckertype="gaussian", bwmethod="normal-reference")$bw
    }
    if (is.null(bw12)) {
      if (ruleofthumb!=1) bw12=npregbw(ydat=y0sq, xdat=phi0, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
      if (ruleofthumb==1) bw12=npudensbw(dat=phi0, ckertype="gaussian", bwmethod="normal-reference")$bw
    }
    lambda1a=npreg(bws=bw11, tydat=y1sq, txdat=phi1, exdat=grid, regtype=rtype, ckertype="gaussian")$mean
    lambda0a=npreg(bws=bw12, tydat=y0sq, txdat=phi0, exdat=grid, regtype=rtype, ckertype="gaussian")$mean

    wgttemp=sqrt(lambda1a/fphi1a+lambda0a/fphi0a);wgttemp[wgttemp<trim]=trim
    wgta=grid/wgttemp
    ind=is.na(wgta)
    wgta=wgta[ind==0]
    cx=sum(wgta)
    latetemp=sum(late[ind==0]*(wgta/(cx)))

    temp=lm(y1sq~phi1); temp2=lm(y0sq~phi0);

    lambda0para=(cbind(1,grid)%*%coef(temp2))
    lambda1para=(cbind(1,grid)%*%coef(temp))
    wgttemppara=sqrt(lambda1para/fphi1para+lambda0para/fphi0para); wgttemppara[wgttemppara<trim]=trim
    wgtpara=grid/wgttemppara
    cxpar=sum(wgtpara)
    latepartemp=sum(latepar*(wgtpara/(cxpar)))
  }

  if (wgtfct!=1 & wgtfct!=3){
    wgttemp=sqrt(1/fphi1a+1/fphi0a);wgttemp[wgttemp<trim]=trim
    wgta=grid/wgttemp
    cx=sum(wgta)
    latetemp=sum(late*(wgta/(cx)))

    wgttemppara=sqrt(1/fphi1para+1/fphi0para); wgttemppara[wgttemppara<trim]=trim
    wgtpara=grid/wgttemppara
    cxpar=sum(wgtpara)
    latepartemp=sum(latepar*(wgtpara/(cxpar)))
  }

  if (wgtfct==3){
    latetemp=median(late[is.na(late)==0])
    latepartemp=median(latepar[is.na(latepar)==0])
  }

  itt=sum(y[r==1]*z1[r==1]/Pz1)/sum(z1[r==1]/Pz1)-sum(y[r==1]*(1-z1[r==1])/(1-Pz1))/sum((1-z1[r==1])/(1-Pz1));
  Pco=sum(d[r==1]*z1[r==1]/Pz1)/sum(z1[r==1]/Pz1)-sum(d[r==1]*(1-z1[r==1])/(1-Pz1))/sum((1-z1[r==1])/(1-Pz1)); n=length(r)
  latenaive=itt/Pco

  list(late=latetemp, latepar=latepartemp, latenaive=latenaive, latenaivepar=latenaive, phi1=phi1, phi1par=phi1par, phi0=phi0, phi0par=phi0par, bw1=bw1, bw2=bw2, bw3=bw3, bw4=bw4, bw5=bw5, bw6=bw6)
}



latenonrespxx=function(y,d,r,z1,z2, x=NULL, xpar=NULL, bres1=NULL, bres0=NULL, bwyz1=NULL, bwdz1=NULL, bwyz0=NULL, bwdz0=NULL, bwps=NULL, bwcox1=NULL, bwcox2=NULL, bw1=NULL, bw2=NULL, bw3=NULL, bw4=NULL, bw5=NULL, bw6=NULL, bw7=NULL, bw8=NULL, bw9=NULL, bw10=NULL, bw11=NULL, bw12=NULL, ruleofthumb=1, wgtfct=2, rtype="ll", numresprob=100, estlate=TRUE, trim=0.01){
  if (ncol(data.frame(x))>1 | is.null(ncol(x))==0) {xd1=x[d==1,]; xd0=x[d==0,]; xz11=x[z1==1,]; xz10=x[z1==0,]; xz11r1=x[z1==1 & r==1,]; xz10r1=x[z1==0 & r==1,]}
  if ( (ncol(data.frame(x))==1) & is.null(ncol(x))  ) {xd1=x[d==1]; xd0=x[d==0]; xz11=x[z1==1]; xz10=x[z1==0]; xz11r1=x[z1==1 & r==1]; xz10r1=x[z1==0 & r==1] }

  if (is.null(bres1)){
    if (ruleofthumb!=1) bres1=npregbw(ydat=factor(r[z1==1 ]), xdat=xz11, regtype="lc", ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bres1=npudensbw(dat=xz11, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  pres1=npreg(bws=bres1, tydat=factor(r[z1==1 ]), txdat=xz11, exdat=x, regtype="lc", ckertype="gaussian")$mean-1
  prespara1=pnorm(cbind(1,xpar)%*%coef(glm(formula=r[z1==1 ]~cbind(xpar)[z1==1,],family=binomial(probit))))
  if (is.null(bres0)){
    if (ruleofthumb!=1) bres0=npregbw(ydat=factor(r[z1==0 ]), xdat=xz11, regtype="lc", ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bres0=npudensbw(dat=xz10, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  pres0=npreg(bws=bres0, tydat=factor(r[z1==0 ]), txdat=xz10, exdat=x, regtype="lc", ckertype="gaussian")$mean-1
  prespara0=pnorm(cbind(1,xpar)%*%coef(glm(formula=r[z1==0 ]~cbind(xpar)[z1==0,],family=binomial(probit))))
  pres=z1*pres1+(1-z1)*pres0; prespara=z1*prespara1+(1-z1)*prespara0

  if (is.null(bwyz1)){
    if (ruleofthumb!=1) bwyz1<-npregbw(ydat=y[z1==1 & r==1], xdat=xz11r1, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bwyz1=npudensbw(dat=xz11r1, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  eyz1=npreg(bws=bwyz1, tydat=y[z1==1 & r==1], txdat=xz11r1, exdat=x, regtype=rtype, ckertype="gaussian")$mean
  eyz1par<-cbind(1,xpar)%*%coef(lm(y[z1==1 & r==1]~cbind(xpar)[z1==1 & r==1,]))

  if (is.null(bwyz0)){
    if (ruleofthumb!=1) bwyz0<-npregbw(ydat=y[z1==0 & r==1], xdat=xz10r1, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bwyz0=npudensbw(dat=xz10r1, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  eyz0=npreg(bws=bwyz0, tydat=y[z1==0 & r==1], txdat=xz10r1, exdat=x, regtype=rtype, ckertype="gaussian")$mean
  eyz0par<-cbind(1,xpar)%*%coef(lm(y[z1==0 & r==1]~cbind(xpar)[z1==0 & r==1,]))

  if (is.null(bwdz1)){
    if (ruleofthumb!=1) bwdz1<-npregbw(ydat=factor(d[z1==1 ]), xdat=xz11, regtype="lc", ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bwdz1=npudensbw(dat=xz11, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  edz1=npreg(bws=bwdz1, tydat=factor(d[z1==1 ]), txdat=xz11, exdat=x, regtype="lc", ckertype="gaussian")$mean-1
  edz1par<-pnorm(cbind(1,xpar)%*%coef(glm(d[z1==1 ]~cbind(xpar)[z1==1 ,], family = binomial(probit))))

  if (is.null(bwyz0)){
    if (ruleofthumb!=1) bwyz0<-npregbw(ydat=y[z1==0 ], xdat=xz10, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bwyz0=npudensbw(dat=xz10, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  edz0=npreg(bws=bwyz0, tydat=factor(d[z1==0 ]), txdat=xz10, exdat=x, regtype="lc", ckertype="gaussian")$mean-1
  edz0par<-pnorm(cbind(1,xpar)%*%coef(glm(d[z1==0 ]~cbind(xpar)[z1==0 ,], family = binomial(probit))))

  regd1=data.frame(z2[d==1],xd1); regd0=data.frame(z2[d==0],xd0); reg=data.frame(z2,x); n=length(d)
  if (is.null(bwps)) {
    if (ruleofthumb!=1) bwps<-npregbw(ydat=factor(z1), xdat=x, regtype="lc", ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bwps<-npudensbw(dat=x, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  Pz1=npreg(bws=bwps, tydat=factor(z1), txdat=x, regtype="lc", ckertype="gaussian")$mean-1
  Pz1para=fitted(glm(z1~xpar, family=binomial(probit)))
  dz1=d*z1
  if (is.null(bwcox1)) {
    if (ruleofthumb!=1) bwcox1<-npregbw(ydat=factor(dz1), xdat=x, regtype="lc", ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bwcox1<-npudensbw(dat=x, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  edz1=npreg(bws=bwcox1, tydat=factor(dz1), txdat=x, regtype="lc", ckertype="gaussian")$mean-1
  edz1para=fitted(glm(dz1~xpar, family=binomial(probit)))
  if (is.null(bwcox2)) {
    if (ruleofthumb!=1) bwcox2<-npregbw(ydat=d, xdat=x, regtype="lc", ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bwcox2<-npudensbw(dat=x, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  ed=npreg(bws=bwcox2, tydat=factor(d), txdat=x, regtype="lc", ckertype="gaussian")$mean-1
  edpara=fitted(glm(d~xpar, family=binomial(probit)))

  yphi=r*(z1-Pz1); yphi2=(z1-Pz1);
  yphipara=r*(z1-Pz1para); yphi2para=(z1-Pz1para);
  Pco=sum(d*z1/Pz1)/sum(z1/Pz1)-sum(d*(1-z1)/(1-Pz1))/sum((1-z1)/(1-Pz1));
  Pcopara=sum(d*z1/Pz1para)/sum(z1/Pz1para)-sum(d*(1-z1)/(1-Pz1para))/sum((1-z1)/(1-Pz1para));

  if (is.null(bw1)) {
    if (ruleofthumb!=1) bw1<-npregbw(ydat=yphi[d==1], xdat=regd1, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw1<-npudensbw(dat=regd1, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  phi1a<-(npreg(bws=bw1, tydat=yphi[d==1], txdat=regd1, exdat=reg, regtype=rtype, ckertype="gaussian")$mean)
  phi1apar<-cbind(1,z2,xpar)%*%coef(lm(yphipara[d==1]~z2[d==1]+cbind(xpar)[d==1,]))

  if (is.null(bw2)) {
    if (ruleofthumb!=1) bw2<-npregbw(ydat=yphi2[d==1], xdat=regd1, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw2<-bw1
  }
  phi1b<-(npreg(bws=bw2, tydat=yphi2[d==1], txdat=regd1, exdat=reg, regtype=rtype, ckertype="gaussian")$mean)
  phi1bpar<-cbind(1,z2,xpar)%*%coef(lm(yphi2para[d==1]~z2[d==1]+cbind(xpar)[d==1,]))
  phi1=phi1a/phi1b; phi1par=phi1apar/phi1bpar

  if (is.null(bw3)) {
    if (ruleofthumb!=1) bw3<-npregbw(ydat=yphi[d==0], xdat=regd0, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw3<-npudensbw(dat=regd0, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  phi0a<-(npreg(bws=bw3, tydat=yphi[d==0], txdat=regd0, exdat=reg, regtype=rtype, ckertype="gaussian")$mean)
  phi0apar<-cbind(1,z2,xpar)%*%coef(lm(yphipara[d==0]~z2[d==0]+cbind(xpar)[d==0,]))

  if (is.null(bw4)) {
    if (ruleofthumb!=1) bw4<-npregbw(ydat=yphi2[d==0], xdat=regd0, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw4<-bw3
  }
  phi0b<-(npreg(bws=bw4, tydat=yphi2[d==0], txdat=regd0, exdat=reg, regtype=rtype, ckertype="gaussian")$mean)
  phi0bpar<-cbind(1,z2,xpar)%*%coef(lm(yphi2para[d==0]~z2[d==0]+cbind(xpar)[d==0,]))

  phi0=phi0a/phi0b; phi0par=phi0apar/phi0bpar

  yrd=y*r*d*(z1-Pz1);yr1_d=y*r*(1-d)*(z1-Pz1)
  yrdpar=y*r*d*(z1-Pz1para);yr1_dpar=y*r*(1-d)*(z1-Pz1para)
  minselprob=max(min(phi1[is.na(phi1)==0]),min(phi0[is.na(phi0)==0]),0.01); maxselprob=min(max(phi1[is.na(phi1)==0]),max(phi0[is.na(phi0)==0]),1)
  regphi1=data.frame(phi1,x); regphi0=data.frame(phi0,x); varx=x

  if (is.null(bw5)) {
    if (ruleofthumb!=1) bw5<-npregbw(ydat=yrd, xdat=regphi1, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw5<-npudensbw(dat=regphi1, ckertype="gaussian", bwmethod="normal-reference")$bw
  }
  predval=seq(from=minselprob,to=maxselprob,length.out=numresprob)

  eyrd=c(); eyrdpar=c(); coef1=lm(yrdpar~phi1par+xpar)$coef
  for (i in 1:numresprob){
    eyrd<-cbind(eyrd,npreg(bws=bw5, tydat=yrd, txdat=regphi1, exdat=data.frame(rep(predval[i],n),varx), regtype=rtype, ckertype="gaussian")$mean)
    eyrdpar<-cbind(eyrdpar,cbind(1,predval[i],xpar)%*%coef1)
  }

  if (is.null(bw6)) {
    if (ruleofthumb!=1) bw6<-npregbw(ydat=yr1_d, xdat=regphi0, regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
    if (ruleofthumb==1) bw6<-npudensbw(dat=regphi0, ckertype="gaussian", bwmethod="normal-reference")$bw
  }

  eyr1_d=c(); eyr1_dpar=c(); coef2=lm(yr1_dpar~phi0par+xpar)$coef
  for (i in 1:numresprob){
    eyr1_d<-cbind(eyr1_d, npreg(bws=bw6, tydat=yr1_d, txdat=regphi0, exdat=data.frame(rep(predval[i],n),varx), regtype=rtype, ckertype="gaussian")$mean)
    eyr1_dpar<-cbind(eyr1_dpar, cbind(1,predval[i],xpar)%*%coef2)
  }

  if (estlate!=0) {
    denom=(Pz1*(1-Pz1)*Pco); denom[denom<trim]=trim;
    denompar=(Pz1para*(1-Pz1para)*Pcopara);  denompar[denompar<trim]=trim
  }
  if (estlate==0) {denom=( (edz1-ed*Pz1)); denom[denom<trim]=trim;
  denompar=( (edz1para-edpara*Pz1para)); denompar[denompar<trim]=trim
  }


  late=c(); latepar=c()
  for (i in 1:numresprob){
    if (estlate!=0){
      late=cbind(late, ((eyrd[,i]+eyr1_d[,i])/(denom*predval[i])))
      latepar=cbind(latepar,((eyrdpar[,i]+eyr1_dpar[,i])/(denompar*predval[i])))
    }
    if (estlate==0){
      # respective functions for the ATE
      late=cbind(late, ((eyrd[,i]+eyr1_d[,i])/( denom*predval[i])))
      latepar=cbind(latepar,((eyrdpar[,i]+eyr1_dpar[,i])/( denompar*predval[i])))
    }
  }

  if (is.null(bw7)) {
    if (ruleofthumb!=1) bw7<-npcdensbw(ydat=phi1[d==1], xdat=xd1, ckertype="gaussian", bwmethod="cv.ls")
    if (ruleofthumb==1) bw7<-npcdensbw(ydat=phi1[d==1], xdat=xd1, ckertype="gaussian", bwmethod="normal-reference")
  }


  if (is.null(bw8)) {
    if (ruleofthumb!=1) bw8<-npcdensbw(ydat=phi1par[d==1], xdat=xd1, ckertype="gaussian", bwmethod="cv.ls")
    if (ruleofthumb==1) bw8<-npcdensbw(ydat=phi1par[d==1], xdat=xd1, ckertype="gaussian", bwmethod="normal-reference")
  }

  if (is.null(bw9)) {
    if (ruleofthumb!=1) bw9<-npcdensbw(ydat=phi0[d==0], xdat=xd0, ckertype="gaussian", bwmethod="cv.ls")
    if (ruleofthumb==1) bw9<-npcdensbw(ydat=phi0[d==0], xdat=xd0, ckertype="gaussian", bwmethod="normal-reference")
  }

  if (is.null(bw10)) {
    if (ruleofthumb!=1) bw10<-npcdensbw(ydat=phi0par[d==0], xdat=xd0, ckertype="gaussian", bwmethod="cv.ls")
    if (ruleofthumb==1) bw10<-npcdensbw(ydat=phi0par[d==0], xdat=xd0, ckertype="gaussian", bwmethod="normal-reference")
  }

  grid=c(predval); lgrid=length(grid)
  fphi1para=c(); fphi0para=c(); fphi1a=c(); fphi0a=c();
  for (j in 1:lgrid){
    fphi1para=rbind(fphi1para,npcdens(bws=bw8, tydat=phi1par[d==1], txdat=xd1, eydat=rep(grid[j],n), exdat=x)$condens)
    fphi0para=rbind(fphi0para,npcdens(bws=bw10, tydat=phi0par[d==0], txdat=xd0, eydat=rep(grid[j],n), exdat=x)$condens)
    fphi1a=rbind(fphi1a,npcdens(bws=bw7, tydat=phi1[d==1], txdat=xd1, eydat=rep(grid[j],n), exdat=x)$condens)
    fphi0a=rbind(fphi0a,npcdens(bws=bw9, tydat=phi0[d==0], txdat=xd0, eydat=rep(grid[j],n), exdat=x)$condens)
  }

  if (wgtfct!=1 & wgtfct!=3){
    wtemp=sqrt(1/fphi1a+1/fphi0a); wtemp[wtemp<trim]=trim
    wgta=grid/wtemp
    cx=colSums(wgta)
    latetemp=mean(colSums(t(late)*(wgta/(rep(1,length(grid))%*%t(cx)))))

    wtemppara=sqrt(1/fphi1para+1/fphi0para)
    wtemppara[wtemppara<trim]=trim
    wgtpara=grid/wtemppara
    cxpar=colSums(wgtpara)
    latepartemp=mean(colSums(t(latepar)*(wgtpara/(rep(1,length(grid))%*%t(cxpar)))))
  }

  if (wgtfct==1){
    y1sq=(y^2)*r*d*((z1/Pz1^2)+(1-z1)/((1-Pz1)^2))
    y0sq=(y^2)*r*(1-d)*((z1/Pz1^2)+(1-z1)/((1-Pz1)^2))
    if (is.null(bw11)) {
      if (ruleofthumb!=1) bw11=npregbw(ydat=y1sq, xdat=data.frame(phi1,x), regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
      if (ruleofthumb==1) bw11=npudensbw(dat=data.frame(phi1,x), ckertype="gaussian", bwmethod="normal-reference")$bw
    }
    if (is.null(bw12)) {
      if (ruleofthumb!=1) bw12=npregbw(ydat=y0sq, xdat=data.frame(phi0,x), regtype=rtype, ckertype="gaussian", bwmethod="cv.ls")$bw
      if (ruleofthumb==1) bw12=npudensbw(dat=data.frame(phi0,x), ckertype="gaussian", bwmethod="normal-reference")$bw
    }
    lambda1a=c(); lambda0a=c();
    for (j in 1:lgrid){
      lambda1a=rbind(lambda1a,npreg(bws=bw11, tydat=y1sq, txdat=data.frame(phi1,x), exdat=data.frame(rep(grid[j],n),x), regtype=rtype, ckertype="gaussian")$mean)
      lambda0a=rbind(lambda0a,npreg(bws=bw12, tydat=y0sq, txdat=data.frame(phi0,x), exdat=data.frame(rep(grid[j],n),x), regtype=rtype, ckertype="gaussian")$mean)
    }
    wgttemp=sqrt(lambda1a/fphi1a+lambda0a/fphi0a); wgttemp[wgttemp<trim]=trim
    wgta=grid/wgttemp
    ind=rowSums(is.na(wgta))
    wgta=wgta[ind==0,]
    cx=colSums(wgta)
    latetemp=mean(colSums(t(late[,ind==0])*(wgta/(rep(1,length(grid[ind==0]))%*%t(cx)))))

    y1sqpar=(y^2)*r*d*((z1/Pz1para^2)+(1-z1)/((1-Pz1para)^2))
    y0sqpar=(y^2)*r*(1-d)*((z1/Pz1para^2)+(1-z1)/((1-Pz1para)^2))
    temp=lm(y1sqpar~cbind(phi1par,xpar)); temp2=lm(y0sqpar~cbind(phi0par,xpar));
    lambda0para=c(); lambda1para=c()
    for (j in 1:lgrid){
      lambda0para=rbind(lambda0para, t(cbind(1,rep(grid[j],n),xpar)%*%coef(temp2)))
      lambda1para=rbind(lambda1para, t(cbind(1,rep(grid[j],n),xpar)%*%coef(temp)))
    }
    wgttemppara=sqrt(lambda1para/fphi1para+lambda0para/fphi0para); wgttemppara[wgttemppara<trim]=trim
    wgtpara=grid/wgttemppara
    cxpar=colSums(wgtpara)
    latepartemp=mean(colSums(t(latepar)*(wgtpara/(rep(1,length(grid))%*%t(cxpar)))))
  }

  if (wgtfct==3){
    latetemp=c(); latepartemp=c()
    for (j in 1:n){
      latetemp=c(latetemp, median(late[j,is.na(late[j,])==0]))
      latepartemp=c(latepartemp, median(latepar[j,is.na(latepar[j,])==0]))
    }
    latetemp=mean(latetemp); latepartemp=mean(latepartemp)
  }

  catenaive=(eyz1-eyz0)/(edz1-edz0)
  catenaivepar=(eyz1par-eyz0par)/(edz1par-edz0par)

  if (estlate!=0){
    itt=sum(r*y *z1 /(pres*Pz1) )/sum(r*z1 /(pres*Pz1) )-sum(r*y *(1-z1 )/(pres*(1-Pz1) ))/sum(r*(1-z1 )/(pres*(1-Pz1) ));
    ittpara=sum(r*y *z1 / (prespara*Pz1para) )/sum(r*z1 / (prespara*Pz1para) )-sum(r*y *(1-z1 )/ (prespara*(1-Pz1para )))/sum(r*(1-z1 )/(prespara*(1-Pz1para) ));
    Pco=sum(d *z1 /Pz1 )/sum(z1 /Pz1 )-sum(d*(1-z1 )/(1-Pz1 ))/sum((1-z1 )/(1-Pz1 )); n=length(r)
    Pcopara=sum(d *z1 /Pz1para )/sum(z1 /Pz1para )-sum(d *(1-z1 )/(1-Pz1para ))/sum((1-z1 )/(1-Pz1para ));
    latenaive=itt/Pco; latenaivepar=ittpara/Pcopara;
  }
  if (estlate==0){
    latenaive=mean(catenaive*r/pres)
    latenaivepar=mean(catenaivepar*r/prespara)
    latenaive=mean(catenaive)
    latenaivepar=mean(catenaivepar)
  }
  list(late=latetemp, latepar=latepartemp, latenaive=latenaive, latenaivepar=latenaivepar, phi1=phi1, phi1par=phi1par, phi0=phi0, phi0par=phi0par, bres1=bres1, bres0=bres0,  bwyz1=bwyz1, bwdz1=bwdz1, bwyz0=bwyz0, bwdz0=bwdz0, bwps=bwps, bwcox1=bwcox1, bwcox2=bwcox2, bw1=bw1, bw2=bw2, bw3=bw3, bw4=bw4, bw5=bw5, bw6=bw6)
}


latenonrespxxfct<-function(y,d,r,z1,z2, x=NULL, xpar=NULL, bres1=NULL, bres0=NULL,  bwyz1=NULL, bwdz1=NULL, bwyz0=NULL, bwdz0=NULL, bwps=NULL, bwcox1=NULL, bwcox2=NULL, bw1=NULL, bw2=NULL, bw3=NULL, bw4=NULL, bw5=NULL, bw6=NULL, bw7=NULL, bw8=NULL, bw9=NULL, bw10=NULL, bw11=NULL, bw12=NULL, ruleofthumb=1, wgtfct=2, rtype="ll", numresprob=100, estlate=TRUE, trim=0.01){
  if ((is.null(x)==0) & (is.null(xpar)==0)) out=latenonrespxx(y=y,d=d,r=r,z1=z1,z2=z2, x=x, xpar=xpar, bres1=bres1, bres0=bres0,  bwyz1=bwyz1, bwdz1=bwdz1, bwyz0=bwyz0, bwdz0=bwdz0, bwps=bwps, bwcox1=bwcox1, bwcox2=bwcox2, bw1=bw1, bw2=bw2, bw3=bw3, bw4=bw4, bw5=bw5, bw6=bw6,  bw7=bw7, bw8=bw8, bw9=bw9, bw10=bw10, bw11=bw11, bw12=bw12, ruleofthumb=ruleofthumb, wgtfct=wgtfct, rtype=rtype, numresprob=numresprob, estlate=estlate, trim=trim)
  if (is.null(x) | is.null(xpar)) out=latenonresp(y=y,d=d,r=r,z1=z1,z2=z2, bw1=bw1, bw2=bw2, bw3=bw3, bw4=bw4, bw5=bw5, bw6=bw6,   bw7=bw7, bw8=bw8, bw9=bw9, bw10=bw10, bw11=bw11, bw12=bw12, ruleofthumb=ruleofthumb, wgtfct=wgtfct, rtype=rtype, numresprob=numresprob, trim=trim)
  results=c(out$late, out$latepar, out$latenaive, out$latenaivepar)
  list(results=results, bres1=out$bres1, bres0=out$bres0,  bwyz1=out$bwyz1, bwdz1=out$bwdz1, bwyz0=out$bwyz0, bwdz0=out$bwdz0, bwps=out$bwps, bwcox1=out$bwcox1, bwcox2=out$bwcox2, bw1=out$bw1, bw2=out$bw2, bw3=out$bw3, bw4=out$bw4, bw5=out$bw5, bw6=out$bw6, phi1=out$phi1, phi1par=out$phi1par, phi0=out$phi0, phi0par=out$phi0par)
}

bootstrap.late.nr<-function(y,d,r,z1,z2, x=NULL, xpar=NULL, bres1=NULL, bres0=NULL,  bwyz1=NULL, bwdz1=NULL, bwyz0=NULL, bwdz0=NULL, bwps=NULL, bwcox1=NULL, bwcox2=NULL, bw1=NULL, bw2=NULL, bw3=NULL, bw4=NULL, bw5=NULL, bw6=NULL, bw7=NULL, bw8=NULL, bw9=NULL, bw10=NULL, bw11=NULL, bw12=NULL, ruleofthumb=1, wgtfct=2, rtype="ll", numresprob=100, boot=1999, estlate=TRUE, trim=0.01){
  mc=c(); temp=c(); j=1
  while(j<=boot){
    sboot<-sample(1:length(d),length(d),TRUE)
    z1b<-z1[sboot]; z2b<-z2[sboot]; db<-d[sboot]; yb=y[sboot]; rb=r[sboot]
    if ((is.null(x)==0) & (is.null(xpar)==0)) {
      if (length(x)==length(d)) xb<-x[sboot]
      if (length(x)!=length(d)) xb<-x[sboot,]
      if (length(xpar)==length(d)) xparb<-xpar[sboot]
      if (length(xpar)!=length(d)) xparb<-xpar[sboot,]
      temp<-latenonrespxxfct(y=yb,d=db,r=rb,z1=z1b,z2=z2b, x=xb, xpar=xparb, bres1=bres1, bres0=bres0,  bwyz1=bwyz1, bwdz1=bwdz1, bwyz0=bwyz0, bwdz0=bwdz0, bwps=bwps, bwcox1=bwcox1, bwcox2=bwcox2, bw1=bw1, bw2=bw2, bw3=bw3, bw4=bw4, bw5=bw5, bw6=bw6,  bw7=bw7, bw8=bw8, bw9=bw9, bw10=bw10, bw11=bw11, bw12=bw12, ruleofthumb=ruleofthumb, wgtfct=wgtfct, rtype=rtype, numresprob=numresprob, estlate=estlate, trim=trim)$results
    }
    if (is.null(x) | is.null(x)) {
      temp<-latenonrespxxfct(y=yb,d=db,r=rb,z1=z1b,z2=z2b, x=NULL, xpar=NULL, bwps=NULL, bwcox1=NULL, bwcox2=NULL, bw1=bw1, bw2=bw2, bw3=bw3, bw4=bw4, bw5=bw5, bw6=bw6, bw7=bw7, bw8=bw8, bw9=bw9, bw10=bw10, bw11=bw11, bw12=bw12, ruleofthumb=ruleofthumb, wgtfct=wgtfct, rtype=rtype, numresprob=numresprob, trim=trim)$results
    }
    if (is.na(sum(temp))==0){
      mc<-rbind(mc,temp)
      j=j+1
    }
  }
  mc
}

# function for RDD kernel regression with covariates
rdd.x.est=function(y,z,x, bw0, bw1, regtype, bwz){
  d=1*(z>=0)
  xz=data.frame(x,z)
  xzcutoff=data.frame(x,rep(0,length(d)))
  xz0=xz[d==0,]; xz1=xz[d==1,]; d1=d[d==1]; d0=d[d==0]; y1=y[d==1]; y0=y[d==0];
  reg0<-npreg(bws=bw0, tydat=y0, txdat=xz0, exdat=xzcutoff, ckertype="epanechnikov", regtype=regtype)$mean
  reg1<-npreg(bws=bw1, tydat=y1, txdat=xz1, exdat=xzcutoff, ckertype="epanechnikov", regtype=regtype)$mean
  kernwgt=npksum(bws=bwz, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
  mu2=0.1; mu1=3/16
  kernwgt=(mu2-mu1*z)*kernwgt;
  effect=(sum((reg1-reg0)*kernwgt))/(sum(kernwgt))
  effect
}

# bootstrap function for RDD kernel regression with covariates
rdd.x.boot<-function(y,z,x, bw0, bw1, bwz, boot=1999, regtype){
  obs<-length(y)
  mc=c(); i=1
  while(i<=boot){
    sboot<-sample(1:obs,obs,TRUE)
    yb=y[sboot]
    zb<-z[sboot]
    if (length(x)==length(y)) xb<-x[sboot]
    if (length(x)!=length(y)) xb<-x[sboot,]

    est<-c(rdd.x.est(y=yb,z=zb,x=xb, bw0=bw0, bw1=bw1, regtype=regtype, bwz=bwz))
    if (sum(is.na(est))==0) mc<-c(mc, est)
    i=i+1
  }
  mc
}
