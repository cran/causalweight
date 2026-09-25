hdtreat_newscore = function(y, d, x, MLmethod = "lasso", k = 3, zeta, seed){
  ybin <- 1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  if (length(d) < k*3) stop("estimation needs at least 3*k observations for k-fold cross-fitting")
  x <- data.frame(x)
  set.seed(seed)
  idx <- sample(length(d), replace=FALSE)
  folds <- split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))
  score <- c()
  # cross-fitting procedure that splits sample in training and testing data
  for (i in 1:k){
    tesample <- folds[[i]]
    trsample=idx[!(idx %in% tesample)]
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
  # se <- sqrt((mean(((tscores-cscores)^2-effect)^2)+var(scorestreat[,4]))/length(tscores))
  se <- sqrt((mean(((tscores-cscores)^2+scorestreat[,4]-effect)^2))/length(tscores))
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
  set.seed(seed)
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
  randomf <- ranger(scorestr~., data = data1, num.trees = 200, importance = "impurity_corrected")
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

#' @exportS3Method NULL
#' @noRd
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

hdmed=function(y,d,m,x,k=3, trim=0.05, normalized=TRUE, MLmethod="lasso"){
  if (length(d)<k*3) stop("estimator needs at least 3*k observations")
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)
  set.seed(1); idx = sample(length(d), replace=FALSE)
  folds = split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))
  score=c(); selall=c()
  # crossfitting procedure that splits sample in training an testing data
  for (i in 1:k){
    tesample=folds[[i]]
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
    
    if (MLmethod=="lasso"){
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
    }
    if (MLmethod!="lasso") {
      xtr=data.frame(xtr); xte=data.frame(xte); xtr1=data.frame(xtr1); xtr0=data.frame(xtr0); xtr11=data.frame(xtr11); xtr10=data.frame(xtr10); xtr01=data.frame(xtr01); xtr00=data.frame(xtr00)
      # predict Pr(M=1|D=1,X) in test data
      pm1=MLfunct(y=mtr1, x=xtr1, MLmethod=MLmethod,  ybin=1)
      pm1te=predict(pm1, xte, onlySL = TRUE)$pred
      # predict Pr(M=1|D=0,X) in test data
      pm0=MLfunct(y=mtr0, x=xtr0, MLmethod=MLmethod,  ybin=1)
      pm0te=predict(pm0, xte, onlySL = TRUE)$pred
      # predict Pr(D=1|X) in test data
      pd=MLfunct(y=dtr, x=xtr, MLmethod=MLmethod,  ybin=1)
      pdte=predict(pd, xte, onlySL = TRUE)$pred
      # predict E(Y| D=1, M=1, X) in test data
      eymx11=MLfunct(y=ytr11, x=xtr11, MLmethod=MLmethod,  ybin=ybin)
      eymx11te=predict(eymx11, xte, onlySL = TRUE)$pred
      # predict E(Y| D=0, M=1, X) in test data
      eymx01=MLfunct(y=ytr01, x=xtr01, MLmethod=MLmethod,  ybin=ybin)
      eymx01te=predict(eymx01, xte, onlySL = TRUE)$pred
      # predict E(Y| D=1, M=0, X) in test data
      eymx10=MLfunct(y=ytr10, x=xtr10, MLmethod=MLmethod,  ybin=ybin)
      eymx10te=predict(eymx10, xte, onlySL = TRUE)$pred
      # predict E(Y| D=0, M=0, X) in test data
      eymx00=MLfunct(y=ytr00, x=xtr00, MLmethod=MLmethod,  ybin=ybin)
      eymx00te=predict(eymx00, xte, onlySL = TRUE)$pred
      #  predict E(Y|D=1, X) in test data
      eyx1=MLfunct(y=ytr1, x=xtr1, MLmethod=MLmethod,  ybin=ybin)
      eyx1te=predict(eyx1, xte, onlySL = TRUE)$pred
      #  predict E(Y|D=0, X) in test data
      eyx0=MLfunct(y=ytr0, x=xtr0, MLmethod=MLmethod,  ybin=ybin)
      eyx0te=predict(eyx0, xte, onlySL = TRUE)$pred
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
hdmedalt=function(y,d,m,x, trim=0.05, fewsplits=FALSE, normalized=TRUE, MLmethod="lasso", debiasfirst=TRUE){
  if (length(d)<9) stop("estimator needs at least 9 observations")   
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)
  xm=cbind(x,m)
  set.seed(1); idx = sample(length(d))
  folds3 = split(idx, cut(seq_along(idx), breaks = 3, labels = FALSE))
  sample1 = folds3[[1]]; sample2 = folds3[[2]]; sample3 = folds3[[3]]
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
    if (MLmethod!="lasso") {
      xm=data.frame(xm); x=data.frame(x)
      pmx=MLfunct(y=d[trsample], x=xm[trsample,], MLmethod=MLmethod,  ybin=1)
      pmxte=predict(pmx, xm[tesample,], onlySL = TRUE)$pred
      px=MLfunct(y=d[trsample], x=x[trsample,], MLmethod=MLmethod,  ybin=1)
      pxte=predict(px, x[tesample,], onlySL = TRUE)$pred
      eymx1=MLfunct(y=y[musample[d[musample]==1]], x=xm[musample[d[musample]==1],], MLmethod=MLmethod,  ybin=ybin)
      # predict E(Y|M,X,D=1) in test data
      eymx1te=predict(eymx1, xm[tesample,], onlySL = TRUE)$pred
      # predict E(Y|M,X,D=1) in delta sample
      eymx1trte=predict(eymx1, xm[deltasample,], onlySL = TRUE)$pred
      dtrte=d[deltasample]; xtrte=data.frame(x)[deltasample,]
      regweymx1=MLfunct(y=eymx1trte[dtrte==0], x=xtrte[dtrte==0,], MLmethod=MLmethod, ybin=0)
      # predict E[E(Y|M,X,D=1)|D=0,X] in the test data
      regweymx1te=predict(regweymx1, x[tesample,], onlySL = TRUE)$pred
      
      if (debiasfirst==TRUE){
        # q(S)=f0(S)/f1(S) via Bayes' rule, with pi(X) and p(S) refit on musample
        # only and predicted out-of-sample onto deltasample
        pxa=MLfunct(y=d[musample], x=x[musample,], MLmethod=MLmethod, ybin=1)
        pxadelta=predict(pxa, x[deltasample,], onlySL = TRUE)$pred
        pmxa=MLfunct(y=d[musample], x=xm[musample,], MLmethod=MLmethod, ybin=1)
        pmxadelta=predict(pmxa, xm[deltasample,], onlySL = TRUE)$pred
        # floor/cap propensities so qhat cannot explode
        pxadelta=pmin(pmax(pxadelta,trim),1-trim); pmxadelta=pmin(pmax(pmxadelta,trim),1-trim)
        qhatdelta=(pxadelta/(1-pxadelta))*((1-pmxadelta)/pmxadelta)
        # doubly robust pseudo-outcome for tau, D=1 units of delta sample
        regweymx1d1=predict(regweymx1, x[deltasample[dtrte==1],], onlySL = TRUE)$pred
        ydr1=regweymx1d1 + qhatdelta[dtrte==1]*(y[deltasample[dtrte==1]]-eymx1trte[dtrte==1])
        tau1=MLfunct(y=ydr1, x=xtrte[dtrte==1,], MLmethod=MLmethod, ybin=0)
        tau1te=predict(tau1, x[tesample,], onlySL = TRUE)$pred
      }
      # debiasing off: tau1te collapses to the plain nested regression
      if (debiasfirst!=TRUE) tau1te=regweymx1te
      
      #  fit E(Y|X,D=1) in total of training data with D=1 by running Y~X
      temp=MLfunct(y=y[trsample[d[trsample]==1]], x=x[trsample[d[trsample]==1],], MLmethod=MLmethod, ybin=ybin)
      # predict E(Y|X,D=1) in the test data
      eyx1te=predict(temp, x[tesample,], onlySL = TRUE)$pred
      # fit E(Y|M,X,D=0) in first training data
      eymx0=MLfunct(y=y[musample[d[musample]==0]], x=xm[musample[d[musample]==0],], MLmethod=MLmethod, ybin=ybin)
      # predict E(Y|M,X,D=0) in test data
      eymx0te=predict(eymx0, xm[tesample,], onlySL = TRUE)$pred
      # predict E(Y|M,X,D=0) in delta sample
      eymx0trte=predict(eymx0, xm[deltasample,], onlySL = TRUE)$pred
      regweymx0=MLfunct(y=eymx0trte[dtrte==1], x=xtrte[dtrte==1,], MLmethod=MLmethod, ybin=0)
      regweymx0te=predict(regweymx0, x[tesample,], onlySL = TRUE)$pred
      
      # mirror of the tau1 debiasing above, using 1/qhat for the D=0 side
      if (debiasfirst==TRUE){
        regweymx0d0=predict(regweymx0, x[deltasample[dtrte==0],], onlySL = TRUE)$pred
        ydr0=regweymx0d0 + (1/qhatdelta[dtrte==0])*(y[deltasample[dtrte==0]]-eymx0trte[dtrte==0])
        tau0=MLfunct(y=ydr0, x=xtrte[dtrte==0,], MLmethod=MLmethod, ybin=0)
        tau0te=predict(tau0, x[tesample,], onlySL = TRUE)$pred
      }
      if (debiasfirst!=TRUE) tau0te=regweymx0te
      
      temp=MLfunct(y=y[trsample[d[trsample]==0]], x=x[trsample[d[trsample]==0],], MLmethod=MLmethod, ybin=ybin)
      # predict E(Y|X,D=0) in the test data
      eyx0te=predict(temp, x[tesample,], onlySL = TRUE)$pred
    }
    if (MLmethod=="lasso") {
      # fit Pr(D=1|M,X) in total of training data and predict Pr(D=1|M,X) in test data
      pmx=rlassologit(d[trsample]~xm[trsample,])
      pmxte=predict(pmx, xm[tesample,], type="response")
      # fit Pr(D=1|X) in total of training data and predict Pr(D=1|X) in test data
      px=rlassologit(d[trsample]~x[trsample,])
      pxte=predict(px, x[tesample,], type="response")
      if (ybin!=1){
        eymx1=rlasso(y[musample[d[musample]==1]]~xm[musample[d[musample]==1],])
        eymx1te=predict(eymx1, xm[tesample,])
        eymx1trte=predict(eymx1, xm[deltasample,])
      }
      if (ybin==1){
        eymx1=rlassologit(y[musample[d[musample]==1]]~xm[musample[d[musample]==1],])
        eymx1te=predict(eymx1, xm[tesample,], type="response")
        eymx1trte=predict(eymx1, xm[deltasample,], type="response")
      }
      dtrte=d[deltasample]; xtrte=x[deltasample,]
      regweymx1=rlasso(eymx1trte[dtrte==0]~xtrte[dtrte==0,])
      regweymx1te=predict(regweymx1, x[tesample,])
      
      # q(S)=f0(S)/f1(S) via Bayes' rule; pi(X), p(S) refit on musample only and
      # predicted out-of-sample onto deltasample
      if (debiasfirst==TRUE){
        pxa=rlassologit(d[musample]~x[musample,])
        pxadelta=predict(pxa, x[deltasample,], type="response")
        pmxa=rlassologit(d[musample]~xm[musample,])
        pmxadelta=predict(pmxa, xm[deltasample,], type="response")
        pxadelta=pmin(pmax(pxadelta,trim),1-trim); pmxadelta=pmin(pmax(pmxadelta,trim),1-trim)
        qhatdelta=(pxadelta/(1-pxadelta))*((1-pmxadelta)/pmxadelta)
        regweymx1d1=predict(regweymx1, x[deltasample[dtrte==1],])
        ydr1=regweymx1d1 + qhatdelta[dtrte==1]*(y[deltasample[dtrte==1]]-eymx1trte[dtrte==1])
        tau1=rlasso(ydr1~xtrte[dtrte==1,])
        tau1te=predict(tau1, x[tesample,])
      }
      if (debiasfirst!=TRUE) tau1te=regweymx1te
      
      if (ybin!=1){
        #  fit E(Y|X,D=1) in total of training data with D=1 by running Y~X
        temp=rlasso(y[trsample[d[trsample]==1]]~x[trsample[d[trsample]==1],])
        eyx1te=predict(temp, x[tesample,])
        eymx0=rlasso(y[musample[d[musample]==0]]~xm[musample[d[musample]==0],])
        eymx0te=predict(eymx0, xm[tesample,])
        eymx0trte=predict(eymx0, xm[deltasample,])
      }
      if (ybin==1){
        temp=rlassologit(y[trsample[d[trsample]==1]]~x[trsample[d[trsample]==1],])
        eyx1te=predict(temp, x[tesample,], type="response")
        eymx0=rlassologit(y[musample[d[musample]==0]]~xm[musample[d[musample]==0],])
        eymx0te=predict(eymx0, xm[tesample,], type="response")
        eymx0trte=predict(eymx0, xm[deltasample,], type="response")
      }
      
      regweymx0=rlasso(eymx0trte[dtrte==1]~xtrte[dtrte==1,])
      regweymx0te=predict(regweymx0, x[tesample,])
      
      # mirror for the D=0 side, using 1/qhat
      if (debiasfirst==TRUE){
        regweymx0d0=predict(regweymx0, x[deltasample[dtrte==0],])
        ydr0=regweymx0d0 + (1/qhatdelta[dtrte==0])*(y[deltasample[dtrte==0]]-eymx0trte[dtrte==0])
        tau0=rlasso(ydr0~xtrte[dtrte==0,])
        tau0te=predict(tau0, x[tesample,])
      }
      if (debiasfirst!=TRUE) tau0te=regweymx0te
      
      if (ybin!=1){
        #  fit E(Y|X,D=0) in total of training data with D=0 by running Y~X
        temp=rlasso(y[trsample[d[trsample]==0]]~x[trsample[d[trsample]==0],])
        eyx0te=predict(temp, x[tesample,])
      }
      if (ybin==1){
        temp=rlassologit(y[trsample[d[trsample]==0]]~x[trsample[d[trsample]==0],])
        eyx0te=predict(temp, x[tesample,], type="response")
      }
    }
    # select observations satisfying trimming restriction
    sel= 1*((((1-pmxte)*pxte)>=trim) & ((1-pxte)>=trim)  & (pxte>=trim) &  (((pmxte*(1-pxte)))>=trim)   )
    
    # select elements of the score functions (tau1te, tau0te appended as columns 11-12)
    score=rbind(score, cbind(dte, pmxte, pxte, yte, eymx0te, regweymx0te, eyx0te, eymx1te, regweymx1te, eyx1te, tau1te, tau0te)[sel==1,])
    
    # collect selection dummies
    selall=c(selall,sel)
  }
  # compute scores for potential outcomes
  # y1m0/y0m1 use the debiased tau (score[,11]/score[,12]) plus the extra QR
  # correction term D*(tau_n-tau)/pi(X); with debiasfirst=FALSE these collapse
  # exactly to the original formulas since tau1te=regweymx1te, tau0te=regweymx0te
  if (normalized==FALSE) {
    y0m1=((1-score[,1])*score[,2]/((1-score[,2])*score[,3])*(score[,4]-score[,5])+score[,1]/score[,3]*(score[,5]-score[,6])+(1-score[,1])*(score[,6]-score[,12])/(1-score[,3])+score[,12])
    y0m0=(score[,7] + (1-score[,1])*(score[,4]-score[,7])/(1-score[,3]))
    y1m0=(score[,1]*(1-score[,2])/(score[,2]*(1-score[,3]))*(score[,4]-score[,8])+(1-score[,1])/(1-score[,3])*(score[,8]-score[,9])+score[,1]*(score[,9]-score[,11])/score[,3]+score[,11])
    y1m1=(score[,10] + score[,1]*(score[,4]-score[,10])/score[,3])
  }
  if (normalized!=FALSE) {
    nobs=nrow(score)
    sumscore1=sum((1-score[,1])*score[,2]/((1-score[,2])*score[,3]))
    sumscore2=sum(score[,1]/score[,3])
    sumscore3=sum((1-score[,1])/(1-score[,3]))
    sumscore4=sum(score[,1]*(1-score[,2])/(score[,2]*(1-score[,3])))
    y0m1=(nobs*(1-score[,1])*score[,2]/((1-score[,2])*score[,3])*(score[,4]-score[,5]))/sumscore1+(nobs*score[,1]/score[,3]*(score[,5]-score[,6]))/sumscore2+(nobs*(1-score[,1])*(score[,6]-score[,12])/(1-score[,3]))/sumscore3+score[,12]
    y0m0=score[,7] + (nobs*(1-score[,1])*(score[,4]-score[,7])/(1-score[,3]))/sumscore3
    y1m0=(nobs*score[,1]*(1-score[,2])/(score[,2]*(1-score[,3]))*(score[,4]-score[,8]))/sumscore4+(nobs*(1-score[,1])/(1-score[,3])*(score[,8]-score[,9]))/sumscore3+(nobs*score[,1]*(score[,9]-score[,11])/score[,3])/sumscore2+score[,11]
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

hddyntreat=function(y2,d1,d2,x0,x1, s=NULL, trim=0.05, MLmethod="lasso", fewsplits=fewsplits, debiasfirst=debiasfirst){
  if (length(y2)<9) stop("estimator needs at least 9 observations")   
  ybin=1*(length(unique(y2))==2 & min(y2)==0 & max(y2)==1)  # check if binary outcome
  x0=data.frame(x0); x0x1=data.frame(x0,x1); d1x0x1=data.frame(d1,x0x1);
  # crossfitting procedure that splits sample in training an testing data
  set.seed(1); idx = sample(length(y2))
  folds3 = split(idx, cut(seq_along(idx), breaks = 3, labels = FALSE))
  sample1 = folds3[[1]]; sample2 = folds3[[2]]; sample3 = folds3[[3]]
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
    
    if (debiasfirst==TRUE){
      p2a  = MLfunct(y=d2[trsample1], x=x0x1[trsample1,], d1=d1[trsample1], MLmethod=MLmethod, ybin=1)  #period-2 propensity, fit on trsample1 only, predicted out-of-sample onto trsample2
      p2tr2 = predict(p2a, x0x1[trsample2,], onlySL = TRUE)$pred
      p2tr2=pmax(p2tr2,trim)
      y2d1d2tr2dr = y2d1d2tr2 + d2[trsample2]*(y2[trsample2]-y2d1d2tr2)/p2tr2 # debiased pseudo-outcome 
      y1d1   = MLfunct(y=y2d1d2tr2dr, x=x0[trsample2,], d1=d1[trsample2], MLmethod=MLmethod, ybin=0)
    }
    if (debiasfirst!=TRUE) y1d1=MLfunct(y=y2d1d2tr2, x=x0[trsample2,], d1=d1[trsample2], MLmethod=MLmethod, ybin=0)
    
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
  if (length(d)<k*3) stop("estimator needs at least 3*k observations")  
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  x=data.frame(x)
  set.seed(1); idx = sample(length(d), replace=FALSE)
  folds = split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))
  score=c();
  # crossfitting procedure that splits sample in training an testing data
  for (i in 1:k){
    tesample=folds[[i]]
    trsample=idx[!(idx %in% tesample)]
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
  if (k < 3) {
    warning("k must be at least 3. Resetting k = 3.")
    k <- 3
  }
  ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  x=data.frame(x)
  dx=data.frame(d,x)
  score=c();
  # crossfitting procedure that splits sample in training an testing data
  nobs <- length(y); set.seed(1); idx <- sample(nobs)
  folds <- split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))
  for (i in 1:k) {
    # Rotate fold indices
    fold_roles <- c(i:k, if (i > 1) 1:(i - 1) else integer(0))
    # Assign roles
    tesample <- folds[[fold_roles[1]]]
    k_remaining <- k - 1
    n1 <- max(1, floor(k_remaining / 2.1))
    n2 <- k_remaining - n1
    trsample1_indices <- fold_roles[2:(1 + n1)]
    trsample2_indices <- fold_roles[(2 + n1):k]
    trsample1 <- unlist(folds[trsample1_indices], use.names = FALSE)
    trsample2 <- unlist(folds[trsample2_indices], use.names = FALSE)
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

MLmean = function(y, x, d, MLmethod = "lasso", k = 3, zeta, seed){
  ybin <- 1*(length(unique(y))==2 & min(y)==0 & max(y)==1)  # check if binary outcome
  if (length(d) < k*3) stop("estimation needs at least 3*k observations for k-fold cross-fitting")
  x <- data.frame(x)
  set.seed(seed)
  idx <- sample(length(d), replace=FALSE)
  folds <- split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))
  score <- c()
  # cross-fitting procedure that splits sample in training and testing data
  for (i in 1:k){
    tesample <- folds[[i]]
    trsample=idx[!(idx %in% tesample)]
    eydx <- MLfunct(y = y[trsample], x = x[trsample,], MLmethod = MLmethod, ybin = ybin)
    eydxte <- predict(eydx, x[tesample,], onlySL = TRUE)$pred  #predict conditional outcome in test data
    score <- rbind(score, cbind(eydxte,zeta[tesample]))
  }
  score <- score[order(idx),]
  score
}


hdtest = function(y1, y0, d, x, trim = 0.01, MLmethod = "lasso", k = 3) {
  
  # Check if the outcome is binary
  ybin = 1 * (length(unique(y1)) == 2 & min(y1) == 0 & max(y1) == 1 &
                length(unique(y0)) == 2 & min(y0) == 0 & max(y0) == 1)
  if (length(d) < k*3) stop("estimator needs at least 3*k observations for k-fold cross-fitting") 
  ## Convert vectors to data frames for processing
  # Dataframe with covariates for did nuisance parameters
  x = data.frame(x)
  
  # Dataframe with covariates and outcome in period 0 for selobs nuisance parameters
  xy0 = data.frame(x, y0)
  
  # Counterfactual treatment status
  d0 = 1 - d
  
  # Difference in outcome period 0 to period 1
  y10 = y1 - y0
  
  # Split sample into k folds
  set.seed(1)
  idx = sample(length(d), replace = FALSE)
  folds <- split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))
  
  # Initialize an empty vector to store nuisance parameters
  nuisance = c()
  
  # Cross-fitting procedure to split sample into training and testing data
  for (i in 1:k) {
    tesample <- folds[[i]]
    trsample=idx[!(idx %in% tesample)]
    
    # Estimate OOB-propensity scores
    # Selobs: train propensity score learner
    p = MLfunct(y = d[trsample], x = xy0[trsample, ], MLmethod = MLmethod, ybin = 1)
    
    # Selobs: predict OOB propensity scores
    pte = predict(p, xy0[tesample, ], onlySL = TRUE)$pred
    
    # Did: train propensity score learner
    pi = MLfunct(y = d[trsample], x = x[trsample, ], MLmethod = MLmethod, ybin = 1)
    
    # Did: predict OOB propensity scores
    pite = predict(pi, x[tesample, ], onlySL = TRUE)$pred
    
    # Estimate OOB-conditional outcome
    # Selobs: train conditional outcome learner
    mu = MLfunct(y = y1[trsample], x = xy0[trsample, ], d1 = d0[trsample], MLmethod = MLmethod, ybin = ybin)
    
    # Selobs: predict OOB conditional outcome learner
    mute = predict(mu, xy0[tesample, ], onlySL = TRUE)$pred
    
    # Did: train conditional outcome learner
    m = MLfunct(y = y10[trsample], x = x[trsample, ], d1 = d0[trsample], MLmethod = MLmethod, ybin = ybin)
    
    # Selobs: predict OOB conditional outcome learner
    mte = predict(m, x[tesample, ], onlySL = TRUE)$pred
    
    # Find observations not satisfying trimming restriction
    trimmed = 1 * ((pte > 1 - trim) | (pite > 1 - trim))
    nuisance = rbind(nuisance, cbind(d[tesample], y1[tesample], y0[tesample], pte, pite, mute, mte, trimmed))
  }
  
  nuisance = nuisance[order(idx), ]
  return(nuisance)
}


## ===========================================================================
## Helper functions for medqteDML (natural direct/indirect QTE via DML),
## implementing Hsu, Huber, and Yen (2026): "Estimation of Direct and
## Indirect Quantile Treatment Effects with Double Machine Learning".
## ===========================================================================

## -----------------------------------------------------------------------
## utilities (renamed from function_med_qte.R for package naming style;
## behavior unchanged)
## -----------------------------------------------------------------------

## finds the tau-th quantile of a random variable from a grid of (a, P(Y<=a))
## pairs, using the rearrangement of Chernozhukov, Fernandez-Val, Galichon
## (2010) to enforce monotonicity of the estimated cdf before inversion
findqx=function(px, ax, taux){
  stopifnot(length(px) == length(ax))
  stopifnot(!is.unsorted(ax))
  stopifnot(!is.unsorted(taux))
  stopifnot(all(px >= 0 & px <= 1), all(taux > 0 & taux < 1))
  px=sort(px)                                      # rearrangement: sort the cdf values
  datax=data.frame(p = px, a = ax)
  result=aggregate(a ~ p, FUN = min, data = datax)  # smallest a for each p
  result1=approx(result$p, result$a, xout = taux, rule = 2)
  return(sort(result1$y))
}

## generates k-fold sample index blocks (n must exceed kfold)
kfoldind=function(n, kfold){
  stopifnot(n > kfold)
  n_pred=floor(n/kfold)
  stopifnot(n_pred >= 1)
  ind_st=seq(1, n, by = n_pred)
  ind_end=seq(n_pred, n, by = n_pred)
  if(length(ind_end) < length(ind_st)) ind_end=c(ind_end, n)
  ind_pred=cbind(ind_st, ind_end)
  while(nrow(ind_pred) > kfold){
    ind_pred=ind_pred[-nrow(ind_pred), , drop = FALSE]
    ind_pred[nrow(ind_pred), 2]=n
  }
  stopifnot(nrow(ind_pred) == kfold)
  fold_sizes=ind_pred[,2] - ind_pred[,1] + 1
  if(any(fold_sizes < 2)) warning("some folds have fewer than 2 observations - check n vs kfold")
  data.frame(ind_pred)
}

## multiplier bootstrap draw of a cdf across a grid, given the point estimate
## p, the per-unit influence function matrix psi (n x length(a)), and a length-n
## multiplier vector xi; monotonized via sort() as in findqx
mbootp=function(p, psi, xi){
  psib=(t(psi)-p)%*%xi
  pb_raw=p+psib/length(xi)
  pb=pmax(pmin(1, pb_raw),0)
  sort(pb)
}

## pointwise CI via the (unscaled) percentile/reflection method
pciper=function(datab, est, alpha){
  stopifnot(is.matrix(datab), nrow(datab) == length(est))
  stopifnot(alpha > 0, alpha < 1)
  std_b=rep(1, times = length(est))
  result=scale(t(datab), center = est, scale = std_b)
  z_up=apply(result, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  z_low=apply(result, 2, quantile, probs = alpha/2, na.rm = TRUE)
  pci=cbind(est - std_b*z_up, est - std_b*z_low)
  pci=apply(pci, 1, sort)
  t(pci)
}

## uniform confidence band via the rescaled quantile spread and bootstrap
## Kolmogorov-Smirnov max-t-statistic (Chernozhukov, Fernandez-Val, Melly 2013;
## Section 2.5, eq. 22-23 of Hsu, Huber and Yen 2026)
sciqs=function(datab, est, q = 0.1, alpha){
  stopifnot(is.matrix(datab), nrow(datab) == length(est))
  stopifnot(q>0, q<0.5)
  stopifnot(alpha > 0, alpha < 1)
  q_up=apply(datab, 1, quantile, probs = 1 - q, na.rm = TRUE)
  q_low=apply(datab, 1, quantile, probs = q, na.rm = TRUE)
  std_b=sqrt((q_up - q_low)^2/(qnorm(1 - q) - qnorm(q))^2)
  stopifnot("all std_b values are NA or zero - check bootstrap draws" = any(std_b != 0, na.rm = TRUE))
  if(any(is.na(std_b))) warning(sum(is.na(std_b)), " grid point(s) have NA std_b; check bootstrap draws.")
  std_b[std_b==0 | is.na(std_b)]=min(std_b[std_b!=0 & !is.na(std_b)])
  result=abs(scale(t(datab), center = est, scale = std_b))
  result_z_max=apply(result, 1, max, na.rm = TRUE)
  zx=as.numeric(quantile(result_z_max, probs = 1 - alpha))
  sci=cbind(est - std_b*zx, est + std_b*zx)
  sci=apply(sci,1,sort)
  t(sci)
}

## lower partial mean (retained from function_med_qte.R; not used internally
## by medqteDML but kept available as a utility)
lpm=function(q, x){
  indx=x<=q
  mean(x*indx)
}

## -----------------------------------------------------------------------
## plassoselectDD1: fold-level lasso covariate selection for the D|X and
## D|M,X models (ports plasso_model_selection_D_D1.R)
## -----------------------------------------------------------------------
plassoselectDD1=function(Ds, Ms, xDs){
  xD1s=as.matrix(cbind(Ms, xDs))
  modD=rlassologit(x = xDs, y = Ds)
  selD=!as.vector(modD$beta==0)
  modD1=rlassologit(x = xD1s, y = Ds)
  selD1=!as.vector(modD1$beta==0)
  selD1[-1]|selD                      # selX0: union of selected covariates, excluding M's own coefficient
}

## -----------------------------------------------------------------------
## plassoDD1: post-lasso fit of pi(X)=P(D=1|X) and p(S)=P(D=1|M,X) on the
## fold-level selected covariates selX, predicted onto the test fold
## (ports plasso_D_D1.R). Also returns which test-fold units had either
## propensity winsorized (see medqteplassofit for why winsorizing, not
## exclusion, is used here).
## -----------------------------------------------------------------------
plassoDD1=function(Ds, Dp, Ms, Mp, xDs, xDp, selX, trim=0.02){
  data_D=data.frame(D=Ds, xDs[,selX,drop=FALSE])
  new_data_D=data.frame(D=Dp, xDp[,selX,drop=FALSE])
  if(sum(selX)!=0){
    modD=glm(D~., data=data_D, family=binomial(link="logit"))
    pD=predict(modD, newdata=new_data_D, type="response")
  } else {
    modD=glm(Ds~1, family=binomial(link="logit"))
    pD=rep(predict(modD, type="response")[1], length(Dp))
  }
  trimD=(pD<trim)|(pD>1-trim)
  pD=pmin(pmax(pD,trim),1-trim)
  pD1=pD; pD0=1-pD
  
  data_D1=data.frame(D=Ds, M=Ms, xDs[,selX,drop=FALSE])
  new_data_D1=data.frame(D=Dp, M=Mp, xDp[,selX,drop=FALSE])
  modD1=glm(D~., data=data_D1, family=binomial(link="logit"))
  pD1mi=predict(modD1, newdata=new_data_D1, type="response")
  trimD1mi=(pD1mi<trim)|(pD1mi>1-trim)
  pD1mi=pmin(pmax(pD1mi,trim),1-trim)
  pD0mi=1-pD1mi
  
  list(pD1=pD1, pD0=pD0, pD1mi=pD1mi, pD0mi=pD0mi, trimmed=trimD|trimD1mi)
}

## -----------------------------------------------------------------------
## qteinputimpute: splits the (Ic_k)-fold training data into two halves and
## fits the first-stage distribution regression model P(Y<=a|D,M,X), then
## predicts it out-of-sample onto the other half (for the nested regression
## in qteregimpute) and onto the test fold (used directly in the score)
## (ports qte_plasso_input_regression_imputation.R). The direction of the
## split alternates across folds, generalizing the original script's
## hardcoded "s>=4" (for kfold=5) to s > ceiling(kfold/2) for general kfold.
## -----------------------------------------------------------------------
qteinputimpute=function(indYs, Ds, Ms, xDs, Dp, Mp, xDp, selX, s, kfold){
  n_s=nrow(xDs)
  if(s <= ceiling(kfold/2)){
    half_size=ceiling(n_s/2)
    idx1=1:half_size; idx2=(half_size+1):n_s
  } else {
    half_size=floor(n_s/2)
    idx2=1:half_size; idx1=(half_size+1):n_s
  }
  indYs1=indYs[idx1]; Ds1=Ds[idx1]; Ms1=Ms[idx1]; xDs1=xDs[idx1,,drop=FALSE]
  Ds2=Ds[idx2]; Ms2=Ms[idx2]; xDs2=xDs[idx2,,drop=FALSE]
  
  new_data_DMX1=data.frame(D=Ds2, M=Ms2, MD=Ms2*Ds2, xDs2[,selX,drop=FALSE])
  new_data_DMX =data.frame(D=Dp,  M=Mp,  MD=Mp*Dp,   xDp[,selX,drop=FALSE])
  
  data_indYs1=data.frame(y=indYs1, D=Ds1, M=Ms1, MD=Ms1*Ds1, xDs1[,selX,drop=FALSE])
  mod_indYs1=glm(y~., data=data_indYs1, family=binomial(link="logit"))
  
  ## imputed pseudo-outcomes for the nested regression (out-of-sample, half 2)
  mod_data_Y1Mi=new_data_DMX1; mod_data_Y1Mi$D=1; mod_data_Y1Mi$MD=mod_data_Y1Mi$D*mod_data_Y1Mi$M
  py1mi_hat=as.numeric(predict(mod_indYs1, newdata=mod_data_Y1Mi, type="response"))
  mod_data_Y0Mi=new_data_DMX1; mod_data_Y0Mi$D=0; mod_data_Y0Mi$MD=mod_data_Y0Mi$D*mod_data_Y0Mi$M
  py0mi_hat=as.numeric(predict(mod_indYs1, newdata=mod_data_Y0Mi, type="response"))
  
  ## counterfactual F(a|d,M,X) predicted directly onto the test fold, for the score
  new_data_Y1Mi=new_data_DMX; new_data_Y1Mi$D=1; new_data_Y1Mi$MD=new_data_Y1Mi$D*new_data_Y1Mi$M
  py1mi=as.numeric(predict(mod_indYs1, newdata=new_data_Y1Mi, type="response"))
  new_data_Y0Mi=new_data_DMX; new_data_Y0Mi$D=0; new_data_Y0Mi$MD=new_data_Y0Mi$D*new_data_Y0Mi$M
  py0mi=as.numeric(predict(mod_indYs1, newdata=new_data_Y0Mi, type="response"))
  
  list(py1mi_hat=py1mi_hat, py0mi_hat=py0mi_hat, py1mi=py1mi, py0mi=py0mi, Ds2=Ds2, xDs2=xDs2)
}

## -----------------------------------------------------------------------
## qteregimpute: nested regression imputation of g_{d,d',a}(X) (ports
## qte_plasso_regression_imputation.R)
## -----------------------------------------------------------------------
qteregimpute=function(py1mi_hat, py0mi_hat, Ds2, xDs2, indYs, Ds, xDs, Dp, xDp, selX){
  data_py1mi=data.frame(y=py1mi_hat, D=Ds2, xDs2[,selX,drop=FALSE])
  data_py0mi=data.frame(y=py0mi_hat, D=Ds2, xDs2[,selX,drop=FALSE])
  data_YDX  =data.frame(y=indYs,     D=Ds,  xDs[,selX,drop=FALSE])
  
  new_data_DX=data.frame(D=Dp, xDp[,selX,drop=FALSE])
  new_data_py1mi=new_data_DX; new_data_py1mi$D=1
  new_data_py0mi=new_data_DX; new_data_py0mi$D=0
  
  mod_py1mi=lm(y~., data=data_py1mi)
  mod_py0mi=lm(y~., data=data_py0mi)
  mod_YDX=glm(y~., data=data_YDX, family=binomial(link="logit"))
  
  g11a_x=as.numeric(predict(mod_YDX, newdata=new_data_py1mi, type="response"))    # = F(a|1,X), same-world
  g10a_x=as.numeric(predict(mod_py1mi, newdata=new_data_py0mi))                    # cross-world, nested
  g10a_x=pmin(pmax(g10a_x,0),1)
  g00a_x=as.numeric(predict(mod_YDX, newdata=new_data_py0mi, type="response"))    # = F(a|0,X), same-world
  g01a_x=as.numeric(predict(mod_py0mi, newdata=new_data_py1mi))                    # cross-world, nested
  g01a_x=pmin(pmax(g01a_x,0),1)
  
  list(g11a_x=g11a_x, g10a_x=g10a_x, g00a_x=g00a_x, g01a_x=g01a_x)
}

## -----------------------------------------------------------------------
## medqteplassofit: the k-fold x grid-point double loop, i.e. the
## cross-fitting estimator theta_hat_{d,d',a} of (15)-(16) and the
## per-unit efficient influence functions needed for the multiplier
## bootstrap (ports med_qte_plasso_db_fit.R). Returns the cdf estimates of
## the four potential outcomes across the grid a, plus their influence
## function matrices (n x length(a)) in the ORIGINAL row order of the data.
## -----------------------------------------------------------------------
medqteplassofit=function(y, d, m, x, a, kfold=5, trim=0.02){
  n=length(y)
  set.seed(1); idx=sample(n)     # randomize row order before sequential fold-blocking
  ind_pred=kfoldind(n=n, kfold=kfold)
  
  xD=as.matrix(x)
  if(is.null(colnames(xD))) colnames(xD)=paste0("X",1:ncol(xD))   # avoid data.frame()'s deparse-based
  xY=as.matrix(cbind(d, m, m*d, xD))                              # column naming when exactly one column
  # is selected out of an unnamed matrix
  
  px11=matrix(0,kfold,length(a)); px10=matrix(0,kfold,length(a))
  px00=matrix(0,kfold,length(a)); px01=matrix(0,kfold,length(a))
  if11=matrix(0,n,length(a)); if10=matrix(0,n,length(a))
  if00=matrix(0,n,length(a)); if01=matrix(0,n,length(a))
  trimmed=matrix(FALSE,n,length(a))    # TRUE where a unit's D|X or D|M,X propensity was winsorized
  
  for(s in 1:nrow(ind_pred)){
    ind_st=ind_pred$ind_st[s]; ind_end=ind_pred$ind_end[s]
    tesample=idx[ind_st:ind_end]; trsample=idx[-c(ind_st:ind_end)]
    
    Ys=y[trsample]; Yp=y[tesample]
    Ms=m[trsample]; Mp=m[tesample]
    Ds=d[trsample]; Dp=d[tesample]
    xDs=xD[trsample,,drop=FALSE]; xDp=xD[tesample,,drop=FALSE]
    xYs=xY[trsample,,drop=FALSE]
    
    selX0=plassoselectDD1(Ds=Ds, Ms=Ms, xDs=xDs)
    
    ind_d1=as.numeric(Dp==1); ind_d0=as.numeric(Dp==0)
    
    for(i in 1:length(a)){
      indYs=as.numeric(Ys<=a[i])
      indYp=as.numeric(Yp<=a[i])
      
      mod_indYs=rlassologit(x=xYs, y=indYs)
      sel_indYs=!as.vector(mod_indYs$beta==0)
      selX=sel_indYs[-c(1:3)]|selX0    # union with the fold-level D|X, D|M,X selection
      
      pD_fit=plassoDD1(Ds=Ds, Dp=Dp, Ms=Ms, Mp=Mp, xDs=xDs, xDp=xDp, selX=selX, trim=trim)
      inp=qteinputimpute(indYs=indYs, Ds=Ds, Ms=Ms, xDs=xDs, Dp=Dp, Mp=Mp, xDp=xDp, selX=selX, s=s, kfold=kfold)
      gimp=qteregimpute(py1mi_hat=inp$py1mi_hat, py0mi_hat=inp$py0mi_hat, Ds2=inp$Ds2, xDs2=inp$xDs2,
                        indYs=indYs, Ds=Ds, xDs=xDs, Dp=Dp, xDp=xDp, selX=selX)
      
      tr11=ind_d1/pD_fit$pD1*(indYp-gimp$g11a_x)+gimp$g11a_x
      tr10=ind_d1/(pD_fit$pD0*pD_fit$pD1mi)*pD_fit$pD0mi*(indYp-inp$py1mi)
      tr10=tr10+ind_d0/pD_fit$pD0*(inp$py1mi-gimp$g10a_x)+gimp$g10a_x
      tr00=ind_d0/pD_fit$pD0*(indYp-gimp$g00a_x)+gimp$g00a_x
      tr01=ind_d0/(pD_fit$pD1*pD_fit$pD0mi)*pD_fit$pD1mi*(indYp-inp$py0mi)
      tr01=tr01+ind_d1/pD_fit$pD1*(inp$py0mi-gimp$g01a_x)+gimp$g01a_x
      
      px11[s,i]=pmax(pmin(1,mean(tr11)),0); px10[s,i]=pmax(pmin(1,mean(tr10)),0)
      px00[s,i]=pmax(pmin(1,mean(tr00)),0); px01[s,i]=pmax(pmin(1,mean(tr01)),0)
      
      if11[tesample,i]=tr11; if10[tesample,i]=tr10
      if00[tesample,i]=tr00; if01[tesample,i]=tr01
      trimmed[tesample,i]=pD_fit$trimmed
    }
  }
  
  list(px11=apply(px11,2,mean,na.rm=TRUE), px10=apply(px10,2,mean,na.rm=TRUE),
       px00=apply(px00,2,mean,na.rm=TRUE), px01=apply(px01,2,mean,na.rm=TRUE),
       if11=if11, if10=if10, if00=if00, if01=if01, trimmed=trimmed)
}


# Internal helpers backing detectIV(). Not exported.
#
# test_conditional_independence() handles both binary and non-binary
# conditioning variables (the "instrument" Z in detectIV). For a binary Z it applies the score of eq. (7) of Apfel et
# al. (2025); for a non-binary Z it applies the multi-partition score of
# eq. (18) / Appendix C of the same paper. The two are structurally the same
# formula — eq. (7) is simply the L=1 special case of eq. (18) — so a single
# function covers both.
#
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# MLmean2
#
# K-fold cross-fitted estimation of a single conditional mean E[y | x] using
# MLfunct2(). Used by detectIV for the first-stage PLR nuisance regressions.
# Automatically detects binary outcomes and switches to binomial family.
#
# Returns an n-vector of out-of-fold predictions.
# ---------------------------------------------------------------------------
MLmean2 <- function(y, x, MLmethod = "lasso", k = 3, seed = 123) {
  ybin     <- as.integer(length(unique(y)) == 2 && min(y) == 0 && max(y) == 1)
  x        <- data.frame(x)
  n        <- length(y)
  stepsize <- ceiling(n / k)
  
  set.seed(seed)
  idx   <- sample(n, replace = FALSE)
  preds <- numeric(n)
  
  for (i in seq_len(k)) {
    te_pos  <- ((i - 1) * stepsize + 1) : min(i * stepsize, n)
    te_idx  <- idx[te_pos]
    tr_idx  <- idx[!(idx %in% te_idx)]
    
    fit        <- MLfunct2(y = y[tr_idx], x = x[tr_idx, , drop = FALSE],
                           MLmethod = MLmethod, ybin = ybin)
    preds[te_idx] <- predict(fit, newdata = x[te_idx, , drop = FALSE],
                             onlySL = TRUE)$pred
  }
  
  preds
}


# ---------------------------------------------------------------------------
# estimate_nuisance_cf
#
# K-fold cross-fitted estimation of three nuisance functions for a single
# binary partition indicator `Z_bin` (0/1):
#
#   mu1(cond_vars) = E[Y | cond_vars, Z_bin = 1]   (fitted on the Z_bin=1 rows)
#   mu0(cond_vars) = E[Y | cond_vars, Z_bin = 0]   (fitted on the Z_bin=0 rows)
#   p(cond_vars)   = Pr(Z_bin = 1 | cond_vars)
#
# cond_vars = cbind(D, X) when called from test_conditional_independence.
# For detectIV binary: Z_bin is the candidate IV bin.
# For detectIV continuous: called once per partition bin l inside
#   test_conditional_independence, with Z_bin = 1(Z in Z_l).
# ---------------------------------------------------------------------------
estimate_nuisance_cf <- function(Y, Z_bin, cond_vars,
                                 MLmethod = "lasso", K = 5) {
  n <- length(Y)
  
  mu1_hat <- numeric(n)
  mu0_hat <- numeric(n)
  p_hat   <- numeric(n)
  
  cv <- as.matrix(cond_vars)
  colnames(cv) <- paste0("V", seq_len(ncol(cv)))
  
  fold_id <- sample(rep(seq_len(K), length.out = n))
  
  for (kk in seq_len(K)) {
    test_idx  <- which(fold_id == kk)
    train_idx <- which(fold_id != kk)
    
    Y_tr   <- Y[train_idx]
    Z_tr   <- Z_bin[train_idx]
    cv_tr  <- cv[train_idx, , drop = FALSE]
    cv_te  <- cv[test_idx,  , drop = FALSE]
    
    fit_p <- MLfunct2(y = Z_tr, x = cv_tr, MLmethod = MLmethod, ybin = 1)
    p_hat[test_idx] <- stats::predict(fit_p,
                                      newdata = data.frame(cv_te),
                                      onlySL  = TRUE)$pred
    
    in1_tr <- which(Z_tr == 1)
    fit_mu1 <- MLfunct2(y = Y_tr[in1_tr],
                        x = cv_tr[in1_tr, , drop = FALSE],
                        MLmethod = MLmethod, ybin = 0)
    mu1_hat[test_idx] <- stats::predict(fit_mu1,
                                        newdata = data.frame(cv_te),
                                        onlySL  = TRUE)$pred
    
    in0_tr <- which(Z_tr == 0)
    fit_mu0 <- MLfunct2(y = Y_tr[in0_tr],
                        x = cv_tr[in0_tr, , drop = FALSE],
                        MLmethod = MLmethod, ybin = 0)
    mu0_hat[test_idx] <- stats::predict(fit_mu0,
                                        newdata = data.frame(cv_te),
                                        onlySL  = TRUE)$pred
  }
  
  list(mu1 = mu1_hat, mu0 = mu0_hat, p = p_hat)
}


# ---------------------------------------------------------------------------
# test_conditional_independence
#
# Tests H0: E[Y | Z_bin=1, D, X] = E[Y | Z_bin=0, D, X] using the doubly
# robust score of Apfel et al. (2025).
#
# When Z_bin is binary (2 unique values), the score of eq. (7) is applied
# directly.  When Z_bin is non-binary it is first discretised into L bins and
# the multi-partition score of eq. (18) / Appendix C is applied, summing the
# squared and linear score terms over all L bins.  Equation (7) is the L=1
# special case of eq. (18), so both code paths implement the same formula.
#
# Arguments
# ---------
# Y         : numeric outcome vector
# Z         : candidate IV — binary or non-binary 
# D         : treatment 
# X         : numeric matrix of additional covariates
# MLmethod  : passed to MLfunct2 / estimate_nuisance_cf
# K         : number of cross-fitting folds
# epsilon   : trimming threshold for propensity scores
# L         : number of partition bins for non-binary Z (ignored when Z is
#             binary).  Variables with fewer than 15 unique values use their
#             natural categories regardless of L.
# ---------------------------------------------------------------------------
test_conditional_independence <- function(Y, Z, D, X,
                                          MLmethod = "lasso", K = 5,
                                          epsilon = 0.05, L = 4) {
  
  cond_vars <- cbind(D, X)   # conditioning set shared by mu and p
  
  num_unique <- length(unique(Z))
  
  ## ---- Binary Z: eq. (7) — single propensity score, single pair of mus ----
  if (num_unique == 2) {
    
    Z_bin <- as.numeric(Z == max(Z))   # ensure {0,1} coding
    
    nuisance <- estimate_nuisance_cf(Y, Z_bin, cond_vars,
                                     MLmethod = MLmethod, K = K)
    mu1   <- nuisance$mu1
    mu0   <- nuisance$mu0
    p_hat <- nuisance$p
    
    keep <- which(p_hat >= epsilon & p_hat <= 1 - epsilon)
    
    delta_mu <- mu1 - mu0
    r1       <- Y - mu1
    r0       <- Y - mu0
    ipw_term <- (r1 * Z_bin / p_hat) - (r0 * (1 - Z_bin) / (1 - p_hat))
    
    psi <- delta_mu^2 + 2 * delta_mu * ipw_term + delta_mu + ipw_term
    
    theta_hat <- mean(psi[keep])
    sigma_hat <- stats::sd(psi[keep])
    n_eff     <- length(keep)
    se        <- sigma_hat / sqrt(n_eff)
    p_value   <- 2 * stats::pnorm(-abs(theta_hat / se))
    
    return(list(teststat = theta_hat, se = se, pval = p_value, n_eff = n_eff))
  }
  
  ## ---- Non-binary Z: eq. (18) — sum over L partition bins -----------------
  
  ## Build partition indicators
  if (num_unique < 15) {
    ## Discrete variable with few levels: use natural categories as bins
    unique_vals <- sort(unique(Z))
    Z_part      <- factor(match(Z, unique_vals))
    L_eff       <- length(unique_vals)
  } else {
    ## Continuous: bin by quantiles
    qbreaks <- unique(quantile(Z, probs = seq(0, 1, length.out = L + 1),
                               na.rm = TRUE))
    if (length(qbreaks) < 3) {
      warning("test_conditional_independence: too few unique quantile breaks; ",
              "falling back to median split.")
      qbreaks <- c(min(Z), stats::median(Z), max(Z))
    }
    Z_part <- factor(cut(Z, breaks = qbreaks,
                         include.lowest = TRUE, labels = FALSE))
    L_eff  <- length(levels(Z_part))
  }
  Z_bin <- stats::model.matrix(~ Z_part - 1)   # n x L_eff indicator matrix
  
  ## Estimate nuisance once per bin, store results, then trim and score.
  n          <- length(Y)
  nuisance_l <- vector("list", L_eff)
  
  for (l in seq_len(L_eff)) {
    nuisance_l[[l]] <- estimate_nuisance_cf(Y, Z_bin[, l], cond_vars,
                                            MLmethod = MLmethod, K = K)
  }
  
  ## Trim: drop observations with any bin propensity outside [epsilon, 1-epsilon]
  ps_mat <- sapply(nuisance_l, `[[`, "p")           # n x L_eff matrix
  keep   <- which(apply(ps_mat, 1,
                        function(r) all(r >= epsilon & r <= 1 - epsilon)))
  
  ## Accumulate score over bins
  psi <- numeric(n)
  for (l in seq_len(L_eff)) {
    in_bin_l <- Z_bin[, l]
    mu1_l    <- nuisance_l[[l]]$mu1
    mu0_l    <- nuisance_l[[l]]$mu0
    p_l      <- nuisance_l[[l]]$p
    
    delta_l <- mu1_l - mu0_l
    ipw_l   <- ((Y - mu1_l) * in_bin_l       / p_l) -
      ((Y - mu0_l) * (1 - in_bin_l) / (1 - p_l))
    
    psi <- psi + delta_l^2 + 2 * delta_l * ipw_l + delta_l + ipw_l
  }
  
  theta_hat <- mean(psi[keep])
  sigma_hat <- stats::sd(psi[keep])
  n_eff     <- length(keep)
  se        <- sigma_hat / sqrt(n_eff)
  p_value   <- 2 * stats::pnorm(-abs(theta_hat / se))
  
  list(teststat = theta_hat, se = se, pval = p_value, n_eff = n_eff)
}


#####
# helpers for CATE heterogeneity
#####

catehet_MLmean <- function(y, x, MLmethod = "lasso", k = 3) {
  n    <- length(y)
  ybin <- 1 * (length(unique(y)) == 2 && all(sort(unique(y)) == c(0,1)))  # check if binary outcome
  X    <- data.frame(x)
  out  <- numeric(n)
  idx  <- sample(n); foldsize <- ceiling(n / k)
  # cross-fitting procedure that splits sample in training and testing data
  for (fold in seq_len(k)) {
    te <- idx[((fold - 1) * foldsize + 1):min(fold * foldsize, n)]
    tr <- setdiff(idx, te)
    fit <- MLfunct(y = y[tr], x = X[tr, , drop = FALSE], MLmethod = MLmethod, ybin = ybin)
    out[te] <- predict(fit, newdata = X[te, , drop = FALSE], onlySL = TRUE)$pred
  }
  out
}

# encode the study indicator: discrete z gives one study per value with the last one dropped,
# continuous z is cut into L bins of equal size, all of which are kept
catehet_zdummies <- function(z, L = 4) {
  categorical <- is.factor(z) || is.character(z) || (is.numeric(z) && all(z == floor(z), na.rm = TRUE))
  if (categorical) {
    Zp <- droplevels(factor(z))
    L  <- nlevels(Zp) - 1
    if (L < 1) stop("'z' must take at least two distinct values.")
    Zb <- model.matrix(~ Zp - 1)[, 1:L, drop = FALSE]
    colnames(Zb) <- levels(Zp)[1:L]
  }
  if (!categorical) {
    Zp <- cut(z, breaks = quantile(z, seq(0, 1, length.out = L + 1), na.rm = TRUE), include.lowest = TRUE)
    Zb <- model.matrix(~ Zp - 1)
    colnames(Zb) <- levels(Zp)
  }
  Zb
}

# score function: for each study, Delta is the difference between its CATE and the CATE in the
# remaining studies, R the doubly robust correction, and psi = sum (Delta^2 + 2*Delta*R + Delta + R)
catehet_psi <- function(y, mu11_mat, mu10_mat, mu01_mat, mu00_mat,
                        p_1_z, p_0_z, p_1_zm, p_0_zm,
                        d, Zb, trim = 0.05, normalized = TRUE) {
  n <- length(y); L <- ncol(Zb)
  psi <- numeric(n); w_full <- numeric(n); n_eff_z <- numeric(L)
  
  for (l in seq_len(L)) {
    mu11 <- mu11_mat[, l]; mu10 <- mu10_mat[, l]; mu01 <- mu01_mat[, l]; mu00 <- mu00_mat[, l]
    p11  <- p_1_z[, l];    p01  <- p_0_z[, l];    p11m <- p_1_zm[, l];   p01m <- p_0_zm[, l]
    
    # observations with sufficient overlap in all four treatment-by-study cells
    in_support <- which(!is.na(p11) & !is.na(p01) & !is.na(p11m) & !is.na(p01m) &
                          p11 >= trim & p11 <= 1 - trim & p01  >= trim & p01  <= 1 - trim &
                          p11m >= trim & p11m <= 1 - trim & p01m >= trim & p01m <= 1 - trim)
    if (length(in_support) == 0) next
    
    I_1_z  <- in_support[ Zb[in_support, l] == 1 & d[in_support] == 1 ]
    I_0_z  <- in_support[ Zb[in_support, l] == 1 & d[in_support] == 0 ]
    I_1_zm <- in_support[ Zb[in_support, l] == 0 & d[in_support] == 1 ]
    I_0_zm <- in_support[ Zb[in_support, l] == 0 & d[in_support] == 0 ]
    
    w1_r  <- if (length(I_1_z))  1 / p11[I_1_z]   else numeric(0)
    w0_r  <- if (length(I_0_z))  1 / p01[I_0_z]   else numeric(0)
    w1m_r <- if (length(I_1_zm)) 1 / p11m[I_1_zm] else numeric(0)
    w0m_r <- if (length(I_0_zm)) 1 / p01m[I_0_zm] else numeric(0)
    
    # normalize the inverse probability weights within each cell
    if (normalized) {
      w1  <- if (length(w1_r))  w1_r  / mean(w1_r)  else numeric(0)
      w0  <- if (length(w0_r))  w0_r  / mean(w0_r)  else numeric(0)
      w1m <- if (length(w1m_r)) w1m_r / mean(w1m_r) else numeric(0)
      w0m <- if (length(w0m_r)) w0m_r / mean(w0m_r) else numeric(0)
    }
    if (!normalized) { w1 <- w1_r; w0 <- w0_r; w1m <- w1m_r; w0m <- w0m_r }
    
    r11 <- y - mu11; r10 <- y - mu10; r01 <- y - mu01; r00 <- y - mu00
    
    R_l <- numeric(n)
    if (length(I_1_z))  R_l[I_1_z]  <- (r11[I_1_z]  * w1)
    if (length(I_0_z))  R_l[I_0_z]  <- R_l[I_0_z]  - (r10[I_0_z]  * w0)
    if (length(I_1_zm)) R_l[I_1_zm] <- R_l[I_1_zm] - (r01[I_1_zm] * w1m)
    if (length(I_0_zm)) R_l[I_0_zm] <- R_l[I_0_zm] + (r00[I_0_zm] * w0m)
    
    Delta_l <- mu11 - mu10 - (mu01 - mu00)
    psi <- psi + Delta_l^2 + 2 * Delta_l * R_l + Delta_l + R_l
    
    if (length(I_1_z))  w_full[I_1_z]  <- abs(w1)
    if (length(I_0_z))  w_full[I_0_z]  <- abs(w0)
    if (length(I_1_zm)) w_full[I_1_zm] <- abs(w1m)
    if (length(I_0_zm)) w_full[I_0_zm] <- abs(w0m)
    
    all_w <- c(w1, w0, w1m, w0m)   # effective sample size for this study
    n_eff_z[l] <- if (length(all_w)) sum(all_w)^2 / sum(all_w^2) else 0
  }
  
  names(n_eff_z) <- colnames(Zb)
  list(psi = psi, w_full = w_full, n_eff_z = n_eff_z)
}

# =============================================================================
# SuperLearner-free replacement for the nuisance learners in causalweight
#
# Replicates, line by line, what SuperLearner::SuperLearner() and
# SuperLearner:::predict.SuperLearner(..., onlySL = TRUE) do in SuperLearner
# 2.0-42 (CRAN, 2026-09-14) for the call pattern used in causalweight:
#   SuperLearner(Y = y, X = x, family = gaussian()/binomial(), SL.library = lib)
# i.e. all defaults: method.NNLS, V = 10 shuffled (unstratified) folds,
# no id, unit observation weights, no screening, saveFitLibrary = TRUE.
#
# The learner wrappers below are verbatim ports of SL.glmnet, SL.ranger,
# SL.xgboost, SL.ksvm, SL.glm, SL.lm and their predict.SL.* methods with
# their default tuning parameters. Given the same RNG state, fitted values,
# NNLS weights, predictions and the RNG state after the call are identical
# to SuperLearner's.
#
# Remaining dependencies: glmnet, ranger, xgboost, kernlab, nnls
# (nnls is only needed for the ensemble weights; see note at .cw_nnls_coef).
# =============================================================================


# -----------------------------------------------------------------------------
# Learner wrappers: fit on (Y, X), predict on newX; return list(pred, fit).
# Same signature and semantics as the SL.* wrappers.
# -----------------------------------------------------------------------------
.cw_learners <- list(
  
  # SL.glmnet: alpha = 1, nfolds = 10, nlambda = 100, useMin = TRUE, loss = "deviance"
  glmnet = function(Y, X, newX, family, obsWeights) {
    if (!is.matrix(X)) {
      X    <- stats::model.matrix(~ -1 + ., X)
      newX <- stats::model.matrix(~ -1 + ., newX)
    }
    fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                               lambda = NULL, type.measure = "deviance",
                               nfolds = 10, family = family$family,
                               alpha = 1, nlambda = 100)
    pred <- stats::predict(fitCV, newx = newX, type = "response", s = "lambda.min")
    fit  <- structure(list(object = fitCV), class = "cw_glmnet")
    list(pred = pred, fit = fit)
  },
  
  # SL.ranger: num.trees = 500, mtry = floor(sqrt(p)), min.node.size = 5 (gaussian) / 1,
  # replace = TRUE, sample.fraction = 1, probability forest for binomial, 1 thread
  ranger = function(Y, X, newX, family, obsWeights) {
    if (family$family == "binomial") Y <- as.factor(Y)
    if (is.matrix(X)) X <- data.frame(X)
    fit <- ranger::ranger(`_Y` ~ ., data = cbind("_Y" = Y, X),
                          num.trees = 500,
                          mtry = floor(sqrt(ncol(X))),
                          min.node.size = ifelse(family$family == "gaussian", 5, 1),
                          replace = TRUE,
                          sample.fraction = 1,
                          case.weights = obsWeights,
                          write.forest = TRUE,
                          probability = family$family == "binomial",
                          num.threads = 1,
                          verbose = TRUE)
    pred <- stats::predict(fit, data = newX)$predictions   # draws one runif(), as in SL.ranger
    if (family$family == "binomial") pred <- pred[, "1"]
    list(pred = pred, fit = structure(list(object = fit, verbose = TRUE), class = "cw_ranger"))
  },
  
  # SL.xgboost: ntrees = 1000, max_depth = 4, shrinkage = 0.1, minobspernode = 10, 1 thread.
  # Two code paths, exactly as SuperLearner (new xgboost() interface for xgboost > 3.0).
  xgboost = function(Y, X, newX, family, obsWeights) {
    if (utils::packageVersion("xgboost") < "0.6") stop("xgboost version >= 0.6 required")
    if (utils::packageVersion("xgboost") > "3.0") {
      if (family$family == "gaussian") {
        model <- xgboost::xgboost(x = X, y = Y, weights = obsWeights,
                                  objective = "reg:squarederror", nrounds = 1000,
                                  max_depth = 4, min_child_weight = 10,
                                  learning_rate = 0.1, verbosity = 0, nthreads = 1)
      }
      if (family$family == "binomial") {
        model <- xgboost::xgboost(x = X, y = as.factor(Y), weights = obsWeights,
                                  objective = "binary:logistic", nrounds = 1000,
                                  max_depth = 4, min_child_weight = 10,
                                  learning_rate = 0.1, verbosity = 0, nthreads = 1,
                                  eval_metric = "logloss")
      }
      pred <- stats::predict(model, newdata = newX, type = "response")
      return(list(pred = pred, fit = structure(list(object = model), class = "cw_xgboost")))
    }
    # legacy path (xgboost <= 3.0); can be dropped if DESCRIPTION requires xgboost (> 3.0)
    if (!is.matrix(X)) X <- stats::model.matrix(~ . - 1, X)
    xgmat <- xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
    args  <- list(data = xgmat, nrounds = 1000, max_depth = 4, min_child_weight = 10,
                  eta = 0.1, verbose = 0, nthread = 1, params = list(), save_period = NULL)
    if (family$family == "gaussian") {
      args$objective <- if (utils::packageVersion("xgboost") >= "1.1.1.1") "reg:squarederror" else "reg:linear"
    }
    if (family$family == "binomial") {
      args$objective   <- "binary:logistic"
      args$eval_metric <- "logloss"
    }
    model <- do.call(xgboost::xgboost, args)
    if (!is.matrix(newX)) newX <- stats::model.matrix(~ . - 1, newX)
    pred <- stats::predict(model, newdata = newX)
    list(pred = pred, fit = structure(list(object = model), class = "cw_xgboost"))
  },
  
  # SL.ksvm: rbfdot kernel, kpar = "automatic" (sigest, random), C = 1, epsilon = 0.1,
  # scaled = TRUE, Platt-scaled probabilities (prob.model) for binomial
  ksvm = function(Y, X, newX, family, obsWeights) {
    if (!is.matrix(X)) {
      X <- stats::model.matrix(~ ., data = X)
      X <- X[, -1]                       # sic: no drop = FALSE in SL.ksvm
    }
    if (family$family == "binomial") {
      Y <- as.factor(Y); predict_type <- "probabilities"
    } else {
      predict_type <- "response"
    }
    model <- kernlab::ksvm(X, Y, scaled = TRUE, type = NULL, kernel = "rbfdot",
                           kpar = "automatic", C = 1, nu = 0.2, epsilon = 0.1,
                           prob.model = family$family == "binomial",
                           class.weights = NULL)
    if (!is.matrix(newX)) {
      newX <- stats::model.matrix(~ ., data = newX)
      newX <- newX[, -1, drop = FALSE]
    }
    pred <- kernlab::predict(model, newX, predict_type)
    if (family$family == "binomial") pred <- pred[, 2]
    list(pred = pred, fit = structure(list(object = model, family = family), class = "cw_ksvm"))
  },
  
  # SL.glm (binary outcomes under MLmethod = "parametric")
  glm = function(Y, X, newX, family, obsWeights) {
    if (is.matrix(X)) X <- as.data.frame(X)
    fit.glm <- stats::glm(Y ~ ., data = X, family = family, weights = obsWeights, model = TRUE)
    if (is.matrix(newX)) newX <- as.data.frame(newX)
    pred <- stats::predict(fit.glm, newdata = newX, type = "response")
    list(pred = pred, fit = structure(list(object = fit.glm), class = "cw_glm"))
  },
  
  # SL.lm (continuous outcomes under MLmethod = "parametric")
  lm = function(Y, X, newX, family, obsWeights) {
    if (is.matrix(X)) X <- as.data.frame(X)
    fit <- stats::lm(Y ~ ., data = X, weights = obsWeights, model = TRUE)
    if (is.matrix(newX)) newX <- as.data.frame(newX)
    pred <- stats::predict(fit, newdata = newX, type = "response")
    if (family$family == "binomial") pred <- pmin(pmax(pred, 0), 1)
    list(pred = pred, fit = structure(list(object = fit, family = family), class = "cw_lm"))
  }
)


# -----------------------------------------------------------------------------
# Prediction from a stored learner fit (ports of predict.SL.*)
# -----------------------------------------------------------------------------
.cw_predict_learner <- function(fit, newdata, family) {
  switch(class(fit)[1],
         
         cw_glmnet = {
           if (!is.matrix(newdata)) newdata <- stats::model.matrix(~ -1 + ., newdata)
           original_cols <- rownames(fit$object$glmnet.fit$beta)
           extra_cols <- setdiff(colnames(newdata), original_cols)
           if (length(extra_cols) > 0) {
             warning(paste("Removing extra columns in prediction data:",
                           paste(extra_cols, collapse = ", ")))
             newdata <- newdata[, !colnames(newdata) %in% extra_cols, drop = FALSE]
           }
           missing_cols <- setdiff(original_cols, colnames(newdata))
           if (length(missing_cols) > 0) {
             warning(paste("Adding missing columns in prediction data:",
                           paste(missing_cols, collapse = ", ")))
             new_cols <- matrix(0, nrow = nrow(newdata), ncol = length(missing_cols))
             colnames(new_cols) <- missing_cols
             newdata <- cbind(newdata, new_cols)[, original_cols, drop = FALSE]
           }
           stats::predict(fit$object, newx = newdata, type = "response", s = "lambda.min")
         },
         
         cw_ranger = {
           pred <- stats::predict(fit$object, data = newdata, verbose = fit$verbose,
                                  num.threads = 1)$predictions
           if (family$family == "binomial") pred <- pred[, "1"]
           pred
         },
         
         cw_xgboost = {
           if (utils::packageVersion("xgboost") > "3.0") {
             stats::predict(fit$object, newdata = newdata)
           } else {
             if (!is.matrix(newdata)) newdata <- stats::model.matrix(~ . - 1, newdata)
             stats::predict(fit$object, newdata = newdata)
           }
         },
         
         cw_ksvm = {
           if (!is.matrix(newdata)) {
             newdata <- stats::model.matrix(~ ., data = newdata)
             newdata <- newdata[, -1, drop = FALSE]
           }
           predict_type <- if (family$family == "binomial") "probabilities" else "response"
           pred <- kernlab::predict(fit$object, newdata, predict_type, coupler = "minpair")
           if (family$family == "binomial") pred <- pred[, 2]
           pred
         },
         
         cw_glm = {
           if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
           stats::predict(fit$object, newdata = newdata, type = "response")
         },
         
         cw_lm = {
           if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
           pred <- stats::predict(fit$object, newdata = newdata, type = "response")
           if (fit$family$family == "binomial") pred <- pmin(pmax(pred, 0), 1)
           pred
         },
         
         stop("unknown learner class: ", class(fit)[1])
  )
}


# -----------------------------------------------------------------------------
# method.NNLS$computeCoef. nnls::nnls (Lawson-Hanson) is what SuperLearner
# uses; keeping it guarantees bit-identical ensemble weights. For a single
# learner the weight is 1 if the NNLS coefficient is positive and 0 otherwise.
# -----------------------------------------------------------------------------
.cw_nnls_coef <- function(Z, Y, obsWeights) {
  fit.nnls <- nnls::nnls(sqrt(obsWeights) * Z, sqrt(obsWeights) * Y)
  initCoef <- fit.nnls$x
  initCoef[is.na(initCoef)] <- 0.0
  if (sum(initCoef) > 0) {
    initCoef / sum(initCoef)
  } else {
    warning("All algorithms have zero weight", call. = FALSE)
    initCoef
  }
}


# -----------------------------------------------------------------------------
# cw_superlearner(): replica of SuperLearner() for the causalweight call pattern
# -----------------------------------------------------------------------------
cw_superlearner <- function(Y, X, family = stats::gaussian(), learners, V = 10L) {
  N <- dim(X)[1L]
  k <- length(learners)
  if (sum(is.na(X)) > 0 | sum(is.na(Y)) > 0) {
    stop("missing data is currently not supported. Check Y, X, and newX for missing values")
  }
  if (!is.numeric(Y)) stop("the outcome Y must be a numeric vector")
  obsWeights <- rep(1, N)
  
  # CVFolds() with shuffle = TRUE, stratifyCV = FALSE, id = NULL
  validRows <- split(sample(1:N), rep(1:V, length = N))
  
  # cross-validated library predictions (folds outer, learners inner, as in SuperLearner)
  Z <- matrix(NA, N, k)
  for (valid in validRows) {
    for (s in seq_len(k)) {
      testAlg <- try(.cw_learners[[learners[s]]](Y = Y[-valid],
                                                 X = X[-valid, , drop = FALSE],
                                                 newX = X[valid, , drop = FALSE],
                                                 family = family,
                                                 obsWeights = obsWeights[-valid]))
      if (inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", learners[s],
                      "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      } else {
        Z[valid, s] <- testAlg$pred
      }
    }
  }
  errorsInCVLibrary <- apply(Z, 2, anyNA)
  if (sum(errorsInCVLibrary) > 0) Z[, errorsInCVLibrary] <- 0
  if (all(Z == 0)) stop("All algorithms dropped from library")
  coef <- .cw_nnls_coef(Z, Y, obsWeights)
  
  # full-sample fits; SuperLearner predicts on newX = X here, which matters for
  # the RNG stream (ranger's predict draws a seed), so this is kept
  fitLibrary <- vector("list", k)
  predY <- matrix(NA, nrow = N, ncol = k)
  for (s in seq_len(k)) {
    testAlg <- try(.cw_learners[[learners[s]]](Y = Y, X = X, newX = X, family = family,
                                               obsWeights = obsWeights))
    if (inherits(testAlg, "try-error")) {
      warning(paste("Error in algorithm", learners[s], " on full data",
                    "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
    } else {
      predY[, s] <- testAlg$pred
      fitLibrary[[s]] <- testAlg$fit
    }
  }
  errorsInLibrary <- apply(predY, 2, anyNA)
  if (sum(errorsInLibrary) > 0) {
    if (sum(coef[errorsInLibrary]) > 0) {
      warning("Re-running estimation of coefficients removing failed algorithm(s)")
      Z[, errorsInLibrary] <- 0
      if (all(Z == 0)) stop("All algorithms dropped from library")
      coef <- .cw_nnls_coef(Z, Y, obsWeights)
    } else {
      warning("Coefficients already 0 for all failed algorithm(s)")
    }
  }
  names(coef) <- learners
  structure(list(fitLibrary = fitLibrary, coef = coef, family = family,
                 learners = learners, cvRisk = colMeans((Z - Y)^2)),
            class = "cw_superlearner")
}


# -----------------------------------------------------------------------------
# predict method: replica of predict.SuperLearner(..., onlySL = TRUE).
# Register in NAMESPACE as S3method(predict, cw_superlearner) (roxygen:
# @exportS3Method stats::predict); then all existing call sites
#   predict(fit, newdata, onlySL = TRUE)$pred
# work unchanged. onlySL is accepted and ignored (always TRUE semantics).
# -----------------------------------------------------------------------------


#' Predictions from the internal nuisance learner
#'
#' @param object Fitted object of class \code{cw_superlearner}.
#' @param newdata Covariates for prediction.
#' @param onlySL Ignored; kept for compatibility with existing call sites.
#' @param ... Not used.
#' @return A list with elements \code{pred} and \code{library.predict}.
#' @method predict cw_superlearner
#' @export
#' @keywords internal
predict.cw_superlearner <- function(object, newdata, onlySL = TRUE, ...) {
  k <- length(object$learners)
  predY <- matrix(0, nrow = nrow(newdata), ncol = k)
  for (mm in which(object$coef > 0)) {
    predY[, mm] <- .cw_predict_learner(object$fitLibrary[[mm]], newdata, object$family)
  }
  coef <- object$coef
  if (sum(coef != 0) == 0) {
    warning("All metalearner coefficients are zero, predictions will all be equal to 0", call. = FALSE)
    pred <- rep(0, nrow(predY))
  } else {
    pred <- crossprod(t(predY[, coef != 0, drop = FALSE]), coef[coef != 0])
  }
  list(pred = pred, library.predict = predY)
}


# -----------------------------------------------------------------------------
# Replacements for MLfunct and MLfunct2 in functions.R
# -----------------------------------------------------------------------------
MLfunct <- function(y, x, d1 = NULL, d2 = NULL, MLmethod = "lasso", ybin = 0) {
  if (is.null(d1) == 0 & is.null(d2) == 0) { y <- y[d1 == 1 & d2 == 1]; x <- x[d1 == 1 & d2 == 1, ] }
  if (is.null(d1) == 0 & is.null(d2) == 1) { y <- y[d1 == 1]; x <- x[d1 == 1, ] }
  lib <- switch(MLmethod,
                lasso        = "glmnet",
                randomforest = "ranger",
                xgboost      = "xgboost",
                svm          = "ksvm",
                ensemble     = c("glmnet", "xgboost", "ksvm", "ranger"),
                parametric   = if (ybin == 1) "glm" else "lm",
                stop("unknown MLmethod: ", MLmethod))
  fam <- if (ybin == 1) stats::binomial() else stats::gaussian()
  cw_superlearner(Y = y, X = x, family = fam, learners = lib)
}

MLfunct2 <- function(y, x, MLmethod = "lasso", ybin = 0) {
  lib <- switch(MLmethod,
                lasso        = "glmnet",
                randomforest = "ranger",
                xgboost      = "xgboost",
                svm          = "ksvm",
                ensemble     = c("glmnet", "xgboost", "ranger"),
                parametric   = if (ybin == 1) "glm" else "lm",
                stop("unknown MLmethod: ", MLmethod))
  fam <- if (ybin == 1) stats::binomial() else stats::gaussian()
  cw_superlearner(Y = y, X = data.frame(x), family = fam, learners = lib)
}


