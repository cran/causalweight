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
    bsamples=matrix( ,boot,6)
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
    bsamples=matrix( ,boot,4)
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
    bsamples=matrix( ,boot,4)
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
  m.dist=npcdist(bws=c(bwreg,bwm), tydat=m, txdat=data.frame(zm,x), ckertype="gaussian")$condist
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
    bsamples=matrix( ,boot,8)
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
    bsamples=matrix( ,boot,6)
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
    bsamples=matrix( ,boot,8)
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
  means=matrix(,3,1)
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
    bsamples=matrix( ,boot,1)
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



hdmed=function(y,d,m,x,k=4, trim=0.05, order=1){
  if (order>1) {x=Generate.Powers(cbind(x),lambda=order); x=as.matrix(x,nrow(x),ncol(x))}
  stepsize=ceiling((1/k)*length(d))
  y1m0=c();y1m1=c();y0m0=c(); y0m1=c(); selall=c()
  # crossfitting procedure that splits sample in training an testing data
  for (i in 1:k){
    tesample=c(((i-1)*stepsize+1):(min((i)*stepsize,length(d))))
    dtr=d[-tesample]; dte=d[tesample]; ytr=y[-tesample]; yte=y[tesample]; ytr1= ytr[dtr==1]; ytr0= ytr[dtr==0];
    mtr=m[-tesample]; mte=m[tesample]; mtr1=mtr[dtr==1]; mtr0=mtr[dtr==0];
    if (is.null(ncol(x)) | ncol(x)==1) {
      xtr=x[-tesample]; xte=x[tesample]; xtr1=xtr[dtr==1]; xtr0=xtr[dtr==0]; xtr11=xtr[dtr==1 & mtr==1]; xtr10=xtr[dtr==1 & mtr==0]; xtr01=xtr[dtr==0 & mtr==1]; xtr00=xtr[dtr==0 & mtr==0]
    }
    if (is.null(ncol(x))==0 & ncol(x)>1) {
      xtr=x[-tesample,]; xte=x[tesample,]; xtr1=xtr[dtr==1,]; xtr0=xtr[dtr==0,]; xtr11=xtr[dtr==1 & mtr==1,]; xtr10=xtr[dtr==1 & mtr==0,]; xtr01=xtr[dtr==0 & mtr==1,]; xtr00=xtr[dtr==0 & mtr==0,]
    }
    ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)
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
    sel= 1*(((pdte*pm1te)>=trim) & ((1-pdte)>=trim)  & (pdte>=trim) &  (((1-pdte)*pm0te)>=trim)   )
    temp=dte*pm0te/(pdte*pm1te)*(yte-eymx1te)+(1-dte)/(1-pdte)*(eymx1te- eta10 )+eta10
    y1m0=c(y1m0, temp[sel==1])
    # predict score functions for E(Y(1,M(1))) in the test data
    temp=eyx1te + dte*(yte-eyx1te)/pdte
    y1m1=c(y1m1,temp[sel==1])
    # predict score functions for E(Y(0,M(1))) in the test data
    eta01=(eymx01te*pm1te+eymx00te*(1-pm1te))
    temp=(1-dte)*pm1te/((1-pdte)*pm0te)*(yte-eymx0te)+dte/pdte*(eymx0te- eta01 )+eta01
    y0m1=c(y0m1, temp[sel==1])
    # predict score functions for E(Y0,M(0)) in the test data
    temp=eyx0te + (1-dte)*(yte-eyx0te)/(1-pdte)
    y0m0=c(y0m0, temp[sel==1])
    selall=c(selall,sel)
  }
  # average over the crossfitting steps
  my1m1=mean(y1m1); my0m1=mean(y0m1); my1m0=mean(y1m0); my0m0=mean(y0m0)
  # compute effects
  tot=my1m1-my0m0; dir1=my1m1-my0m1; dir0=my1m0-my0m0; indir1=my1m1-my1m0; indir0=my0m1-my0m0;
  #compute variances
  vtot=mean((y1m1-y0m0-tot)^2); vdir1=mean((y1m1-y0m1-dir1)^2); vdir0=mean((y1m0-y0m0-dir0)^2);
  vindir1=mean((y1m1-y1m0-indir1)^2); vindir0=mean((y0m1-y0m0-indir0)^2); vcontrol=mean((y0m0-my0m0)^2)
  c(tot, dir1, dir0, indir1, indir0, my0m0, vtot, vdir1, vdir0, vindir1, vindir0, vcontrol, sum(selall))
}

# function for mediation with high dimensional covariates based on Bayes rule
hdmedalt=function(y,d,m,x, trim=0.05, order=1, fewsplits=FALSE){
  #generate higher order terms for lasso
    xm=cbind(x,m)
    if (order>1) {x=Generate.Powers(cbind(x),lambda=order); x=as.matrix(x,nrow(x),ncol(x))}
    if (order>1) xm=Generate.Powers(xm,lambda=order); xm=as.matrix(xm,nrow(xm),ncol(xm))

  stepsize=ceiling((1/4)*length(d))
  y1m0=c();y1m1=c();y0m0=c(); y0m1=c(); selall=c()
  # crossfitting procedure that splits sample in training an testing data
  sample1=c(1:stepsize); sample2=c( (stepsize+1):(2*stepsize)); sample3=c( (2*stepsize+1):(3*stepsize));
  sample4=c(((3)*stepsize+1):(min((4)*stepsize,length(d))))
  for (i in 1:4){
    if (i==1) {tesample=sample1; musample=sample2; psample=sample3; deltasample=sample4}
    if (i==2) {tesample=sample4; musample=sample1; psample=sample2; deltasample=sample3}
    if (i==3) {tesample=sample3; musample=sample4; psample=sample1; deltasample=sample2}
    if (i==4) {tesample=sample2; musample=sample3; psample=sample4; deltasample=sample1}
    trsample=c(musample,psample,deltasample); dte=d[tesample]; yte=y[tesample]
    # in case that fewsplits is one, psample and musample are merged
    if (fewsplits==1){musample=c(musample,psample);psample=musample}
      x=as.matrix(x,nrow(x),ncol(x)); xm=as.matrix(xm,nrow(xm),ncol(xm))
      ybin=1*(length(unique(y))==2 & min(y)==0 & max(y)==1)
      pmx=rlassologit(d[psample]~xm[psample,])
      # predict Pr(D=1|M,X) in test data
      pmxte=predict(pmx, xm[tesample,], type="response")
      # predict Pr(D=1|M,X) in delta sample
      pmxtrte=predict(pmx, xm[deltasample,], type="response")
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
      #compute E(Y|M,X,D=1)*(1-Pr(D=1|M,X)) based on predictions in delta sample
      weymx1trte=eymx1trte*(1-pmxtrte)
      # fit E[E(Y|M,X,D=1)*(1-Pr(D=1|M,X)|X] in delta sample
      regweymx1=rlasso(weymx1trte~x[deltasample,])
      # predict E[E(Y|M,X,D=1)*(1-Pr(D=1|M,X)|X] in the test data
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
      #compute E(Y|M,X,D=0)*Pr(D=1|M,X) based on predictions in second training data
      weymx0trte=eymx0trte*(pmxtrte)
      # fit E[E(Y|M,X,D=0)*Pr(D=1|M,X)|X] in delta sample
      regweymx0=rlasso(weymx0trte~x[deltasample,])
      # predict E[E(Y|M,X,D=0)*Pr(D=1|M,X)|X] in the test data
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
    # predict E(Y0,M(1)) in the test data
    temp=((1-dte)*pmxte/((1-pmxte)*pxte)*(yte-eymx0te)+dte/pxte*(eymx0te-regweymx0te/pxte)+regweymx0te/pxte)
    y0m1=c(y0m1,temp[sel==1])
    # predict E(Y0,M(0)) in the test data
    temp=(eyx0te + (1-dte)*(yte-eyx0te)/(1-pxte))
    y0m0=c(y0m0,temp[sel==1])
    # predict E(Y1,M(0)) in the test data
    temp=(dte*(1-pmxte)/(pmxte*(1-pxte))*(yte-eymx1te)+(1-dte)/(1-pxte)*(eymx1te-regweymx1te/(1-pxte))+regweymx1te/(1-pxte))
    y1m0=c(y1m0,temp[sel==1])
    # predict E(Y1,M(1)) in the test data
    temp=(eyx1te + dte*(yte-eyx1te)/pxte)
    y1m1=c(y1m1,temp[sel==1])
    # collect selection dummies
    selall=c(selall,sel)
  }
  # average over the crossfitting steps
  my1m1=mean(y1m1); my0m1=mean(y0m1); my1m0=mean(y1m0); my0m0=mean(y0m0)
  # compute effects
  tot=my1m1-my0m0; dir1=my1m1-my0m1; dir0=my1m0-my0m0; indir1=my1m1-my1m0; indir0=my0m1-my0m0;
  #compute variances
  vtot=mean((y1m1-y0m0-tot)^2); vdir1=mean((y1m1-y0m1-dir1)^2); vdir0=mean((y1m0-y0m0-dir0)^2);
  vindir1=mean((y1m1-y1m0-indir1)^2); vindir0=mean((y0m1-y0m0-indir0)^2); vcontrol=mean((y0m0-my0m0)^2)
  # report effects, mean of Y(0,M(0)), variances, number of non-trimmed observations
  c(tot, dir1, dir0, indir1, indir0, my0m0, vtot, vdir1, vdir0, vindir1, vindir0, vcontrol, sum(selall))
}

