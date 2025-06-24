artless<-function(dat,z,x,xm=NULL,near=NULL,fine=NULL,
                  ncontrols=1,rnd=2,solver="rlemon"){
  #
  # Check input
  #
  stopifnot((solver=="rlemon")|(solver=="rrelaxiv"))
  stopifnot(is.vector(z))
  stopifnot(all((z==0)|(z==1)))
  stopifnot(is.matrix(dat)|is.data.frame(dat))
  stopifnot(is.matrix(x)|is.data.frame(x))
  stopifnot(length(z)==(dim(dat)[1]))
  stopifnot(length(z)==(dim(x)[1]))
  if (!is.null(xm)){
    stopifnot(is.matrix(xm)|is.data.frame(xm))
    stopifnot(length(z)==(dim(xm)[1]))
    stopifnot(is.numeric(xm))
  }
  if (!is.null(near)) {
    stopifnot(is.matrix(near)|is.data.frame(near)|is.vector(near))
    stopifnot(is.numeric(near))
    if (is.vector(near)) stopifnot(length(near)==length(z))
    else stopifnot((dim(near)[1])==length(z))
  }
  if (!is.null(fine)) {
    stopifnot(is.matrix(fine)|is.data.frame(fine)|is.vector(fine))
    stopifnot(is.numeric(fine))
    if (is.vector(fine)) stopifnot(length(fine)==length(z))
    else stopifnot((dim(fine)[1])==length(z))
  }
  stopifnot(is.vector(ncontrols)&(length(ncontrols)==1))
  stopifnot(ncontrols>=1)
  stopifnot(ncontrols==round(ncontrols))
  if ((sum(z)*ncontrols)>sum(1-z)){
    stop(paste("You have ",sum(z)," treated and ",sum(1-z),
          " controls, so you cannot match each treated to ",
                ncontrols," controls.",sep=""))
  }
  #
  # Set up match
  #
  if (is.null(rownames(dat))) rownames(dat)<-1:(dim(dat)[1])
  rndat<-rownames(dat)
  names(z)<-rownames(dat)
  prop<-stats::glm.fit(x,z,family=stats::binomial())
  pr<-prop$fitted.values
  dat<-cbind(dat,pr)
  pr6<-as.integer(cut(pr,c(0,stats::quantile(pr,c(1/6,2/6,3/6,4/6,5/6)),1)))
  pr3<-as.integer(cut(pr,c(0,stats::quantile(pr,c(1/3,2/3)),1)))
  left<-iTOS::startcost(z)
  right<-iTOS::startcost(z)
  if (!is.null(xm)) left<-iTOS::addMahal(left,z,xm)
  if (!is.null(near)) {
    if (is.vector(near)) left<-iTOS::addNearExact(left,z,near)
    else {
      for (j in 1:(dim(near)[2])) left<-iTOS::addNearExact(left,z,near[,j])
    }
  }
  left<-iTOS::addinteger(left,z,pr3,penalty=5)
  left<-iTOS::addinteger(left,z,pr6,penalty=2)
  right<-iTOS::addinteger(right,z,pr3,penalty=250)
  right<-iTOS::addinteger(right,z,pr6,penalty=25)
  if (!is.null(fine)) {
    if (is.vector(fine)) right<-iTOS::addNearExact(right,z,fine,penalty=50)
    else {
      for (j in 1:(dim(fine)[2])) right<-iTOS::addNearExact(right,z,fine[,j],penalty=50)
    }
  }
  cc<-100*(pr[z==0]<min(pr[z==1]))
  m<-iTOS::makematch(dat,left,right,ncontrols=ncontrols,controlcosts=cc,solver=solver)
  propensityScore<-pr
  v<-cbind(propensityScore,x)
  if (!is.null(xm)) v<-cbind(v,xm)
  if (!is.null(fine)) v<-cbind(v,fine)
  if (!is.null(near)) v<-cbind(v,near)
  rownames(v)<-rownames(dat)
  treated<-apply(v[z==1,],2,mean)
  treateds<-apply(v[z==1,],2,stats::sd)
  controlB<-apply(v[z==0,],2,mean)
  controlBs<-apply(v[z==0,],2,stats::sd)
  control<-apply(v[(z==0)&(is.element(rownames(v),rownames(m))),],2,mean)
  st<-sqrt(((treateds^2)+(controlBs^2))/2)
  balance<-data.frame(treated,control,controlB,(treated-control)/st,
                      (treated-controlB)/st)
  for (j in 1:(dim(v)[2])) colnames(v)[j]<-paste(j,":",colnames(v)[j],sep="")
  rownames(balance)<-colnames(v)
  colnames(balance)<-c("Treated","Matched.Control","All.Control",
                       "StDif.After","StDif.Before")
  datu<-dat[!is.element(rownames(dat),rownames(m)),]
  mset<-rep(NA,dim(datu)[1])
  datu<-cbind(datu,mset)
  matched<-c(rep(TRUE,dim(m)[1]),rep(FALSE,dim(datu)[1]))
  m<-rbind(m,datu)
  m<-cbind(m,matched)
  list(match=m,balance=round(balance,rnd))
  }
