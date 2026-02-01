alittleArt <-
function(dat,z,x=NULL,pr=NULL,xm=NULL,near=NULL,fine=NULL,
               xinteger=NULL,xbalance=NULL,ncontrols=1,rnd=2,solver="rlemon",
               min.penalty=c(10,1,0.05),pr.penalty=c(2,5,25,250),
               near.penalty=1000,fine.penalty=50,integer.penalty=20){
  #
  # Check input
  #
  if (is.null(x)&is.null(pr)) stop("If pr is NULL, then x cannot be NULL.")
  if ((!is.null(x))&(!is.null(pr)))
    warning("Because pr is not NULL, x is not used to estimate a propensity score.")
  stopifnot((solver=="rlemon")|(solver=="rrelaxiv"))
  stopifnot(is.matrix(dat)|is.data.frame(dat))
  # Check the treatment assignment vector z
  stopifnot(is.vector(z))
  stopifnot(all((z==0)|(z==1)))
  stopifnot(length(z)==(dim(dat)[1]))
  # Did the user supply the propensity score (or another score)?
  if (!is.null(x)) {
    stopifnot(is.matrix(x)|is.data.frame(x))
    stopifnot(length(z)==(dim(x)[1]))
    x<-as.matrix(x)
    if (!is.numeric(x)) stop("The columns of x must be numeric.")
  }
  else stopifnot (is.vector(pr)&(length(pr)==length(z)))
  # Check the matrix xm
  if (!is.null(xm)){
    stopifnot(is.matrix(xm)|is.data.frame(xm))
    stopifnot(length(z)==(dim(xm)[1]))
    stopifnot(is.numeric(xm))
  }
  # Check near and near.penalty
  if (!is.null(near)) {
    stopifnot(is.matrix(near)|is.data.frame(near)|is.vector(near))
    stopifnot(is.numeric(near))
    if (is.vector(near)) {
      stopifnot(length(near)==length(z))
      stopifnot(length(near.penalty)==1)
    }
    else {
      stopifnot((dim(near)[1])==length(z))
      stopifnot((length(near.penalty)==1)|(length(near.penalty)==(dim(near)[2])))
      if (length(near.penalty)==1)
        near.penalty<-rep(near.penalty,(dim(near)[2]))
    }
    stopifnot(all(near.penalty>=0))
  }
  # Check fine and fine.penalty
  if (!is.null(fine)) {
    stopifnot(is.matrix(fine)|is.data.frame(fine)|is.vector(fine))
    stopifnot(is.numeric(fine))
    if (is.vector(fine)) {
      stopifnot(length(fine)==length(z))
      stopifnot(length(fine.penalty)==1)
    }
    else {
      stopifnot((dim(fine)[1])==length(z))
      stopifnot((length(fine.penalty)==1)|(length(fine.penalty)==(dim(fine)[2])))
      if (length(fine.penalty)==1) {
        fine.penalty<-rep(fine.penalty,(dim(fine)[2]))
      }
      stopifnot(all(fine.penalty>=0))
    }
  }
  # Check xinteger and integer.penalty
  if (!is.null(xinteger)) {
    stopifnot(is.matrix(xinteger)|is.data.frame(xinteger)|is.vector(xinteger))
    stopifnot(is.numeric(xinteger))
    if (is.vector(xinteger)) {
      stopifnot(length(xinteger)==length(z))
      stopifnot(length(integer.penalty)==1)
    }
    else {
      stopifnot((dim(xinteger)[1])==length(z))
      stopifnot((length(integer.penalty)==1)|(length(integer.penalty)==(dim(xinteger)[2])))
      if (length(integer.penalty)==1) {
        integer.penalty<-rep(integer.penalty,(dim(xinteger)[2]))
      }
      stopifnot(all(integer.penalty>=0))
    }
  }
  # Check xbalance
  if (!is.null(xbalance)) {
    stopifnot(is.matrix(xbalance)|is.data.frame(xbalance)|is.vector(xbalance))
    stopifnot(is.numeric(xbalance))
    if (is.vector(xbalance)) xbalance<-matrix(xbalance,length(xbalance),1)
    stopifnot(length(z)==(dim(xbalance)[1]))
  }
  else if (!is.null(x)) xbalance<-x
  else if (is.null(x)) xbalance<-matrix(pr,length(pr),1)
  # Check the number of controls, ncontrols
  stopifnot(is.vector(ncontrols)&(length(ncontrols)==1))
  stopifnot(ncontrols>=1)
  stopifnot(ncontrols==round(ncontrols))
  if ((sum(z)*ncontrols)>sum(1-z)){
    stop(paste("You have ",sum(z)," treated and ",sum(1-z),
               " controls, so you cannot match each treated to ",
               ncontrols," controls.",sep=""))
  }
  # Check the propensity penalties
  stopifnot(length(pr.penalty)==4)
  stopifnot(all(pr.penalty>=0))
  stopifnot(pr.penalty[1]<=pr.penalty[2])
  stopifnot(pr.penalty[3]<=pr.penalty[4])
  # Check the min.penalty
  stopifnot(length(min.penalty)==3)
  stopifnot(all(min.penalty>=0))
  stopifnot((min.penalty[3]>0)&(min.penalty[3]<1))
  #
  # Set up match
  #
  if (is.null(rownames(dat))) rownames(dat)<-1:(dim(dat)[1])
  rndat<-rownames(dat)
  names(z)<-rownames(dat)
  if (is.null(pr)){
    constant<-rep(1,dim(x)[1])
    prop<-stats::glm.fit(cbind(constant,x),z,family=stats::binomial())
    pr<-prop$fitted.values
    dat<-cbind(dat)
  }
  pr6<-as.integer(cut(pr,stats::quantile(pr,c(0,1/6,2/6,3/6,4/6,5/6,1)),
                      include.lowest=TRUE))
  pr3<-as.integer(cut(pr,stats::quantile(pr,c(0,1/3,2/3,1)),
                      include.lowest=TRUE))
  left<-iTOS::startcost(z)
  right<-iTOS::startcost(z)
  if (!is.null(xm)) left<-iTOS::addMahal(left,z,xm)
  if (!is.null(near)) {
    if (is.vector(near)) left<-iTOS::addNearExact(left,z,near,
                                      penalty=near.penalty)
    else {
      for (j in 1:(dim(near)[2]))
        left<-iTOS::addNearExact(left,z,near[,j],penalty=near.penalty[j])
    }
  }
  left<-iTOS::addinteger(left,z,pr3,penalty=pr.penalty[2])
  left<-iTOS::addinteger(left,z,pr6,penalty=pr.penalty[1])
  right<-iTOS::addinteger(right,z,pr3,penalty=pr.penalty[4])
  right<-iTOS::addinteger(right,z,pr6,penalty=pr.penalty[3])
  if (!is.null(fine)) {
    if (is.vector(fine)) right<-iTOS::addNearExact(right,z,
                                          fine,penalty=fine.penalty)
    else {
      for (j in 1:(dim(fine)[2]))
        right<-iTOS::addNearExact(right,z,fine[,j],penalty=fine.penalty[j])
    }
  }
  # Add integer penalties to the network
  if (!is.null(xinteger)) {
    if (is.vector(xinteger)) right<-iTOS::addinteger(right,z,
                                xinteger,penalty=integer.penalty)
    else {
      for (j in 1:(dim(xinteger)[2]))
        right<-iTOS::addinteger(right,z,xinteger[,j],penalty=integer.penalty[j])
    }
  }
  # Add minimum propensity score penalties to the network
  # using controlcosts
  cc<-min.penalty[1]*(pr[z==0]<min(pr[z==1]))
  cc<-cc+min.penalty[2]*(pr[z==0]<=quantile(pr[z==1],min.penalty[3]))
  # Use the network to find the match
  m<-iTOS::makematch(dat,left,right,ncontrols=ncontrols,controlcosts=cc,solver=solver)
  propensityScore<-pr
  v<-cbind(propensityScore,xbalance)
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
