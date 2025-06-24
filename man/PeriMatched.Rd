\name{PeriMatched}
\alias{PeriMatched}
\docType{data}
\title{
Matched Periodontal Disease Data
}
\description{
Matched data from NHANES 2009-2010, 2011-2012, 2013-2014 concerning smoking
and periodontal disease.  The matched data were built from the unmatched
data in PeriUnmatched in this package.
}
\usage{data("PeriMatched")}
\format{
  A data frame with 3489 observations on the following 18 variables.
  \describe{
    \item{\code{SEQN}}{NHANES ID number}
    \item{\code{female}}{1=female, 0=male}
    \item{\code{age}}{Age in years, capped at 80 for confidentiality}
    \item{\code{ageFloor}}{Age decade = floor(age/10)}
    \item{\code{educ}}{Education as 1 to 5.  1 is less than 9th grade, 2 at least 9th grade with no high school degree, 3 is a high school degree, 4 is some college, such as a 2-year associates degree, 5 is at least a 4-year college degree.}
    \item{\code{noHS}}{No high school degree.  1 if educ is 1 or 2, 0 if educ is 3 or more}
    \item{\code{income}}{Ratio of family income to the poverty level, capped at 5 for confidenditality}
    \item{\code{nh}}{The specific NHANES survey.  A factor \code{nh0910} < \code{nh1112} < \code{nh1314}}
    \item{\code{cigsperday}}{Number of cigarettes smoked per day.  0 for nonsmokers.}
    \item{\code{z}}{Daily smoker.  1 indicates someone who smokes everyday.  0 indicates a never-smoker who smoked fewer than 100 cigarettes in their life.}
    \item{\code{pd}}{A percent indicating periodontal disease.  See details. }
    \item{\code{prop}}{A propensity score created in the example for PeriUnmatched.  This propensity score decided which smokers would have 1 control and which would have 5 controls.}
    \item{\code{pr}}{A second propensity score used to create matched pairs or matched 1-to-4 sets, after the split based on prop}
    \item{\code{mset}}{Indicator of the matched set, 1, 2, ..., 1425}
    \item{\code{treated}}{The SEQN for the smoker in this matched set.  Contains the same information as mset, but in a different form.}
    \item{\code{pair}}{1 for a matched pair, 0 for a 1-to-4 matched set}
    \item{\code{grp2}}{An ordered factor with the same information as z: S=daily smoker, N=never smoker. \code{S} < \code{N}}
    \item{\code{grp3}}{A factor with the joint information in pair and grp2.  \code{1-1:S} \code{1-1:N} \code{1-4:S} \code{1-4:N}}
}
}
\details{
Measurements were made for up to 28 teeth, 14 upper, 14 lower, excluding 4 wisdom teeth. Pocket depth and loss of attachment are two complementary measures of the degree to which the gums have separated from the teeth; see Wei, Barker and Eke (2013). Pocket depth and loss of attachment are measured at six locations on each tooth, providing the tooth is present. A measurement at a location was taken to exhibit disease if it had either a loss of attachement >=4mm or a pocked depth >=4mm, so each tooth contributes six binary scores, up to 6x28=168 binary scores.  The variable pd is the percent of these binary scores indicating periodontal disease, 0 to 100 percent.

The data from three NHANES surveys (specifically 2009-2010, 2011-2012, and 2013-2014) contain periodontal data and are used as an example in Rosenbaum (2025).  The data from one survey, 2011-2012, were used in Rosenbaum (2016).
The example replicates analyses from Rosenbaum (2025).
}
\note{All analyses below distinguish the 1-to-1 pairs and the 1-to-4 sets, even though the information they provide is often combined.  Alternatively, one can combine analyses of pairs and 1-to-4 sets using methods that take account of the matched blocks of variable sizes.  For instance, for continuous responses, one can use the methods in Rosenbaum (2007) as implemented in the R package sensitivitymv; see also Rosenbaum (2015).  For binary responses, one can use the methods in Rosenbaum and Small (2017) as implemented in the R package sensitivity2x2xk.

In contrast, some care is required in plots and descriptive statistics.  One can straightforwardly plot the pairs, then separately plot the 1-to-4 sets, and one can do the same with descriptive statistics.  Suppose, however, that one merges the two treated groups from pairs and 1-to-4 sets, and merges the two control groups from pairs and 1-to-4 sets; then marginal distributions of outcomes from the pooled treated and control groups are no longer comparable.  See Pimentel, Yoon and Keele (2015).  For instance, in the example, there is exact matching for sex; however, most pairs are men and most 1-to-4 sets are women.  Pool the pairs and the 1-to-4 sets and the pooled control group has proportionately more women than the pooled treated group.  To see this, type:

data("PeriMatched")

tapply(PeriMatched$female,PeriMatched$grp3,mean)

tapply(PeriMatched$female,PeriMatched$grp2,mean)

The simple, often enlightening, solution is to plot pairs and 1-to-4 sets in parallel but separately, and to do the same with descriptive statistics.}

\source{
US National Health and Nutrition Examination Survey (NHANES). https://www.cdc.gov/nchs/nhanes/
}
\references{
Pimentel, S. D., Yoon, F., & Keele, L. (2015) <doi:10.1002/sim.6593> Variable‐ratio matching with fine balance in a study of the Peer Health Exchange. Statistics in Medicine, 34(30), 4070-4082.

Rosenbaum, P. R. (2007) <doi:10.1111/j.1541-0420.2006.00717.x> Sensitivity analysis for m-estimates, tests, and confidence intervals in matched observational studies. Biometrics, 63(2), 456-464.

Rosenbaum, P. R. (2015) <doi:10.1353/obs.2015.0000> Two R packages for sensitivity analysis in observational studies. Observational Studies, 1(2), 1-17.  Available on-line at: muse.jhu.edu/article/793399/summary

Rosenbaum, P. R. (2016) <doi:10.1214/16-AOAS942> Using Scheffe projections for multiple outcomes in an observational study of smoking and periondontal disease. Annals of Applied Statistics, 10, 1447-1471.

Rosenbaum, P. R., & Small, D. S. (2017) <doi:10.1111/biom.12591> An adaptive Mantel–Haenszel test for sensitivity analysis in observational studies. Biometrics, 73(2), 422-430.

Rosenbaum, Paul R. (2025) A Design for Observational Studies in Which Some People Avoid Treatment.  Manuscript.

Tomar, S. L. and Asma, S. (2000). Smoking attributable periodontitis in the United States: Findings from NHANES III. J. Periodont. 71, 743-751.

Wei, L., Barker, L. and Eke, P. (2013). Array applications in determining periodontal disease measurement. SouthEast SAS User's Group. (SESUG2013) Paper CC-15, analytics.ncsu.edu/ sesug/2013/CC-15.pdf.
}
\examples{
data(PeriMatched)

# The analysis in Rosenbaum (2025) is replicated below
#
dm2<-PeriMatched
dm<-PeriMatched[PeriMatched$pair==1,]
dm1<-PeriMatched[PeriMatched$pair==0,]
pd1<-t(matrix(dm$pd,2,dim(dm)[1]/2))
pd4<-t(matrix(dm1$pd,5,dim(dm1)[1]/5))
dm2$mset<-as.integer(dm2$mset)

#
#  Make Figure 1
#
old.par <- par(no.readonly = TRUE)
par(mfrow=c(1,3))

boxplot(dm2$prop~dm2$grp3,names=c(expression(S[1]),expression(N[1]),
                                  expression(S[4]),expression(N[4])),
        las=1,sub="Left is 1-1,  Right is 1-4",cex.sub=.9,cex.axis=1,
        ylab="Propensity Score",xlab="(i) Propensity Score")
#axis(3,at=1:4,lab=round(tapply(dm2$prop,dm2$grp3,mean),2),cex.axis=1)
axis(3,at=1:4,lab=c("0.36","0.34","0.10","0.10"),cex.axis=1) # don't round 0.1

boxplot(dm2$educ~dm2$grp3,names=c(expression(S[1]),expression(N[1]),
                                  expression(S[4]),expression(N[4])),
        las=1,sub="Left is 1-1,  Right is 1-4",cex.sub=.9,cex.axis=1,
        ylab="Education: 1 is <9th, 3 is HS, 5 is BA",xlab="(ii) Education")
#axis(3,at=1:4,lab=round(tapply(dm2$educ,dm2$grp3,mean),1),cex.axis=1)
axis(3,at=1:4,lab=c("3.0","3.1","4.0","4.0"),cex.axis=1)

boxplot(dm2$income~dm2$grp3,names=c(expression(S[1]),expression(N[1]),
                                    expression(S[4]),expression(N[4])),
        las=1,sub="Left is 1-1,  Right is 1-4",cex.sub=.9,cex.axis=1,
        ylab="Income / (Poverty Level)",xlab="(iii) Income")
axis(3,at=1:4,lab=round(tapply(dm2$income,dm2$grp3,mean),1),cex.axis=1)

#
# Make Figure 2
#
par(mfrow=c(1,2))

boxplot(dm2$cigsperday~dm2$grp3,names=c(expression(S[1]),expression(N[1]),
                                        expression(S[4]),expression(N[4])),
        las=1,sub="Left is 1-1,  Right is 1-4",cex.sub=.9,cex.axis=1,
        ylab="Cigarettes Per Day",xlab="(i) Cigarettes Per Day")
axis(3,at=1:4,lab=round(tapply(dm2$cigsperday,dm2$grp3,mean),0),cex.axis=1)


boxplot(dm2$pd~dm2$grp3,names=c(expression(S[1]),expression(N[1]),
                                expression(S[4]),expression(N[4])),
        las=1,sub="Left is 1-1,  Right is 1-4",cex.sub=.9,cex.axis=1,
        ylab="Periodonal Disease",xlab="(ii) Periodontal Disease")
axis(3,at=1:4,lab=round(tapply(dm2$pd,dm2$grp3,mean),0),cex.axis=1)

#
# Make Table 1
#
tb<-NULL
N<-tapply(dm2$female,dm2$grp3,length)
tb<-cbind(tb,N)
rm(N)
Female<-tapply(dm2$female,dm2$grp3,mean)*100
tb<-cbind(tb,Female)
rm(Female)
Age<-tapply(dm2$age,dm2$grp3,mean)
tb<-cbind(tb,Age)
rm(Age)
Income<-tapply(dm2$income,dm2$grp3,mean)
tb<-cbind(tb,Income)
rm(Income)
Income10<-tapply(dm2$income,dm2$grp3,quantile,c(.1))
tb<-cbind(tb,Income10)
rm(Income10)
Income90<-tapply(dm2$income,dm2$grp3,quantile,c(.9))
tb<-cbind(tb,Income90)
rm(Income90)

Education25<-tapply(dm2$educ,dm2$grp3,quantile,c(.25))
tb<-cbind(tb,Education25)
rm(Education25)
Education50<-tapply(dm2$educ,dm2$grp3,quantile,c(.5))
tb<-cbind(tb,Education50)
rm(Education50)
Education75<-tapply(dm2$educ,dm2$grp3,quantile,c(.75))
tb<-cbind(tb,Education75)
rm(Education75)
PropensityMin<-tapply(dm2$prop,dm2$grp3,min)
tb<-cbind(tb,PropensityMin)
rm(PropensityMin)
Propensity<-tapply(dm2$prop,dm2$grp3,median)
tb<-cbind(tb,Propensity)
rm(Propensity)
PropensityMax<-tapply(dm2$prop,dm2$grp3,max)
tb<-cbind(tb,PropensityMax)
rm(PropensityMax)
xtable::xtable(tb,digits=c(NA,0,1,1,1,1,1,0,0,0,2,2,2))

addmargins(table(dm2$z,dm2$prop>.15))
#
# Make Table 2 regarding sensitivity analysis
#
gammas<-c(1:5,5.5,6)
ngamma<-length(gammas)
tabSen<-matrix(NA,ngamma,4)
colnames(tabSen)<-c("Pairs 1-1","Sets 1-4","Fisher","Truncated")
rownames(tabSen)<-gammas
for (i in 1:ngamma) tabSen[i,1]<-weightedRank::wgtRank(pd1,phi="u878",gamma=gammas[i])$pval
for (i in 1:ngamma) tabSen[i,2]<-weightedRank::wgtRank(pd4,phi="u878",gamma=gammas[i])$pval
for (i in 1:ngamma) {
  if (min(tabSen[i,1:2]==0)) tabSen[i,3:4]<-0
  else{
    tabSen[i,3]<-sensitivitymv::truncatedP(tabSen[i,1:2],trunc=1)
    tabSen[i,4]<-sensitivitymv::truncatedP(tabSen[i,1:2],trunc=0.2)
  }
}
# Table 2
xtable::xtable(t(tabSen),digits=4)

# Compare Table 2 to a sensitivity analysis for 1425 pairs-only
# by randomly selecting 1 of 4 controls from the 1-to-4 sets
set.seed(12345)
a<-sample(2:5,(dim(pd4)[1]),replace=TRUE)
pd4r<-rep(NA,(dim(pd4)[1]))
for (i in 1:(dim(pd4)[1])) pd4r[i] <- pd4[i,a[i]]
pd4r<-cbind(pd4[,1],pd4r)
rm(a)

weightedRank::wgtRank(rbind(pd1,pd4r),phi="u878",gamma=4.2)
weightedRank::wgtRank(rbind(pd1,pd4r),phi="quade",gamma=4)
weightedRank::wgtRank(rbind(pd1,pd4r),phi="quade",gamma=3)

#
# Make Table 3 regarding counterfactual risk
#
ctab<-table(dm2$pd>=20,dm2$grp3)
ctab<-ctab[2:1,]
ctab<-rbind(ctab,prop.table(ctab,2)[1,]*100)
ctab<-rbind(ctab,c(ctab[1,1]*ctab[2,2]/(ctab[1,2]*ctab[2,1]),
                   mantelhaen.test(table(dm$pd>=20,dm$z,dm$mset))$estimate,
                   ctab[1,3]*ctab[2,4]/(ctab[1,4]*ctab[2,3]),
                   mantelhaen.test(table(dm1$pd>=20,dm1$z,dm1$mset))$estimate))
xtable::xtable(ctab,digits=1)

#
#  Evidence factors analysis -- cigarettes per day
#
crosscutplot<-function (x, y, ct = 0.25, xlab = "", ylab = "", main = "",
                        ylim = NULL)
{
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(length(x) == length(y))
  stopifnot((ct > 0) & (ct <= 0.5))
  qx1 <- stats::quantile(x, ct)
  qx2 <- stats::quantile(x, 1 - ct)
  qy1 <- stats::quantile(y, ct)
  qy2 <- stats::quantile(y, 1 - ct)
  use <- ((x <= qx1) | (x >= qx2)) & ((y <= qy1) | (y >= qy2))
  if (is.null(ylim))
    graphics::plot(x, y, xlab = xlab, ylab = ylab, main = main,
                   type = "n",las=1,cex.lab=.9,cex.axis=.9,,cex.main=.9)
  else graphics::plot(x, y, xlab = xlab, ylab = ylab, ylim = ylim,,cex.main=.9,
                      main = main, type = "n",las=1,cex.lab=.9,cex.axis=.9)
  graphics::points(x[use], y[use], pch = 16,cex=.6)
  graphics::points(x[!use], y[!use], col = "gray", pch = 16,cex=.6)
  graphics::abline(h = c(qy1, qy2))
  graphics::abline(v = c(qx1, qx2))
}

dCigs1<-dm$cigsperday[dm$z==1]
dCigs4<-dm1$cigsperday[dm1$z==1]
dif1<-pd1[, 1] - pd1[, 2]
dif4<-pd4[,1]-apply(pd4[,2:5],1,median)
par(mfrow=c(1,2))
crosscutplot(dCigs1,dif1,xlab="Cigarettes per Day",ylim=c(-100,100),
             ylab="Periodontal Disease",main="1212 Pairs")
text(70,-80,paste("Odds Ratio =",round(89*135/(84*72),2)),cex=.7)
crosscutplot(dCigs4,dif4,xlab="Cigarettes per Day",ylim=c(-100,100),
             ylab="Periodontal Disease",
             main="213 Matched 1-to-4 Sets")
text(31,-80,paste("Odds Ratio =",round(28*18/(12*9),2)),cex=.7)
DOS2::crosscut(dCigs1,dif1)
DOS2::crosscut(dCigs4,dif4)
tb<-c(as.vector(DOS2::crosscut(dCigs1,dif1)$table),
      as.vector(DOS2::crosscut(dCigs4,dif4)$table))
tb<-array(tb,c(2,2,2))
sensitivity2x2xk::mh(tb,Gamma=1.6)
sensitivity2x2xk::mh(tb[,,1],Gamma=1.375)
sensitivity2x2xk::mh(tb[,,2],Gamma=1.7)


par(old.par)
rm(gammas,ngamma,crosscutplot,tb,i,tabSen,pd4r,old.par,ctab)
}
\keyword{datasets}
