\name{PeriUnmatched}
\alias{PeriUnmatched}
\docType{data}
\title{
Unmatched Periodontal Disease Data
}
\description{
Unmatched data from NHANES 2009-2010, 2011-2012, 2013-2014 concerning smoking
and periodontal disease.
}
\usage{data("PeriUnmatched")}
\format{
  A data frame with 6255 observations on the following 11 variables.
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
  }
}
\details{
Measurements were made for up to 28 teeth, 14 upper, 14 lower, excluding 4 wisdom teeth. Pocket depth and loss of attachment are two complementary measures of the degree to which the gums have separated from the teeth; see Wei, Barker and Eke (2013). Pocket depth and loss of attachment are measured at six locations on each tooth, providing the tooth is present. A measurement at a location was taken to exhibit disease if it had either a loss of attachement >=4mm or a pocked depth >=4mm, so each tooth contributes six binary scores, up to 6x28=168 binary scores.  The variable pd is the percent of these binary scores indicating periodontal disease, 0 to 100 percent.

The data from three NHANES surveys (specifically 2009-2010, 2011-2012, and 2013-2014) contain periodontal data and are used as an example in Rosenbaum (2025).  The data from one survey, 2011-2012, were used in Rosenbaum (2016).
The example uses these unmatched data twice in artless() to create the fused match in Rosenbaum (2025).  The fused match combines some 1-to-1 matched pairs and some 1-to-4 matched sets based on the values of the propensity score.  The data are useful in learning about fused matching, but the example in the documentation for artless() should be used as the main example illustrating artless().
}

\note{An analysis of outcomes should take appropriate account of the matching; see the note in the documentation for PeriMatched.  Often, covariate balance is assessed by comparing the marginal distributions of covariates in treated and control groups after matching; however, some care is required when there are both 1-to-1 pairs and 1-to-4 sets.  One can assess covariate balance for the pairs, and separately assess covariate balance for the 1-to-4 sets.  Alternatively, one can measure covariate balance in the pairs and the 1-to-4 sets separately, perhaps taking the difference in means, and then take a weighted combination of the two differences in means for pairs and 1-to-4 sets, along the lines indicated by Pimentel, Yoon and Keele (2015).  However, one cannot assess covariate balance by pooling the two treated groups from pairs and 1-to-4 sets, pooling the two control groups from pairs and 1-to-4 sets, and comparing the two pooled groups.  In the example, there is exact matching for sex; however, most pairs are men and most 1-to-4 sets are women.  Pool the pairs and the 1-to-4 sets and the pooled control group has proportionately more women than the pooled treated group.  To see this, type:

data("PeriMatched")

tapply(PeriMatched$female,PeriMatched$grp3,mean)

tapply(PeriMatched$female,PeriMatched$grp2,mean)
}
\source{
US National Health and Nutrition Examination Survey (NHANES). https://www.cdc.gov/nchs/nhanes/
}
\references{
Pimentel, S. D., Yoon, F., & Keele, L. (2015) <doi:10.1002/sim.6593> Variable‐ratio matching with fine balance in a study of the Peer Health Exchange. Statistics in Medicine, 34(30), 4070-4082.

Rosenbaum, P. R. (2016) <doi:10.1214/16-AOAS942> Using Scheffe projections for multiple outcomes in an observational study of smoking and periondontal disease. Annals of Applied Statistics, 10, 1447-1471.

Rosenbaum, Paul R. (2025) A Design for Observational Studies in Which Some People Avoid Treatment.  Manuscript.

Tomar, S. L. and Asma, S. (2000). Smoking attributable periodontitis in the United States: Findings from NHANES III. J. Periodont. 71, 743-751.

Wei, L., Barker, L. and Eke, P. (2013). Array applications in determining periodontal disease measurement. SouthEast SAS User's Group. (SESUG2013) Paper CC-15, analytics.ncsu.edu/ sesug/2013/CC-15.pdf.
}
\examples{
# The code below creates the matched data, PeriMatched, from the unmatched
# data PeriUnmatched using the function artless() twice. Individuals
# with prop above 0.15 were matched in pairs.  Individuals with prop of at
# most 0.15 were matched in a 1-to-5 ratio.
data(PeriUnmatched)
\donttest{
# Controls matched for female, age, education, income
d0<-PeriUnmatched
prop<-stats::glm(d0$z~d0$female+d0$age+d0$educ+d0$income,family=binomial)$fitted
d0<-cbind(d0,prop)
rm(prop)

# Pair match for higher propensity individuals
d1<-d0[d0$prop>0.15,]
attach(d1)
ageFloor<-floor(age/10)
lowInc<-1*(income<2)
highInc<-1*(income>=4)
x<-cbind(female,age,educ,income)
xm<-cbind(age,educ,income)
near<-cbind(female,ageFloor)
age60<-1*(age>=60)
fine<-cbind(age60,noHS,lowInc,highInc,female)
# Match does the following: estimates a new propensity score in
# this subpopulation using the covariates in x, uses a
# Mahalanobis distance for the covariates in xm, performs near-exact
# matched for the covariates in near, and performs near-fine balancing
# of the covariates in near.  The solves rlemon is used because it is
# available in R, but rrelaxiv may be a better choice, though it
# requires a separate installation.
m<-artless(d1,z,x,xm=xm,near=near,fine=fine,solver="rlemon")
detach(d1)
# Some clean-up follows
rm(age60)
dm<-m$match
dm<-dm[!is.na(dm$mset),]
rm(x,xm,fine,near,d1,ageFloor,lowInc,highInc)
treated<-as.vector(rbind(dm$SEQN[dm$z==1],dm$SEQN[dm$z==1]))
dm<-cbind(dm,treated)
rm(treated)

# Now match 1-to-4 for low propensity individuals
d1<-d0[d0$prop<=0.15,]
attach(d1)
ageFloor<-floor(age/10)
lowInc<-1*(income<2)
highInc<-1*(income>=4)
x<-cbind(female,age,educ,income)
xm<-cbind(age,educ,income)
near<-cbind(female,ageFloor)
age60<-1*(age>=60)
fine<-cbind(age60,noHS,lowInc,highInc,female)
ncontrols<-4
# Match does the following: estimates a new propensity score in
# this subpopulation using the covariates in x, uses a
# Mahalanobis distance for the covariates in xm, performs near-exact
# matched for the covariates in near, and performs near-fine balancing
# of the covariates in near.  The solves rlemon is used because it is
# available in R, but rrelaxiv may be a better choice, though it
# requires a separate installation.
m1<-artless(d1,z,x,xm=xm,near=near,fine=fine,solver="rlemon",
                     ncontrols=ncontrols)
detach(d1)
# Some clean-up follows
rm(age60)
dm1<-m1$match
dm1<-dm1[!is.na(dm1$mset),]
rm(x,xm,fine,near,d1,ageFloor,lowInc,highInc)
treated1<-dm1$SEQN[dm1$z==1]
treated<-treated1
for (i in 1:(ncontrols)) treated<-rbind(treated,treated1)
treated<-as.vector(treated)
dm1<-cbind(dm1,treated)
rm(treated,treated1,i,ncontrols)

# Pool the two matched sames into one data.frame dm2
pair<-rep(1,dim(dm)[1])
dm<-cbind(dm,pair)
dm$mset<-as.integer(dm$mset)
pair<-rep(0,dim(dm1)[1])
dm1<-cbind(dm1,pair)
dm1$mset<-as.integer(dm1$mset)+max(dm$mset)
dm2<-rbind(dm1,dm)
rm(pair)
grp2<-factor(dm2$z,levels=1:0,labels=c("S","N"),ordered=TRUE)
grp3<-factor(dm2$pair,levels=c(1,0),labels=c("1-1","1-4"),ordered=TRUE):grp2
dm2<-cbind(dm2,grp2,grp3)
rm(grp2,grp3)

# There are 1212 pairs and 213 1-to-4 sets
table(table(dm2$mset))
# Check the balance tables separately for pairs and sets
# Pairs
m$balance
# 1-to-4 sets
m1$balance
}
}
\keyword{datasets}
