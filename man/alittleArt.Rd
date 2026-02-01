\name{alittleArt}
\alias{alittleArt}
\concept{Causal inference}
\concept{Propensity score}
\concept{Matching}
\concept{Observational studies}
\concept{Two-criteria matching}
\concept{Fine balance}
\concept{Optimal matching}
\concept{Multiple controls}
\title{
Artful Optimal Matching
}
\description{
Implements a simple version of multivariate matching using a propensity score, near-exact matching, near-fine balance, and robust Mahalanobis distance matching.  Provides fine control of the penalties used in matching.
}
\usage{
alittleArt(dat, z, x = NULL, pr = NULL, xm = NULL, near = NULL,
  fine = NULL, xinteger = NULL, xbalance = NULL, ncontrols = 1,
  rnd = 2, solver = "rlemon", min.penalty = c(10, 1, 0.05),
  pr.penalty = c(2, 5, 25, 250), near.penalty = 1000,
  fine.penalty = 50, integer.penalty = 20)
}
\arguments{
  \item{dat}{
A dataframe containing the data set that will be matched.  Let N be the
number of rows of dat.
}
  \item{z}{
A binary vector with N coordinates where z[i]=1 if the ith row of dat describes
a treated individual and z[i]=0 if the ith row of dat describes a control.
}
  \item{x}{
x is a numeric matrix with N rows.  If pr is NULL, then the covariates in x are used to estimate a propensity score using a linear logit model that predicts z from x.  An error will stop the program if pr and x are both NULL.  If neither pr nor x is NULL, then a harmless warning message will remind you that your propensity score pr was used in matching and x was not used to estimate the propensity score.  If xbalance is NULL, then the balance table will describe the covariates in x; so, those covariates should be continuous variables or binary variables that can be described by a mean or a proportion, not nominal categories.
}
  \item{pr}{
A vector with N coordinates containing an estimated propensity or similar quantity.  If pr is NULL, then the program estimates the propensity score; see the discussion of x above.
}
  \item{xm}{
xm is a numeric matrix with N rows.  The covariates in xm are used to define
a robust Mahalanobis distance between treated and control individuals.  The covariates in xm may be continuous variables like weight, integer covariates like number of rooms in a home, or binary variables; however, they should not be unordered nominal covariates like 1=New York, 2=Chicago, 3=London, 4=Tokyo.
}
  \item{near}{
A numeric vector of length N or a numeric matrix with N rows.  Each column of near should represent levels of a nominal covariate with two or a few levels.  The variables in near are used in near-exact matching.
}
  \item{fine}{
A numeric vector of length N or a numeric matrix with N rows.  Each column of fine should represent levels of a nominal covariate with two or a few levels.  The variables in fine are used in near-fine balancing.
}
  \item{xinteger}{
A numeric vector of length N or a numeric matrix with N rows.  Each column of xinteger should represent levels of an integer covariate with three or a few levels.  The variables in xinteger are used in near-fine balancing that prefers an imbalance from an adjacent category to an imbalance from a distant category. See the notes.
}
  \item{xbalance}{
If not NULL, xbalance is numeric vector of length N or a numeric matrix with N rows.  If xbalance is not NULL, then the balance table will describe the covariates in xbalance; so, those covariates should be continuous variables or binary variables that can be described by a mean or a proportion, not nominal categories.  See also the discussion of x above and the notes.
}
  \item{ncontrols}{
A positive integer.  ncontrols is the number of controls to be matched to each treated individual.
}
  \item{rnd}{
A nonnegative integer.  The balance table is rounded for display to rnd digits.
}
  \item{solver}{
Either "rlemon" or "rrelaxiv".  The rlemon solver is automatically available without special installation.  The rrelaxiv requires a special installation.  See the note.
}
  \item{min.penalty}{
A vector of three nonnegative coordinates.  The third coordinate must be strictly greater than zero and strictly less than one.  See the notes.
}
  \item{pr.penalty}{
A vector with four nonnegative coordinates that determine aspects of matching for the propensity score.  See the notes.
}
  \item{near.penalty}{
Either one nonnegative number of a vector of nonnegative numbers with one coordinate for each column of near.  See the notes.
}
  \item{fine.penalty}{
Either one nonnegative number of a vector of nonnegative numbers with one coordinate for each column of fine.  See the notes.
}
  \item{integer.penalty}{
Either one nonnegative number of a vector of nonnegative numbers with one coordinate for each column of xinteger.  See the notes.
}
}
\details{
This function builds a matched treated-control sample from an unmatched data set.  It asks you to designate roles for specific covariates, and it does the rest.  Unlike artlessV2(), the function alittleArt() gives you control over the penalties used in matching.  In particular, if in an initial match one covariate, say age, remains out of balance, then you can adjust a penalty specific to age to attempt to improve its balance.
}
\value{
\item{match }{A dataframe containing the matched data set.  match contains the rows of dat in a different order.  match adds two columns to dat, called mset and matched, which identify matched pairs or matched sets.  Specifically, matched is TRUE if a row is in the matched sample and is FALSE otherwise.  Rows of dat that are in the same matched set have the same value of mset.  The rows of match are sorted by mset with the treated individual before the matched controls.  The unmatched controls with matched=FALSE appear as the last rows of match. When you analyze the matched data, you will want to remove rows of match with matched==FALSE.}
\item{balance }{A matrix called the balance table.  The matrix has one row for each covariate in x.  It also has a first row for the propensity score.  There are five columns.
Column 1 is the mean of the covariate in the treated group.  Column 2 is the mean of the covariate in the matched control group.  Column 3 is the mean of the covariate among all controls prior to matching.  Column 4 is the difference between columns 1 and 2 divided by a pooled estimate of the standard deviation of the covariate before matching.  Column 5 is the difference between columns 1 and 3 divided by a pooled estimate of the standard deviation of the covariate before matching.  Notice that columns 4 and 5 have the same denominator, but different numerators.  Tom Love (2002) suggests a graphical display of this information.}
}
\references{
Bertsekas, D. P., Tseng, P. (1988) <doi:10.1007/BF02288322> The Relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.

Bertsekas, D. P. (1990) <doi:10.1287/inte.20.4.133> The auction algorithm for assignment and other network flow problems: A tutorial. Interfaces, 20(4), 133-149.

Bertsekas, D. P., Tseng, P. (1994) <http://web.mit.edu/dimitrib/www/Bertsekas_Tseng_RELAX4_!994.pdf> RELAX-IV: A Faster Version of the RELAX Code for Solving Minimum Cost Flow Problems.

Greifer, N. and Stuart, E.A., (2021). <doi:10.1093/epirev/mxab003> Matching methods for confounder adjustment: an addition to the epidemiologist’s toolbox. Epidemiologic Reviews, 43(1), pp.118-129.

Hansen, B. B. and Klopfer, S. O. (2006) <doi:10.1198/106186006X137047> "Optimal full matching and related designs via network flows". Journal of computational and Graphical Statistics, 15(3), 609-627. ('optmatch' package)

Hansen, B. B. (2007) <https://www.r-project.org/conferences/useR-2007/program/presentations/hansen.pdf> Flexible, optimal matching for observational studies. R News, 7, 18-24. ('optmatch' package)

Lee, K., Small, D.S. and Rosenbaum, P.R. (2018) <doi:10.1111/biom.12884> A powerful approach to the study of moderate effect modification in observational studies. Biometrics,
74:(4)1161-1170.

Love, Thomas E. (2002) Displaying covariate balance after adjustment for selection bias. Joint Statistical Meetings. Vol. 11.
\url{https://chrp.org/love/JSM_Aug11_TLove.pdf}

Niknam, B.A. and Zubizarreta, J.R. (2022). <10.1001/jama.2021.20555> Using cardinality matching to design balanced and representative samples for observational studies. JAMA, 327(2), pp.173-174.

Pimentel, S. D., Yoon, F., & Keele, L. (2015) <doi:10.1002/sim.6593> Variable‐ratio matching with fine balance in a study of the Peer Health Exchange. Statistics in Medicine, 34(30), 4070-4082.

Pimentel, S. D., Kelz, R. R., Silber, J. H. and Rosenbaum, P. R. (2015)
<doi:10.1080/01621459.2014.997879> Large, sparse optimal matching with refined covariate balance in an observational study of the health outcomes produced by new surgeons. Journal of the American Statistical Association, 110, 515-527.

Rosenbaum, P. R. and Rubin, D. B. (1985) <doi:10.1080/00031305.1985.10479383> Constructing a control group using multivariate matched sampling methods that incorporate the propensity score. The American Statistician, 39, 33-38.

Rosenbaum, P. R. (1989) <doi:10.1080/01621459.1989.10478868> Optimal matching for observational studies. Journal of the American Statistical Association, 84(408), 1024-1032.

Rosenbaum, P. R., Ross, R. N. and Silber, J. H. (2007) <doi:10.1198/016214506000001059> Minimum distance matched sampling with fine balance in an observational study of treatment for ovarian cancer. Journal of the American Statistical Association, 102, 75-83.

Rosenbaum, P. R. (2020a) <doi:10.1007/978-3-030-46405-9> Design of Observational Studies (2nd Edition). New York: Springer.

Rosenbaum, P. R. (2020b). <doi:10.1146/annurev-statistics-031219-041058> Modern algorithms for matching in observational studies. Annual Review of Statistics and Its Application, 7(1), 143-176.

Rosenbaum, P. R. and Zubizarreta, J. R. (2023). <doi:10.1201/9781003102670>
Optimization Techniques in Multivariate Matching. Handbook of Matching and Weighting Adjustments for Causal Inference, pp.63-86.  Boca Raton: FL: Chapman and Hall/CRC Press.

Rosenbaum, P. R. (2025) <doi:10.1007/978-3-031-90494-3> Introduction to the Theory of Observational Studies.  New York: Springer.

Rubin, D. B. (1980) <doi:10.2307/2529981> Bias reduction using Mahalanobis-metric matching. Biometrics, 36, 293-298.

Rubin, D. B. (2008) <doi:10.1214/08-AOAS187> For objective causal inference, design trumps analysis. Annals of Applied Statistics, 2, 808-840.

Stuart, E.A., (2010). <doi:10.1214/09-STS313> Matching methods for causal inference: A review and a look forward. Statistical Science, 25(1), 1-21.

Yang, D., Small, D. S., Silber, J. H. and Rosenbaum, P. R. (2012)
<doi:10.1111/j.1541-0420.2011.01691.x> Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes. Biometrics, 68, 628-636.

Yu, Ruoqi, and P. R. Rosenbaum. <doi:10.1111/biom.13098> Directional penalties for optimal matching in observational studies. Biometrics 75, no. 4 (2019): 1380-1390.

Yu, R., Silber, J. H., & Rosenbaum, P. R. (2020) <doi:10.1214/19-STS699> Matching methods for observational studies derived from large administrative databases. Statistical Science, 35(3), 338-355.

Yu, R. (2021) <doi:10.1111/biom.13374> Evaluating and improving a matched comparison of antidepressants and bone density. Biometrics, 77(4), 1276-1288.

Yu R. & Rosenbaum, P. R. (2022) <doi:10.1080/10618600.2022.2058001> Graded matching for large observational studies. Journal of Computational and Graphical Statistics, 31(4):1406-1415.

Yu, R. (2023) <doi:10.1111/biom.13771> How well can fine balance work for covariate balancing? Biometrics. 79(3), 2346-2356.

Zhang, B., D. S. Small, K. B. Lasater, M. McHugh, J. H. Silber, and P. R. Rosenbaum (2023) <doi:10.1080/01621459.2021.1981337> Matching one sample according to two criteria in observational studies. Journal of the American Statistical Association, 118, 1140-1151.

Zubizarreta, J.R., 2012. <doi:10.1080/01621459.2012.703874>Using mixed integer programming for matching in an observational study of kidney failure after surgery. Journal of the American Statistical Association, 107(500), pp.1360-1371.

Zubizarreta, J. R., Reinke, C. E., Kelz, R. R., Silber, J. H. and Rosenbaum, P. R. (2011) <doi:10.1198/tas.2011.11072> Matching for several sparse nominal variables in a case control study of readmission following surgery. The American Statistician, 65(4), 229-238.

Zubizarreta, J.R., Stuart, E.A., Small, D.S. and Rosenbaum, P.R. eds. (2023).
<doi:10.1201/9781003102670> Handbook of Matching and Weighting Adjustments for Causal Inference. Boca Raton: FL: Chapman and Hall/CRC Press.
}
\author{
Paul R. Rosenbaum
}
\seealso{
\code{\link{artlessV2}}, \code{\link[iTOS]{makematch}}
}
\note{The mathematical structure of alittleArt() is a very special implementation of the method in Zhang et al. (2023). The method is also described in Chapters 5 and 6 of Rosenbaum (2025).  alittleArt() calls functions in the iTOS package, where more detail may be found.}

\note{\strong{Penalty Structure}:  near.penalty, fine.penalty and integer.penalty relate to matrices
near, fine and xinteger, respectively, and they have a similar structure.
If any of these penalties is a scalar, that scalar is repeated to form a vector with one coordinate for each column of its matrix.  In the example, near has two columns.  For example, if near has two columns, and near.penalty=1000, then near.penalty becomes c(1000,1000).  The penalties apply to the corresponding columns; so, you can apply different penalties to different covariates; however, by default, all columns have the same penalty.}

\note{\strong{The near Matrix}:  An attempt is made to exactly match for covariates in near.  In the example, near contains two binary covariates, namely female and dontSmoke.  This means that the match will try to match women to women and men to men, nonsmokers to nonsmokers, and smokers to smokers.  If near.penalty=c(1000,500), then a mismatch for female increases cost by 1000, a mismatch for dontSmoke increases cost by 500, a mismatch for both costs 1500, and mismatching two people for dontSmoke costs the same as a single mismatch for female.  A small penalty, say near.penalty=c(2,1), will increase he number of exact matches, but will often be overridden by other considerations.}

\note{\strong{The fine Matrix}:  Fine balance refers to the marginal distributions of a covariate in treated and control groups, not to who is paired with whom.  An attempt is made to balance covariates in fine.  In the example, fine includes a covariate expressing four broad age categories, one low education category (less than high school), and a binary covariate distinguishing daily-smokers from everyone else.  This means that the match will work hard to have the same proportion of people with less-than-high-school education in treated and control groups, but it will not prioritize pairing two people with less-than-high-school education.  As with near.penalty above, fine.penalty can be adjusted to increase or decrease the emphasis on fine balancing, or to increase or decrease the emphasis on one column of fine rather than another column.}


\note{\strong{The xinteger Matrix}:  An attempt is made to balance covariates in xinteger, in a manner similar to the covariates in fine.  The difference is that the covariates in fine are viewed as nominal, but the covariates in xinteger are viewed as integers.  Take ageC in the example.  ageC is an ordered factor that is made into an integer using as.integer(ageC). ageC cuts age into 4 categories at 30, 45, and 60.  If used in fine, the categories <30 and 60+ are nominal categories.  If used in xinteger, <30 is far from 60+, and <30 is closer to 30-45 than to 60+.  If ageC cannot be perfectly balanced, the penalty is smaller for imbalances in nearby categories than for distant categories. If integer.penalty=20, then there is 0 penalty for the same category, 20 for a one-category difference, 40 for a two-category difference, etc.  If used as a nominal covariate in fine, every imbalance for ageC would cost the same, ignoring the fact that <30 is closer to 30-45 than to 60+.}

\note{\strong{The Propensity Score}:  Three separate attempts are made to, first, balance the propensity score in the sense of fine balance, and second to pair closely for the propensity score, and third to avoid controls with propensity scores below all treated individuals.  These attempts are controlled by two parameters, min.penalty and pr.penalty.

With the default min.penalty = c(10, 1, 0.05), there is a penalty of 10 for each control that is used in the match whose propensity score is below the minimum propensity score in the treated group.  Also, there is an additional penalty of 1 for each control that is used whose propensity score is below the 0.05 quantile of the propensity scores in the treated group.  Taking
min.penalty = c(0, 0, 0.05) removes this feature.  Taking
min.penalty = c(10, 0, 0.05) uses the minimum penalty but not the 0.05 quantile penalty, etc.  This is a simple directional penalty of the type in Yu and Rosenbaum (2019).  By construction, the propensity score tends to be low in the control group -- so the needed direction is clear -- but the magnitudes of the penalties -- defaulting to 10 and 1 -- may need adjustment based on boxplots of the propensity scores in matched samples.

The propensity score is made into two integer variables, pr6 and pr3, where pr6 is 1 to 6 and cuts pr at its 1/6, 1/3, 1/2, 2/3, and 5/6 quantiles, while pr3 is 1, 2 or 3 and cuts pr at its 1/3 and 2/3 quantiles.  Note that pr3=1 when pr6=1 or pr6=2.  An attempt is made to pair for pr6 and pr3 and to balance for pr6 and pr3.  At the default pr.penalty=c(2,5,25,250), there is a penalty of 2 for a one-category mismatch for p6, and an additional penalty of 5 for a one category mismatch for p3; moreover, as in xinteger, these are doubled for a two category mismatch, etc.  At the default pr.penalty=c(2,5,25,250), there is a 25 penalty for a one category imbalance in p6 and a 250 penalty for a one category imbalance in p3.

Changing pr.penalty to c(0,0,25,250) will make no attempt to pair for the propensity score, while trying to balance it.  Changing pr.penalty to c(2,5,0,0) will make no attempt to finely balance the propensity score, while trying to pair for it.  Changing pr.penalty to c(2000,0,0,0) will try very hard to match exactly for pr6.
}

\note{\strong{The xm Matrix}: The variables in xm are used to construct a robust rank based covariate distance similar to the Mahalanobis distance; see section 9.3 of Rosenbaum (2020a).  Robustness refers to two problems with the usual Mahalanobis distance when used in matching.  First, in the usual Mahalanobis distance, outliers in a covariate can increase its sample variance, thereby decreasing the importance of a 1-unit difference in the covariate.  First, in the usual Mahalanobis distance, a rare binary covariate has a small variance even though a mismatch is always of size |1-0|=1; so, in a US sample, a mismatch for lives-in-Wyoming is much more important than a mismatch for lives-in-California.  The robust rank based covariate distance fixes both problems: an outlier cannot make a covariate less important, and binary variables, rare or common, are equally important.}

\note{\strong{SOLVER}:
The package uses by default the solver rlemon; it is available in R.  The alternative, rrelaxiv, requires a special installation that will now be described.

With solver="rrelaxiv", the package indirectly uses the callrelax() function in Samuel Pimentel's rcbalance package.  This function was originally intended to call the excellent RELAXIV Fortan code of Bertsekas and Tseng (1988,1994).  Unfortunately, that code has an academic license and is not available from CRAN; so, by default the package calls the rlemon function instead, which is available at CRAN.  If you qualify as an academic, then you may be able to download the RELAXIV code from Github at <https://github.com/josherrickson/rrelaxiv/> and use it in artlessV2 by setting solver="rrelaxiv".
}

\note{
-- The following are some \strong{practical tips} on how to use alittleArt().

-- Perhaps do a first match with the default settings.  Examine the balance table and parallel boxplots of covariates in matched treated and control groups.  Adjust the various penalties, if needed, to fix any covariate imbalances you find.

-- It is harder to match 1 treated individual to 3 controls (with ncontrols = 3) than to match one control to each treated individual.  If you are having difficulty finding a well-balanced 1-to-5 match, try 1-to-3 or 1-to-1.  If it was easy to find a balanced 1-to-1 match, try 1-to-3.  With a better choice of penalties, it may be possible to match 1-to-3, while a worse choice of penalties produces adequate covariate balance only for 1-to-1.

-- Most covariates that you want to balance should be included in the propensiy score, either in pr or in x.

-- The covariates in x could include, say: (i) a quadratic in age,
(age-mean(age))^2, (ii) an interaction, (age-mean(age))*(bmi-mean(bmi)), or
(iii) spline terms computed from age. Alternatively, you can build your own propensity score in pr or substitute a different kind of score, rather than automatically using a linear logit model fitted by maximum likelihood.

-- It is sometimes important to look for effect-modification, meaning that the treatment effect varies systematically with one or more covariates.  If you want to look for effect modification by female-vs-male, then it is useful to include a binary female covariate in the matrix near with a large penalty.  This will ensure that all or almost all pairs will be exactly matched for female; so, the pairs can be split into female pairs and male pairs for separate or comparative analysis.  Various re-pairing methods can be used with finely balanced covariates; see Lee et al. (2018).

-- The matrix xbalance allows you to control the ways covariates are expressed in the balance table separately from the way covariates are expressed in the match.

-- An attempt is made to pair closely for covariates in xm.  A continuous covariate, like age or bmi, might be placed in x and in xm.  A binary covariate like female can also be used.  Covariates in xm are given roughly equal importance; so, do not put unimportant covariates in xm.

-- The match should be finalized before any outcome information is examined; see Rubin (2008).

-- There can exist treated and control groups that cannot be matched.  If all of the treated individuals are under age 20 and all of the controls are over age 50, then there is no way you can match for age.  You could do regression or covariance adjustment for age, but of course it would be silly.  Matching will often stop you from doing silly things, while regression will let you do silly things.
}

\note{\strong{TECHNICAL DETAILS}:  The following details refer to Figure 1 in Zhang et al. (2023) or Figure 5.5 in Rosenbaum (2025).  In particular LEFT refers to treated-control edges on the left side of the network, RIGHT refers to control-treated edges on the right side of the network, and CC refers to the control-control edges in the center of the network.  Various functions from the iTOS package are mentioned; they have detailed documentation in the iTOS package.

xm adds a penalty on the LEFT using addMahal() from the iTOS package.

near adds a penalty on the LEFT using addNearExact() from the iTOS package.

fine adds a penalty on the RIGHT using addNearExact() from the iTOS package.

xinteger adds a penalty on the RIGHT using addinteger() from the iTOS package.

pr does three things:

1.  the first two penalties in pr.penalty are added on the LEFT using addinteger() from the iTOS package.

2.  the last two penalties in pr.penalty are added on the RIGHT using addinteger() from the iTOS package.

3.  min.penalty adds two penalties to CC central edges.
}

\examples{
\donttest{
# The example below uses the binge data from the iTOS package.
# See the documentation for binge in the iTOS package for more information.
#
library(iTOS)
data(binge)
b2<-binge[binge$AlcGroup!="P",] # Match binge drinkers to nondrinkers
z<-1*(b2$AlcGroup=="B") # Treatment/control indicator
b2<-cbind(b2,z)
rm(z)
rownames(b2)<-b2$SEQN
attach(b2)
# Estimate a propensity score
pr<-stats::glm(z~age+female+education+bmi+vigor+
      smokenow+smokeQuit+bpRX,family=binomial)$fitted.values
#
#  Create nominal covariates to include in near or fine
#
smoke<-1*(smokenow==1)
dontSmoke<-1*(smokenow==3)
age50<-1*(age>=50)
bmi30<-1*(bmi>=30)
ed2<-1*(education<=2)
smoke<-1*(smokenow==1)
#
#  near contains covariates to be matched as exactly as possible
#
near<-cbind(female,dontSmoke)
#
# xm contains covariates in the robust Mahalanobis distance
# Includes some continuous covariates.
#
xm<-cbind(age,bmi,vigor,smokenow,education)
#
# fine contains covariate that will be balanced, but not matched
#
fine<-cbind(ed2,smoke,dontSmoke)

# variable to be used in xinteger
ageCi<-as.integer(ageC)
xbalance<-cbind(pr,age,female,education,bmi,vigor,smokenow,smokeQuit,bpRX,
   ageCi,ed2,smoke,dontSmoke,bmi30,smoke,ed2,age50)
b2<-cbind(b2,pr)
rm(bmi30,smoke,ed2,age50,dontSmoke)
detach(b2)

mc<-alittleArt(b2,b2$z,pr=pr,xm=xm,near=near,fine=fine,xinteger=ageCi,
   ncontrols=3,xbalance=xbalance,pr.penalty = c(3, 5, 50, 250))
#
#  Here are the first two 1-to-3 matched sets.
#
mc$match[1:8,]
#
#  You can check that every matched set is exactly matched for
#  female and nonsmoking.  This is from near-exact matching.
#  In some other data set, the number of mismatches might be
#  minimized, not driven to zero.
#
#  The balance table shows that large imbalances in covariates
#  existed before matching, but are much smaller after matching.
#  Look, for example, at the propensity score, female, and
#  the several versions of the smoking variable.
#
mc$balance
m<-mc$match
m<-m[m$matched,] # Remove the unmatched controls
table(m$z)
prop.table(table(m$ageC,m$z),2)
# You could improve this table by setting integer.penalty=500.
# Other things might suffer a bit.  The boxplot of age is good as is.
boxplot(m$age~m$z)
boxplot(m$pr~m$z)
}
}

