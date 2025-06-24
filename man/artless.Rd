\name{artless}
\alias{artless}
\concept{Causal inference}
\concept{Observational Studies}
\concept{Propensity score}
\concept{Matching}
\concept{Fine balance}
\concept{Near exact matching}
\title{
Artless Automatic Matching
}
\description{
Implements a simple version of multivariate matching using a propensity score, near-exact matching, near-fine balance, and robust Mahalanobis distance matching.  You specify the variables, and the program does everything else.  Should you
be artful, not artless?  See the notes.
}
\usage{
artless(dat, z, x, xm = NULL, near = NULL, fine = NULL,
   ncontrols = 1, rnd = 2, solver="rlemon")
}
\arguments{
  \item{dat}{
A dataframe containing the data set that will be matched.  Let N be the
number of rows of dat.
}
  \item{z}{
A binary vector of length N where z[i]=1 if the ith row of dat describes
a treated individual and z[i]=0 if the ith row of dat describes a control.
}
  \item{x}{
x is a numeric matrix with N rows.  The covariates in x are used to estimate
a propensity score using a linear logit model.
}
  \item{xm}{
xm is a numeric matrix with N rows.  The covariates in xm are used to define
a robust Mahalanobis distance between treated and control individuals.
}
  \item{near}{
A numeric vector of length N or a numeric matrix with N rows.  Each column of near should represent levels of a nominal covariate with two or a few levels.  The variables in near are used in near-exact matching.
}
  \item{fine}{
A numeric vector of length N or a numeric matrix with N rows.  Each column of fine should represent levels of a nominal covariate with two or a few levels.  The variables in fine are used in near-fine balancing.
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
}
\details{
This package builds a matched treated-control sample from an unmatched data set.  It asks you to designate roles for specific covariates, and it does the rest.  It is described as ``artless automatic matching'' because it makes decisions by default.  Perhaps you could make better decisions; if so, perhaps try the iTOS package which gives you much more control over
decisions.  The package will often create a reasonable matched sample with little effort; however, it also could be used as a first step in learning the art of constructing a matched sample.  Wittgenstein spoke of a the ``ladder you throw away after you have climbed it,'' and the package can also serve that function.
}
\value{
\item{match }{A dataframe containing the matched data set.  match contains the rows of dat in a different order.  match adds two columns to dat, called mset and matched, which identify matched pairs or matched sets.  Specifically, matched is TRUE if a row is in the matched sample and is FALSE otherwise.  Rows of dat that are in the same matched set have the same value of mset.  The rows of match are sorted by mset with the treated individual before the matched controls.  The unmatched controls with matched=FALSE appear as the last rows of match.  match also adds the estimated propensity score as a probability pr.
When you analyze the matched data, you will want to remove rows of match with matched==FALSE.}
\item{balance }{A matrix called the balance table.  The matrix has one row for each covariate in x, xm, near and fine; so, some covariates may be repeated.  It also has a first row for the propensity score.  There are five columns.
Column 1 is the mean of the covariate in the treated group.  Column 2 is the mean of the covariate in the matched control group.  Column 3 is the mean of the covariate among all controls prior to matching.  Column 4 is the difference between columns 1 and 2 divided by a pooled estimate of the standard deviation of the covariate before matching.  Column 5 is the difference between columns 1 and 3 divided by a pooled estimate of the standard deviation of the covariate before matching.  Notice that columns 4 and 5 have the same denominator, but different numerators.}
}
\references{
Bertsekas, D. P., Tseng, P. (1988) <doi:10.1007/BF02288322> The Relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.

Bertsekas, D. P. (1990) <doi:10.1287/inte.20.4.133> The auction algorithm for assignment and other network flow problems: A tutorial. Interfaces, 20(4), 133-149.

Bertsekas, D. P., Tseng, P. (1994) <http://web.mit.edu/dimitrib/www/Bertsekas_Tseng_RELAX4_!994.pdf> RELAX-IV: A Faster Version of the RELAX Code for Solving Minimum Cost Flow Problems.

Greifer, N. and Stuart, E.A., (2021). <doi:10.1093/epirev/mxab003> Matching methods for confounder adjustment: an addition to the epidemiologist’s toolbox. Epidemiologic Reviews, 43(1), pp.118-129.

Hansen, B. B. and Klopfer, S. O. (2006) <doi:10.1198/106186006X137047> "Optimal full matching and related designs via network flows". Journal of computational and Graphical Statistics, 15(3), 609-627. ('optmatch' package)

Hansen, B. B. (2007) <https://www.r-project.org/conferences/useR-2007/program/presentations/hansen.pdf> Flexible, optimal matching for observational studies. R News, 7, 18-24. ('optmatch' package)

Pimentel, S. D., Yoon, F., & Keele, L. (2015) <doi:10.1002/sim.6593> Variable‐ratio matching with fine balance in a study of the Peer Health Exchange. Statistics in Medicine, 34(30), 4070-4082.

Niknam, B.A. and Zubizarreta, J.R. (2022). <10.1001/jama.2021.20555> Using cardinality matching to design balanced and representative samples for observational studies. JAMA, 327(2), pp.173-174.

Pimentel, S. D., Kelz, R. R., Silber, J. H. and Rosenbaum, P. R. (2015)
<doi:10.1080/01621459.2014.997879> Large, sparse optimal matching with refined covariate balance in an observational study of the health outcomes produced by new surgeons. Journal of the American Statistical Association, 110, 515-527.

Rosenbaum, P. R. and Rubin, D. B. (1985) <doi:10.1080/00031305.1985.10479383> Constructing a control group using multivariate matched sampling methods that incorporate the propensity score. The American Statistician, 39, 33-38.

Rosenbaum, P. R. (1989) <doi:10.1080/01621459.1989.10478868> Optimal matching for observational studies. Journal of the American Statistical Association, 84(408), 1024-1032.

Rosenbaum, P. R., Ross, R. N. and Silber, J. H. (2007) <doi:10.1198/016214506000001059> Minimum distance matched sampling with fine balance in an observational study of treatment for ovarian cancer. Journal of the American Statistical Association, 102, 75-83.

Rosenbaum, P. R. (2020a) <doi:10.1007/978-3-030-46405-9> Design of Observational Studies (2nd Edition). New York: Springer.

Rosenbaum, P. R. (2020b). <doi:10.1146/annurev-statistics-031219-041058> Modern algorithms for matching in observational studies. Annual Review of Statistics and Its Application, 7(1), 143-176.

Rosenbaum, P. R. and Zubizarreta, J. R. (2023). <doi:10.1201/9781003102670>
Optimization Techniques in Multivariate Matching. Handbook of Matching and Weighting Adjustments for Causal Inference, pp.63-86.  Boca Raton: FL: Chapman and Hall/CRC Press.

Rosenbaum, P. R. (2025) Introduction to the Theory of Observational Studies.
New York: Springer.

Rubin, D. B. (1980) <doi:10.2307/2529981> Bias reduction using Mahalanobis-metric matching. Biometrics, 36, 293-298.

Stuart, E.A., (2010). <doi:10.1214/09-STS313> Matching methods for causal inference: A review and a look forward. Statistical Science, 25(1), 1-21.

Yang, D., Small, D. S., Silber, J. H. and Rosenbaum, P. R. (2012) <doi:10.1111/j.1541-0420.2011.01691.x> Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes. Biometrics, 68, 628-636.

Yu, Ruoqi, and P. R. Rosenbaum. <doi:10.1111/biom.13098> Directional penalties for optimal matching in observational studies. Biometrics 75, no. 4 (2019): 1380-1390.

Yu, R., Silber, J. H., & Rosenbaum, P. R. (2020) <doi:10.1214/19-STS699> Matching methods for observational studies derived from large administrative databases. Statistical Science, 35(3), 338-355.

Yu, R. (2021) <doi:10.1111/biom.13374> Evaluating and improving a matched comparison of antidepressants and bone density. Biometrics, 77(4), 1276-1288.

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

\note{
-- The following are some practical tips on how to use artless.

-- Placing a covariate in x means that it is included in the propensity
score.  Most or all covariates that you want to balance should be placed
in x.

-- A limited number of nominal covariates with a few levels can be placed in near or in fine.  Both near and fine covariates are given overriding importance; so, if you place too many covariates in near or fine, or if they have too many levels, they will override everything else, and the match quality will be poor.  The same covariate can appear, perhaps in different forms, in x, xm, near and fine.  In the example, a five-level education variable is in x and xm, and a two-level education variable formed from the five-level education variable is in fine.

-- An attempt is made to exactly match for covariates in near.  In the example, near contains two binary covariates, namely female and dontSmoke.  This means that the match will try whenever possible to match women to women and men to men, nonsmokers to nonsmokers, and smokers to smokers.  Other considerations are subbordinated to this goal.

-- An attempt is made to balance covariates in fine.  In the example, fine includes a covariate expressing four broad age categories, one low education category (less than high school), and a binary covariate distinguishing daily-smokers from everyone else.  This means that the match will work hard to have the same proportion of people with less-than-high-school education in treated and control groups, but it will not prioritize pairing two people with less-than-high-school education.  Although subbordinate to near exact matching, fine balance is given more importance than other considerations.

-- Two separate attempts are made to, first, balance the propensity score in the sense of fine balance and to pair closely for the propensity score.  More emphasis is given to balancing the propensity score, much less to pairing for it.  The match also tries in a limited way to avoid using many controls whose propensity scores are below the minimum propensity score in the treated group.

-- An attempt is made to pair closely for covariates in xm; however, this task has the lowest priority of the several goals.  A continuous covariate, like age or bmi, might be placed in x and in xm.  Covariates in xm are given roughly equal importance, so do not put unimportant covariates in xm.

-- The covariates in x could include, say: (i) a quadratic in age,
(age-mean(age))^2, (ii) an interaction, (age-mean(age))*(bmi-mean(bmi)), or
(iii) spline terms computed from age. The propensity score is fitted as a linear logit model in the covariates in x, but you can fit various nonlinear propensity scores by passing in x various nonlinear transformations of a more limited set of covariates.

-- Usually, the first match you construct is imperfect, and you see this in the balance table or in plots of the matched data.  So, you make small adjustments to x, xm, near and fine to fix the imperfections.
The match should be finalized before any outcome information is examined.  Taking the first match without looking at it and improving it is not artless; it is incompetent.

-- Once you have developed some experience with the artless function, you may want to learn about other artful tactics that can enhance your ability to remove imperfections in a match.  Some of these tactics are implemented in the iTOS package that is called by artless.

-- There are treated and control groups that cannot be matched.  If all of the treated individuals are under age 20 and all of the controls are over age 50, then there is no way you can match for age.  You could do regression or covariance adjustment for age, but of course it would be silly.  Matching will often stop you from doing silly things, while regression will let you do silly things.
}

\note{
Should you be artful rather than artless?  Essentially, the artless() function
is setting priorities by default.  This makes artless() easy to use, but its default priorities might not be your priorities.  An alternative is to set your own priorities by using the matching methods in, say, the iTOS package.  The artless() function calls the functions in the iTOS package, but it sets default priorities when it does this.  There are also many more options in the iTOS package.

What can artful use of iTOS do that artless() cannot?  artless() automatically sets priorities and penalties, but iTOS lets you adjust them.  artless() automatically gives an emphasis to the propensity score, and does this in a particular way, but iTOS lets you decide.  The directional penalties of Yu and Rosenbaum (2019) need to be titrated to produce desired effects; they are in iTOS but not in artless().  Near-exact and near-fine matching are implemented for nominal variables in artless(), but iTOS has other options for ordered categories.  iTOS lets you give more emphasis to one covariate, less to another, but artless() does this only indirectly through the matrices x, xm, near and fine.  In artless() all variables in near are treated as equally important, and all variables in fine are treated as equally important, but iTOS lets you decide.  Caliper matching is possible in iTOS but not in artless().  artless() uses the control-control edge costs in Zhang et al. (2023) to avoid low propensity scores in the control group, but iTOS lets you use this feature any way you prefer.  The iTOS package is associated with Rosenbaum (2025), especially its Chapters 5 and 6.
}

\note{
This note provides some references and detail about what the package is actually doing.  You do not have to read this note to use the package.

Matching using propensity scores and a Mahalanobis distance is discussed in Rosenbaum and Rubin (1985).  The robust Mahalanobis distance is discussed in Section 9.3 of Rosenbaum (2020a) and more briefly in Section 4.1 of Rosenbaum (2020b).

Near-exact matching (also known as almost-exact matching) is an attempt to match exactly for a few nominal covariates, while also matching for other things.  It is described in Sections 10.3 and 10.4 of Rosenbaum (2020a) and more briefly in Section 4.3 of Rosenbaum (2020b).  Near-exact matching is implemented by a large penalty added to a covariate distance: if two people are not exactly matched for a near-exact covariate, then the covariate distance between them is very large.  Near-exact matching minimizes the number of individuals who are not exactly matched.

Fine balance attempts to balance a covariate without pairing for it.  For example, female is balanced if the treated and control groups have the same proportion of females, but female is exactly matched if females are always matched to females.  Fine balance is discussed in Chapter 11 of Rosenbaum (2020a) and more briefly in Section 4.4 of Rosenbaum (2020b).  Fine balance was introduced in Section 3.2 of Rosenbaum (1989), and is further developed in Rosenbaum, Ross and Silber (2007).  If one seeks a match as close as possible to fine balance, then one is doing near-fine balance.  Near-fine balance is often implemented using penalties for imbalances; see Yang et al. (2012), Pimentel et al. (2015) and Zhang et al. (2023).

One can do near-exact matching and fine balancing of the same variable, perhaps leading the proportion of females to be exactly the same in treated and control groups, with pairs matched for female as often as is possible.  See Zubizarreta et al. (2011) for discussion.

artless() uses the control-control edge costs in Zhang et al. (2013) to moderately penalize the use of a control whose propensity score is below the minimum propensity score in the treated group.  This penalty is smaller than the penalty for near-exact matching and for aspects of propensity score balancing, but it is larger than the penalty for each variable in near-fine matching.

This package implements a very specific version of two-criteria matching from Zhang et al. (2023) using functions from the iTOS package.  Two-criteria matching integrates a number of earlier techniques into a single network structure.  The package picks several one-size-fits-all penalties for distances for two-criteria matching.  An artful match might vary penalties in a thoughtful way to achieve a better, closer, more balanced match with a larger value of ncontrols.  The package does not use asymmetric calipers and directional penalties from Yu and Rosenbaum (2019) because these are not easily automated, but the artful use of these techniques can produce a better match.

The package uses optimal matching by minimum cost flow in a network.  See Bertsekas (1990) for an introduction to this optimization technique, and see Rosenbaum (1989) for its application to matching in observational studies.

The package indirectly uses the callrelax() function in Samuel Pimentel's rcbalance package.  This function was originally intended to call the excellent RELAXIV Fortan code of Bertsekas and Tseng (1988,1994).  Unfortunately, that code has an academic license and is not available from CRAN; so, by default it calls the rlemon function instead, which is available at CRAN.  If you qualify as an academic, then you may be able to download the RELAXIV code from Github at <https://github.com/josherrickson/rrelaxiv/> and use it in artless by setting solver="rrelaxiv".

artless() uses a dense network, so it can match moderately large data sets, but not very large data sets.  For very large data sets, see Yu et al. (2020) and Yu's bigmatch package in R.

Network optimization is only one of several optimization techniques that may be used in multivariate matching.  See Niknam and Zubizarreta (2022),
Zubizarreta (2012) and Rosenbaum and Zubizarreta (2023).
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
#
agec<-as.integer(ageC)
#
# x contains the variables in the propensity score
#
x<-data.frame(age,female,education,bmi,vigor,smokenow,smokeQuit,bpRX)
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
fine<-cbind(ageC,ed2,smoke,dontSmoke)
rm(agec,bmi30,smoke,ed2,age50)
detach(b2)

mc<-artless(b2,b2$z,x,xm=xm,near=near,fine=fine,ncontrols=3)
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
}
}

