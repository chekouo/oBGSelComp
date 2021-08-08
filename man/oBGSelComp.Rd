\name{oBGSelCompAlgo}
\alias{oBGSelCompAlgo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	An MCMC algorithm to perform an overlapping Bayesian group selection for compostional response data.
}
\description{
The function implements an MCMC algorithm for a Bayesian model for compositional data using a Dirichlet  model by allowing for straightforward incorporation of available high-dimensional covariate information within a log-linear regression framework. The method imposes prior distributions that account for the overlapping structure between groups.
}

\usage{
oBGSelCompAlgo(outcome=outcome,X=X,Pathway=Pathway,sample=10000, 
burnin=1000,hyperprob=c(1,1),hypersig=c(5,0.01),hyperlamb=c(5,2), seed=1)
}

\arguments{
\item{outcome}{ A compositional multivariate response of dimension \eqn{n\times q}}
\item{X}{ A set of covariate matrix of dimension \eqn{n\times p} (e.g. gene expression)}
\item{Pathway}{ A binary pathway indicator matrix \eqn{P} of dimension \eqn{K\times p} where cell \eqn{(k,j)} is 1 if the  covariate \eqn{j}  belongs to group (or pathway) \eqn{k}, and \eqn{0} otherwise.   }
\item{sample}{ Total number of MCMC draws. It must be larger than burnin. }
\item{burnin}{ Number of draws to discard for burn-in}
	\item{hyperprob}{ Hyperparameters of a beta(a,b) distribution, prior distribution of the probability of selecting groups}
		\item{hypersig}{ Hyperparameters of an inverse gamma  inv.gamma(a,b) distribution, prior distribution of the variance\eqn{\sigma^2_l}}
	\item{hyperlamb}{ Hyperparameters of a gamma(a,b) distribution, prior distribution of the shrinkage parameter \eqn{\lambda_l}}
	\item{seed }{Random number are generated with seed=1 during the MCMC algorithm}
}
\details{
	The function will return two R objects, which can be assigned to a variable. To see the results, use the ``$" operator. 
}
\value{
\item{PostProbGrp}{ Marginal posterior probabilities of the group selection indicators, \eqn{Prob(\gamma_{lk}=1)}}
\item{BetaSample}{ An MCMC sample of the regression effects \eqn{\beta_{lk}} }
}

\references{
 Thierry Chekouo et al (2021), \emph{A Bayesian group selection with compositional responses for analysis of radiologic tumor proportions and their genomic determinants}, \emph{submitted}.
}

\seealso{
\code{\link{SimulateData}}
}

\examples{
library(oBGSelComp)
## overlap groups (10 percent overlap)
dat=SimulateData(n=100,betaAbs=1,sig=.1,r2=0.5,propOverl=10,seed=1)
Res=oBGSelCompAlgo(outcome=dat$y,X=dat$X,Pathway=dat$Pathway,sample=10000, 
burnin=1000,hyperprob=c(1,1),hypersig=c(4,0.01),hyperlamb=c(5,2), seed=1)
str(Res)
library(AUC)
q=4
AUC1=rep(0,q)
for (l in 1:q){
AUC1[l]= auc(roc(Res$PostProbGrp[l,],as.factor(dat$GamKnown[l,])))
}
AUC1
}

