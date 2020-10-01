\name{SimulateData}
\alias{SimulateData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation of simulated data as explained in the reference manuscript.	      
}
\description{ 
	This function generates compositional response and covariate data as  described in the manuscript (see Section 5.2 in the manuscript).
}

\usage{
SimulateData(n=100,betaAbs=1,sig=.1,r2=0.5,propOverl=0,seed=1)
}

\arguments{
\item{n}{ The sumber of subjects/samples}
\item{betaAbs}{ Effect values in absolute values of impaortant features}
\item{sig}{ Value of  the variance of latent values alpha's given the feature set}
\item{r2}{ Variance value that represents correlation of features within pathways
}
\item{propOverl}{ Proportion of overlap between "consecutive" groups. It can only takes three values;0,10 and 50}
\item{seed}{ Seed to generate random numbers}
}

\details{
The function will return three data matrices y, X, Pathway and GamKnown which are respectively the compositional response data, the covariate matrix, the Pathway matrix indicator and the true gamma (gamma(l,k)=1 if pathway k is important for outcome l, and 0 otherwise) 
}
\references{
	   Thierry Chekouo et al (2020), \emph{A Bayesian Group Selection approach for the analysis of volumetric images of brain cancers and their genomic determinants}, \emph{submitted}.
}

\examples{
dat=SimulateData(n=100,betaAbs=1,sig=.1,r2=0.5,propOverl=0,seed=1)
str(dat)
}
