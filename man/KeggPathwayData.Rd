\docType{data}
\name{KeggPathwayData}
\alias{KeggPathwayData}
\title{Data sets of patients with Glioblastoma and Kegg pathway information}
\format{
A list of 4 data frames:
\describe{
  \item{TumorProportion}{ Brain tumor proportions: edema, necrotic, non-enhancing and enhancing}
  \item{mRNA}{ mRNA expression of 364 genes.}
  \item{KeggPathway}{ KEGG pathway information: it's a binary matrix that indicates the membership of a gene within a KEGG pathway.}
  \item{Keggpathwaynames}{ KEGG pathway names.}
}
}

\usage{
KeggPathwayData
}
\description{
A list of 4 data frames that contains tumor proportion, mRNA expression, and KEGG pathway information and their names.
}
\examples{
data("KeggPathwayData")

}
\references{
Thierry Chekouo et al (2021), \emph{A Bayesian group selection with compositional responses for analysis of radiologic tumor proportions and their genomic determinants}, \emph{submitted}.
}
\keyword{datasets}
