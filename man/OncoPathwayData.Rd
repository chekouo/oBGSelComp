\docType{data}
\name{OncoPathwayData}
\alias{OncoPathwayData}
\title{Data sets of patients with Glioblastoma and Oncogenetic pathway information}
\format{A list of 4 data frames:
\describe{
  \item{TumorProportion}{Brain tumor proportions: edema, necrotic, non-enhancing and enhancing}
  \item{mRNA}{mRNA expression of 1032 genes.}
  \item{OncoPathway}{Oncogenetic pathway information: it's a binary matrix that indicates the membership of a gene within an Oncogenetic pathway}
  \item{Oncopathnames}{Oncogenetic signature names}
}}

\usage{
OncoPathwayData
}
\description{
A list of 4 data frames that contains tumor proportion, mRNA expression, and Oncogenetic pathway information and their names.
}
\examples{
data("OncoPathwayData")

}
\references{
Thierry Chekouo et al (2021), \emph{A Bayesian group selection with compositional responses for analysis of radiologic tumor proportions and their genomic determinants}, \emph{submitted}.
}
\keyword{datasets}
