\docType{data}
\name{GBMPathwayData}
\alias{GBMPathwayData}
\title{Data sets of patients with Glioblastoma and GBM pathway information}
\format{A list of 4 data frames:
\describe{
  \item{TumorProportion}{Brain tumor proportions: edema, necrotic, non-enhancing and enhancing}
  \item{mRNA}{mRNA expression of 718 genes.}
  \item{GBMPathway}{GBM pathway information: it's a binary matrix that indicates the membership of a gene within an GBM pathway}
  \item{GBMpathwaynames}{GBM pathway names}
}}

\usage{
GBMPathwayData
}
\description{
A list of 4 data frames that contains tumor proportion, mRNA expression, and GBM pathway information and their names.
}
\examples{
data("GBMPathwayData")

}
\references{
Thierry Chekouo et al (2021), \emph{A Bayesian group selection with compositional responses for analysis of radiologic tumor proportions and their genomic determinants}, \emph{submitted}.
}
\keyword{datasets}
