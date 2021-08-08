\docType{data}
\name{HallmarkPathwayData}
\alias{HallmarkPathwayData}
\title{Data sets of patients with Glioblastoma and Hallmark pathway information}
\format{A list of 4 data frames:
\describe{
  \item{TumorProportion}{Brain tumor proportions: edema, necrotic, non-enhancing and enhancing}
  \item{mRNA}{mRNA expression of 1172 genes.}
  \item{HallmarkPathway}{Hallmark pathway information: it's a binary matrix that indicates the membership of a gene within a Hallmark pathway}
  \item{Hallmarkpathnames}{Hallmark pathway names}
}}

\usage{
HallmarkPathwayData
}
\description{
A list of 4 data frames that contains tumor proportion, mRNA expression, and Hallmark pathway information and thier names.
}
\examples{
data("HallmarkPathwayData")

}
\references{
Thierry Chekouo et al (2021), \emph{A Bayesian group selection with compositional responses for analysis of radiologic tumor proportions and their genomic determinants}, \emph{submitted}.
}
\keyword{HallmarkPathwayData}
