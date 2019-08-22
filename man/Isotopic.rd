\name{Isotopic}
\alias{Isotopic}

\title{
Preparation library Full-MS.
}
\description{
Generates a given a library of MS1 metabolites, a list of adducts and isotopes. Adapted from R package OrgMassSpecR.
}
\usage{
Isotopic(library, name=c("M+H"),adducts)
}

\arguments{
  \item{library}{
   Library of MS1 metabolites.
}
  \item{name}{
  Name of the adduct to generate.
  }
  \item{adducts}{
  List of adducts.
  }}

\author{Manuel D Peris Diaz}
\references{
1. R-MetaboList: a flexible tool for metabolite extraction from high-resolution data-independent acquisition mass spectrometry analysis. Metabolites. Soon

2. A Survey of Orbitrap All Ion Fragmentation Analysis Assessed by an R MetaboList Package to Study Small-Molecule Metabolites. Chromatographia. 2018, 81, 981-994.

}
\examples{
#library(MetaboList)
#Reading the library and indicating type of adduct
#Isotopic(library,name=c("M+H"))



}
