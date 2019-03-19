\name{CE.isolation}
\alias{CE.isolation}

\title{
Separation of MS/MS files with regards collision energy
}
\description{
Isolation of MS/MS events acquired at different collision energies into single files at particular collision energy
}
\usage{
CE.isolation(file, output)
}

\arguments{
  \item{file}{
   A mzXML file from the LC-MS/MS experiment in positive or negative ionization mode.
}
  \item{output}{
Multiple mzXML files separated by collision energies.
  }}

\author{Manuel D Peris Diaz}
\references{
1. R-MetaboList: a flexible tool for metabolite extraction from high-resolution data-independent acquisition mass spectrometry analysis. Metabolites. Soon

2. A Survey of Orbitrap All Ion Fragmentation Analysis Assessed by an R MetaboList Package to Study Small-Molecule Metabolites. Chromatographia. 2018, 81, 981-994.

}
\examples{
library(MetaboList)
#Reading the file.mzXML
#CE.isolation("AIFpos1000-AIF.mzXML","fileposB")


}
