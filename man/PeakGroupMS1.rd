\name{PeakGroupMS1}
\alias{PeakGroupMS1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Peak Grouping for MS1 file
}
\description{
Annotated metabolites from single MS1 files are subjected to peak grouping.
}
\usage{
PeakGroupMS1(annotation_pos)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotation_pos}{
   Result obtained from \link[MetaboList]{FullMS}
}}



\value{
  \item{csv}{A Peakgroup.csv files for each metabolite}
}


\author{Manuel D Peris Diaz}
\references{
1. R-MetaboList: a flexible tool for metabolite extraction from high-resolution data-independent acquisition mass spectrometry analysis. Metabolites. Soon

2. A Survey of Orbitrap All Ion Fragmentation Analysis Assessed by an R MetaboList Package to Study Small-Molecule Metabolites. Chromatographia. 2018, 81, 981-994.


}
\examples{
library(MetaboList)

#fullmsposH<-FullMS(file, databaseH,mzw = 0.0015,rtw=NULL)



#PeakGroupMS1_Qtof<-PeakGroupMS1(fullmsposH)




}
