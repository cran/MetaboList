\name{PeakGroup}
\alias{PeakGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Peak Grouping for multiple DIA files
}
\description{
Annotated metabolites from single DIA files acquired at different collision energies are subjected to peak grouping. The function groups metabolites that are presented along the DIA files.
}
\usage{
PeakGroup(aif1,aif2,aif3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aif1}{
   Result obtained from \link[MetaboList]{AIF}
}
  \item{aif2}{
Result obtained from \link[MetaboList]{AIF}
  }
  \item{aif3}{
Result obtained from \link[MetaboList]{AIF}
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
#Reading the file.mzXML for Full-MS scan and MS/MS.
#The files were previously separated with the SEPC.R function.

#fileMS<-"fullMS.mzXML"
#fileMS2CE5<-"fileMS2CE5.mzXML"

# Separation by collision energy
#CE.isolation("AIFpos1000-AIF.mzXML","fileposB")

#Reading the database.csv file:
# database<- read.csv("C:/database.csv")

#Processing peak-picking and annotation with default parameters

#aif5<-AIF(fileMS,fileMS2CE5,databasepos,CE=5,
#ion_mode = "positive",mzw = 0.005,rtw = 7
#,minintMS2=1,minintMS1=1)
#aif10<-AIF(fileMS,fileMS2CE10,databasepos,CE=10,
#ion_mode = "positive",rtw = 7, mzw = 0.005
#,minintMS2=1,minintMS1=1)
#aif20<-AIF(fileMS,fileMS2CE20,databasepos, CE=20,
#ion_mode = "positive",rtw = 7, mzw = 0.005,
#minintMS2=1,minintMS1=1)

# Peak grouping

#Peakgroup<-PeakGroup(aif5,aif10,aif20)


}
