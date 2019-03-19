\name{FullMS}
\alias{FullMS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Processing and Annotation of LC-MS Full-Scan data
}
\description{
Peak picking of MS data is performed by the enviPick algorithm embedded. Second, it is performed a targeted extraction with a mass tolerance and m/z interval windows constraints for general peak grouping and library interrogation. Retention time might be considered as optional constraints. Library listing needs to follow the format following the exampled attached.
}
\usage{
FullMS(file, database, rtw = 10, mzw = 0.001,
dmzgap = 50,dmzdens = 20, drtgap = 25, drtsmall = 50,
drtdens = 20, drtfill = 5,drttotal = 100, minpeak = 5,
recurs = 3, weight = 1,SB = 2, SN = 1.5, minint = 1000,
maxint = 9e+09, ion_mode = "positive",
ppm = TRUE,ended=6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
   A mzXML file from the LC-MS/MS experiment in positive or negative ionization mode.
}
  \item{database}{
  A file with data arranged in columns as follows: Molecular Formula; Retention time (optional); Neutral mass; Compound name
  }
  \item{rtw}{
  numeric. The difference between the theoretical retention time value and the experimental. Default value=3
  }
  \item{mzw}{
  numeric. The difference between the theoretical m/z value and the experimental (Da). Default value=0.004
  }
  \item{dmzgap}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{dmzdens}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{drtgap}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{drtsmall}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{drtdens}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{drtfill}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{drttotal}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{minpeak}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{recurs}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{weight}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{SB}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{SN}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{minint}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{maxint}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{ion_mode}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{ppm}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{ended}{Arguments to be passed from \link[enviPick]{enviPickwrap}}

}



\value{
  \item{ms}{Annotated metabolites.}
  \item{RawData1}{Raw scans from raw data.}
  \item{PP}{Results obtained throughout enviPick algorithm performed on MS level 1.}
  \item{Peaklist}{Matrix with the peak picked obtained throughout enviPick algorithm.}
}


\author{Manuel D Peris Diaz}
\references{
1. R-MetaboList: a flexible tool for metabolite extraction from high-resolution data-independent acquisition mass spectrometry analysis. Metabolites. Soon

2. A Survey of Orbitrap All Ion Fragmentation Analysis Assessed by an R MetaboList Package to Study Small-Molecule Metabolites. Chromatographia. 2018, 81, 981-994.
}
\examples{
library(MetaboList)
#Reading the file.mzXML
# file<-fullMS.mzXML

#Reading the database.csv file:
# database<- read.csv("C:/FullMS1.csv")

#Processing peak-picking and annotation with default parameters
# FullMS_results<-FullMS(file,database, ion_mode = "positive",)

#Output:
#FullMS_results$ms

}
