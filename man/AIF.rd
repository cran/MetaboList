\name{AIF}
\alias{AIF}
\title{Automatic Metabolite Annotation from LC-MS DIA experiments.}
\description{Analysis and annotation of LC-MS/MS DIA data with the use of in-house mass spectral libraries.}
\usage{AIF(fileMS1,fileMS2,CE=0, database,rtw=7,
mzw=0.05,dmzgap=50,drtdens =20,drtgap=25,
drtsmallMS1=100,drtsmallMS2=30,dmzdensMS1=15,
dmzdensMS2=30,drtfill=5,drttotal=100,minpeakMS1=5,
minpeakMS2=3,recurs=2,weight=2,SB=3, SN=2,
minintMS1=1000,minintMS2=100,maxint=9e+09,
ion_mode="positive",ppm=TRUE,ended=6)}
\arguments{
  \item{fileMS1}{
   A .mzXML file extension with the MS/MS experiment obtained at particular collision energy (CE).
}
  \item{fileMS2}{
   A .mzXML file extension with the MS/MS experiment obtained at particular collision energy (CE).
}
  \item{CE}{
   Collision energy employed for the MS/MS experiment.
}
  \item{database}{
  A csv file with data arranged in columns including the names in the first row: Metabolite; Monoisotopic mass for the precursor; Mass Fragment 1; Mass Fragment 2; Mass Fragment 3... There is no need to used retention times as a constraint.
  }
  \item{rtw}{
  numeric. The difference between the theoretical retention time value and the experimental. Default value=3
  }
  \item{mzw}{
  numeric. The difference between the theoretical m/z value and the experimental. Default value=0.004
  }
  \item{dmzgap}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{dmzdensMS1}{Arguments to be passed from \link[enviPick]{enviPickwrap}for MS1 mode}
    \item{dmzdensMS2}{Arguments to be passed from \link[enviPick]{enviPickwrap}for MS2 mode}

  \item{drtgap}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{drtsmallMS1}{Arguments to be passed from \link[enviPick]{enviPickwrap} for MS1 mode}
    \item{drtsmallMS2}{Arguments to be passed from \link[enviPick]{enviPickwrap} for MS2 mode}
  \item{drtdens}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{drtfill}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{drttotal}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{minpeakMS1}{Mininum number of scans that comprise a peak for MS1 mode.}
  \item{minpeakMS2}{Mininum number of scans that comprise a peak for MS2 mode.}

  \item{recurs}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{weight}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{SB}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{SN}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{minintMS1}{Arguments to be passed from \link[enviPick]{enviPickwrap}. Value for the MS1 mode.}
  \item{minintMS2}{Arguments to be passed from \link[enviPick]{enviPickwrap}. Value for the MS2 mode}

  \item{maxint}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{ion_mode}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{ppm}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
  \item{ended}{Arguments to be passed from \link[enviPick]{enviPickwrap}}
}
\value{
  \item{ms1}{Peak picking for MS1 level}
  \item{ms2}{Peak picking for MS2 level}
  \item{annotation}{Annotated metabolites after interrogation with library. In case a feature appears at more than one retention time, that one with highest intensity is selected.}
  \item{annotationfull}{Matrix with the metabolites annotated.}
  \item{nodetected}{List with the fragments no detected to each metabolite.}
  \item{nodetectedMS2}{List with the fragments no detected to each metabolite considering remaining ion molecular in MS2}
  \item{RawData1}{Raw scans from raw data.}
  \item{RawData2}{Raw scans from raw data.}

}


\author{Manuel D Peris Diaz}
\references{
1. R-MetaboList: a flexible tool for metabolite extraction from high-resolution data-independent acquisition mass spectrometry analysis. Metabolites. Soon

2. A Survey of Orbitrap All Ion Fragmentation Analysis Assessed by an R MetaboList Package to Study Small-Molecule Metabolites. Chromatographia. 2018, 81, 981-994.
}
\examples{
library(MetaboList)

CE.isolation("AIFpos1000-AIF.mzXML","fileposB")

#Reading the database.csv file:
# database<- read.csv("C:/database.csv")

#Processing peak-picking and annotation with default parameters

#aif5<-AIF(fileMS,fileMS2CE5,database,CE=5, ion_mode = "positive")



}
