\name{annotation}
\alias{annotation}
\title{Automatic Metabolite Annotation from LC-MS experiments.}
\description{Automatic metabolite annotation through filters after processed the LC-MS data with enviPick package.}
\usage{annotation(file, database, rtw = 3, mzw = 0.004,
dmzgap = 20, dmzdens = 8, drtgap = 1000,
drtsmall = 5, drtdens = 5, drtfill = 10,
drttotal = 200, minpeak = 2, recurs = 200,
weight = 0.25, SB = 1.5, SN = 0.5,
minint = 10, maxint = 9e+09, ion_mode = "positive",
ppm = TRUE, isobaric = FALSE)}
\arguments{
  \item{file}{
   A mzXML file from the LC-MS experiment in positive or negative ionization.
}
  \item{database}{
  A file with data arranged in columns including the names in the first row: Metabolite; Mass molecular ion; Mass Fragment 1; Mass Fragment 2; Mass Fragment 3...
  }
  \item{rtw}{
  numeric. The difference between the theoretical retention time value and the experimental. Default value=3
  }
  \item{mzw}{
  numeric. The difference between the theoretical m/z value and the experimental. Default value=0.004
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
  \item{isobaric}{Annotation of isobaric ions. The database contains the isobaric ions (ms) and fragments. FALSE by default.}
}
\value{
  \item{annotation}{It contains the metabolites annotated with the following columns: Metabolite name; Ion metabolite; retention time (Rt) (min); Interval Peak (min); Area; Maximum intensity}
  \item{summary}{List with the fragments no detected to each metabolite.}
  \item{ms}{Matrix with the pick-picking performed with enviPick function to LC-MS data}
  \item{ms2}{Matrix with the pick-picking performed with enviPick function to LC-MS/MS data}}
\author{Manuel David Peris Diaz}
\examples{
library(MetaboList)
#Reading the file.mzXML
# file<-Cells_pos.mzXML
#Reading the database.csv file:
# database<- read.csv("C:/CellsPosDatabase.csv")
#Processing peak-picking and annotation with default parameters
# processMS<-annotate(file,database)
#Output:
#processMS$annotation
# processMS$nodetected
}
