
PeakGroup<-function (aif1,aif2,aif3)
{
  if(missing(aif3)){
  file1<-aif1$annotation
  file2<-aif2$annotation
  file1Raw1<-aif1$RawData1
  file1Raw2<-aif1$RawData2
  file2Raw1<-aif2$RawData1
  file2Raw2<-aif2$RawData2
  }else{
    file1<-aif1$annotation
    file2<-aif2$annotation
    file3<-aif3$annotation

    f123<-c()
   f123<-rbind(file1,file2,file3)
   f123<-f123[order(f123$Metabolite),]

for(i in 1:length(levels(f123$Metabolite))){

  name<-levels(f123$Metabolite)[i]

  peakG<-c()
  peakG <- f123[which(f123$Metabolite == name), ]
  peakG<-peakG[,-8]

  ppmdif<-((as.numeric(peakG[,2])-peakG[,13])/peakG[,13])*1000000
  peakGG<-cbind(peakG,ppmdif)
  peakGF<-unique(peakGG)
  write.csv(peakGF,paste0("PeakGroup",name,".csv"))


}

  }


}
