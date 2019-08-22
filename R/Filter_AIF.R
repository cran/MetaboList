

Filter_AIF<-function(aif5,a=1,database,rttol=0.2){
"Rt(min)"=NULL
peakGG<-data.frame()


  Num<-length(levels(aif5$annotationfull$name))

   Annotated<-aif5$annotationfull


peakGGG<-data.frame()
for(i in 1:Num){
  name<-levels(Annotated$name)[i]
  peakG<-c()
  peakG <- Annotated[which(Annotated$name == name), ]
  nameDatabase<-trimws(name)
  DatabaseMet<-database[which(database$Metabolites == (nameDatabase)),]
  peakG$`Rt(min)` <- round(peakG$`Rt(min)`, digits = 2)
  options(digits=9)
  fragments <- as.numeric(DatabaseMet[, 3:length(DatabaseMet[,])])
  fragments <-  fragments[!is.na(fragments)]
  NFragments<-length(fragments)

  MS1<-subset(peakG,peakG$`MS Level`==1)
  MS2<-subset(peakG,peakG$`MS Level`==2)
  rttol<-rttol


  if( nrow(MS2)>0& nrow(MS1)>0){
    if(a==1){
      NMS2_matching_MS1<-1}
    if(a==0){
      NMS2_matching_MS1<-NFragments}

       peakGG<-data.frame()
          for (i in 1:nrow(MS1)){

            matches<-data.frame()
            if (sum(abs(MS1$`Rt(min)`[i] - as.numeric(MS2$`Rt(min)`)) <
                    rttol) >= NMS2_matching_MS1){
              rowsMS2 <-MS2[(abs(MS1$`Rt(min)`[i] - as.numeric(MS2$`Rt(min)`)) <= rttol),]
              rowsMS1 <-MS1[i,]
              matches<-rbind(rowsMS1,rowsMS2) }

            if(nrow(matches)>0){
              peakGG<-rbind(peakGG,matches)}

          }
  }

  peakGGG<-rbind(peakGG, peakGGG)

  }
FinalResult<-peakGGG[!duplicated(peakGGG), ]


  write.csv2(FinalResult,file="Filter_AIF.csv")
  return(length(unique(FinalResult$name)))


}
