

Filter_AIF<-function(aif5,full=TRUE,a=0,database){
"Rt(min)"=NULL
peakGG<-data.frame()
Filtered<-data.frame()
if (full=="TRUE"){
Num<-length(levels(aif5$annotationfull$name))

Annotated<-aif5$annotationfull



for(i in 1:Num){
  name<-levels(Annotated$name)[i]
  peakG<-c()
  peakG <- Annotated[which(Annotated$name == name), ]
  nameDatabase<-trimws(name)
  DatabaseMet<-database[which(database$Metabolites == (nameDatabase)),]
  peakG$`Rt(min)` <- round(peakG$`Rt(min)`, digits = 1)
  options(digits=9)
  fragments <- as.numeric(DatabaseMet[, 3:length(DatabaseMet[,])])
  fragments <-  fragments[!is.na(fragments)]
  NFragments<-length(fragments)

  N=NFragments+1

  SubsetMS2<-subset(peakG,peakG$`MS Level`==2)
  if( nrow(SubsetMS2)>0){

  Filtered <-data.frame(peakG %>% group_by(`Rt(min)`) %>% filter(n() >= N-a))

  if(nrow(Filtered)>0){
   peakGG<-rbind(peakGG,Filtered)
  }

}}



  write.csv2(peakGG,file="Filter_AIF.csv")
  return(length(unique(peakGG$name)))

}else{

  Num<-length(levels(aif5$annotation$Metabolite))

  Annotated<-aif5$annotation



  for(i in 1:Num){
    name<-levels(Annotated$Metabolite)[i]
    peakG<-c()
    peakG <- Annotated[which(Annotated$Metabolite == name), ]
    nameDatabase<-trimws(name)
    DatabaseMet<-database[which(database$Metabolites == (nameDatabase)),]
    peakG$`Rt(min)` <- round(peakG$`Rt(min)`, digits = 1)
    options(digits=9)
    fragments <- as.numeric(DatabaseMet[, 3:length(DatabaseMet[,])])
    fragments <-  fragments[!is.na(fragments)]
    NFragments<-length(fragments)


    N=NFragments+1
    SubsetMS2<-subset(peakG,peakG$`MS Level`==2)
    if( nrow(SubsetMS2)>0){

    Filtered <-data.frame(peakG %>% group_by(`Rt(min)`) %>% filter(n() >= N-a))

    if(nrow(Filtered)>0){
      peakGG<-rbind(peakGG,Filtered)
    }
    }
  }



  write.csv2(peakGG,file="Filter_AIF.csv")
  return(length(unique(peakGG$Metabolite)))















  }








}

