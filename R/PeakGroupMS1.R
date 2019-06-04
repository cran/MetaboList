

PeakGroupMS1<- function (annotation_pos)
{

  f123<-c()
  Lista<-c()
  f123<-annotation_pos
    f123<-f123[order(f123$name),]
    
    for(i in 1:length(unique(f123$name))){
      
      name<-unique(f123$name)[i]
      
      peakG<-c()
      peakG <- f123[which(f123$name == name), ]
      peakGF<-unique(peakG)
      
      Lista[i]<-list(peakGF)
      write.csv2(peakGF,paste0("PeakGroupMS1",name,".csv"))
      

    }
    return(Lista)
    
  }
  
  
