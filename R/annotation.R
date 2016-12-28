

annotation<-function(file,database,rtw=3,mzw=0.004,dmzgap=20, dmzdens=8,
                    drtgap=1000, drtsmall=5, drtdens=5, drtfill=10, drttotal=200, minpeak=2,
                    recurs=200, weight=0.25, SB=1.5, SN=0.5, minint=1E+01, maxint=9E+09,
                    ion_mode="positive",ppm=TRUE,isobaric=FALSE){

  if(is.numeric(rtw)==FALSE){
    stop("Please, establish a numeric retention time window ")
  }
  if(is.numeric(mzw)==FALSE){
    stop("Please, establish a numeric retention time window ")
  }


    ms1<-enviPickwrap(file,MSlevel=c(1), dmzgap, dmzdens,ppm,
                     drtgap, drtsmall, drtdens, drtfill, drttotal, minpeak,
                     recurs, weight, SB, SN, minint, maxint,
                     ion_mode=ion_mode)

    ms2<-enviPickwrap(file,MSlevel=c(2), dmzgap, dmzdens,ppm,
                      drtgap, drtsmall, drtdens, drtfill, drttotal, minpeak,
                      recurs, weight, SB, SN, minint, maxint,
                      ion_mode=ion_mode)
    database<-database
    dff<-data.frame()
    df2<-data.frame()
    summary2<-list()
    summary<-list()
    b<-ms1[[8]]
    a<-ms2[[8]]
    f<-mzw
    if(isobaric==TRUE){
    for (i in 1:length(database[,2])){
      fullms<-database[i,2]
      fragments<-as.numeric(database[i,3:length(database[i,])])
      name<-database[i,1]
      cc<-as.numeric(database[i,2:length(database[i,])])
      ms<-b[which(b[,1]>fullms-f & b[,1]<fullms+f),]
      x<-data.frame()
      aiflist<-data.frame()
      for (i in 1:length(fragments)){
        x<-a[which(a[,1] >as.numeric(fragments[i]-f) & a[,1]< as.numeric(fragments[i]+f)),]
        aiflist<-rbind(aiflist,x)
      }
      xy<-data.frame()
      xy<-rbind(ms,aiflist)
      xx<-data.frame(xy[,1],xy[,5],xy$minRT,xy$maxRT,xy$max_int,xy$sum_int,xy[,5]/60)
      trrr<-xx
      relative<-data.frame()
      X<-data.frame()
      ions<-data.frame()
      for(i in 1:length(ms[,1])){
        ions<- trrr[which(trrr[i,2]-trrr[,2]>-rtw &trrr[i,2]-trrr[,2]<rtw),]
        relative<-(ions$xy.max_int)/(max(ions$xy.max_int))*100
        X<-data.frame(ions[,1],ions[,2]/60,ions[,4]/60-ions$xy.minRT/60,(ions[,4]-ions[,2])/(ions[,2]-ions$xy.minRT),ions$xy.sum_int,ions$xy.max_int,relative)
        colnames(X) <- c("Ion Metabolite (m/z)", "Rt(min)",  "Interval Peak (min)","Assymetry (B/A)","Area","Maximum Intensity","Relative Intensity (%)")
        ions<-data.frame()
        x<-data.frame()
        ions1bb<-data.frame()
        for (i in 1:length(X[,1])){
          ions<- X[which(X[i,1]-X[,1]>-1 & X[i,1]-X[,1]<1),]
          ions1bb<-data.frame(paste('',name),mean(ions$`Ion Metabolite (m/z)`),mean(ions$`Rt(min)`),
                              mean(ions$`Interval Peak (min)`),sum(ions$Area),sum(ions$`Maximum Intensity`))
          x<-rbind(ions1bb,x)
        }
        colnames(x) <- c("name","Ion Metabolite (m/z)", "Rt(min)",  "Interval Peak (min)","Area","Maximum Intensity")
        ionst<-aggregate(.~`Ion Metabolite (m/z)`,data=x,FUN=mean)
        jk<-ionst$`Ion Metabolite (m/z)`
        jkr<-round(jk, digits=1)
        ccr<-round(cc,digits=1)
        dif<-(setdiff(ccr,jkr))
        xf<-data.frame(paste('',name),ionst)
        xf<-xf[,-3]
        dd<-data.frame(xf)
        dff<-rbind(dd,dff)
        summary2<-list(NoDetected=c(dif))
        summary[[paste('',name)]]<-summary2
      }}
    colnames(dff) <- c("Metabolite","Ion Metabolite (m/z)", "Rt(min)",  "Interval Peak (min)","Area","Maximum Intensity")
    annotation<-dff
    results<-list(ms1=ms1[[8]],ms2=ms2[[8]],annotation=annotation,nodetected=summary)
    return(results)}
  else{
    for (i in 1:length(database[,2])){
      fullms<-database[i,2]
      fragments<-as.numeric(database[i,3:length(database[i,])])
      name<-database[i,1]
      cc<-as.numeric(database[i,2:length(database[i,])])
      ms<-b[which(b[,1]>fullms-f & b[,1]<fullms+f),]
      x<-data.frame()
      aiflist<-data.frame()
      for (i in 1:length(fragments)){
        x<-a[which(a[,1] >fragments[i]-f& a[,1]< fragments[i]+f),]
        aiflist<-rbind(aiflist,x)
        colnames(aiflist) <- c("m/z", "var_m/z",  "max_int","sum_int","RT","minRT","maxRT","part_ID","EIC_ID","peak_ID","Score")
      }
      xy<-data.frame()
      xy<-rbind(ms,aiflist)
      xx<-data.frame(xy[,1],xy[,5],xy$minRT,xy$maxRT,xy$max_int,xy$sum_int,xy[,5]/60)
      trrr<-xx
      ions<- trrr[which(trrr[1,2]-trrr[,2]>-rtw &trrr[1,2]-trrr[,2]<rtw),]
      relative<-(ions$xy.max_int)/(max(ions$xy.max_int))*100
      X<-data.frame(ions[,1],ions[,2]/60,ions[,4]/60-ions$xy.minRT/60,(ions[,4]-ions[,2])/(ions[,2]-ions$xy.minRT),ions$xy.sum_int,ions$xy.max_int,relative)
      colnames(X) <- c("Ion Metabolite (m/z)", "Rt(min)",  "Interval Peak (min)","Assymetry (B/A)","Area","Maximum Intensity","Relative Intensity (%)")
      ions<-data.frame()
      x<-data.frame()
      ions1bb<-data.frame()
      for (i in 1:length(X[,1])){
        ions<- X[which(X[i,1]-X[,1]>-1 & X[i,1]-X[,1]<1),]
        ions1bb<-data.frame(paste('',name),mean(ions$`Ion Metabolite (m/z)`),mean(ions$`Rt(min)`),
                            mean(ions$`Interval Peak (min)`),sum(ions$Area),sum(ions$`Maximum Intensity`))
        x<-rbind(ions1bb,x)
      }
      colnames(x) <- c("name","Ion Metabolite (m/z)", "Rt(min)",  "Interval Peak (min)","Area","Maximum Intensity")
      ionst<-aggregate(.~`Ion Metabolite (m/z)`,data=x,FUN=mean)
      jk<-ionst$`Ion Metabolite (m/z)`
      jkr<-round(jk, digits=1)
      ccr<-round(cc,digits=1)
      dif<-(setdiff(ccr,jkr))
      xf<-data.frame(paste('',name),ionst)
      xf<-xf[,-3]
      dd<-data.frame(xf)
      dff<-rbind(dd,dff)
      summary2<-list(NoDetected=c(dif))
      summary[[paste('',name)]]<-summary2
    }

    colnames(dff) <- c("Metabolite","Ion Metabolite (m/z)", "Rt(min)",  "Interval Peak (min)","Area","Maximum Intensity")
    annotationn<-dff
    results<-list(ms1=ms1[[8]],ms2=ms2[[8]],annotation=annotationn,nodetected=summary)
    return(results)
  }
  }
