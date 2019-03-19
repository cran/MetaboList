

ScoresMS1<-function(fullmspos,fullmsposb,ID1,ID2){


  annotated<-fullmspos$ms
  annotatedb<-fullmsposb$ms


  peak1 <- annotated[which(annotated[,8] == ID1), ]
  peak2 <- annotatedb[which(annotatedb[,8] == ID2), ]
  i<-peak1$name
  duplicated<-rbind(peak1,peak2)
  #duplicated<-annotated[annotated$name %in% annotated$name[duplicated(annotated$name)],]

  # max(as.numeric(as.vector(duplicated[,2])))-min(as.numeric(as.vector(duplicated[,2])))

   last<-duplicated[which.max(as.numeric(as.vector(duplicated[,2]))),]
   first<-duplicated[which.min(as.numeric(as.vector(duplicated[,2]))),]

    IntensityRatio<- as.numeric(as.vector(first$`Maximum Intensity`))/as.numeric(as.vector(last$`Maximum Intensity`))
    AssymetriRatio<- first$`Peak Assymetry`/last$`Peak Assymetry`

    Peakmax<-duplicated[which.max(duplicated$`Ion Metabolite (m/z)`),]
    Peakmin<-duplicated[which.min(duplicated$`Ion Metabolite (m/z)`),]

    peakIDmin<-Peakmin$peakID
    peakIDmax<-Peakmax$peakID

    peakmin <- fullmspos$RawData[which(fullmspos$RawData[,7] == peakIDmin), ]
    peakmax <- fullmsposb$RawData[which(fullmsposb$RawData[,7] == peakIDmax), ]

     peakmin<-as.data.frame( peakmin)
    peakmax<-as.data.frame(peakmax)
    npeakmin<-peakmin[order(peakmin$RT),]
    npeakmax<-peakmax[order(peakmax$RT),]

    fitmin<-smooth.spline(npeakmin$RT/60, npeakmin$intensity)
    chrommin<- data.frame(int = fitmin$y, RT = npeakmin$RT/60)

    fitmax<-smooth.spline(npeakmax$RT/60, npeakmax$intensity)
    chromax<- data.frame(int = fitmax$y, RT = npeakmax$RT/60)

    #### plot coelution


    p3 <- ggplot(data =chromax,aes(x=chromax[,2]  ,y=(chromax[,1])))+

      geom_area(data=chrommin,aes(x=chrommin[,2],y=(chrommin[,1])),fill="#F8766D",alpha = 0.75,color="tomato3")+

    geom_area(data=chromax,aes(x=chromax[,2],y=(chromax[,1])),fill="#00BFC4",alpha = 0.5,color="darkcyan")

     p3b<-p3+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
    plot<-p3b+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)
    ggsave(plot,filename=paste(i,"Coelution.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')



    chrommin<- data.frame(int = fitmin$y, RT = npeakmin$RT)
    chromax<- data.frame(int = fitmax$y, RT = npeakmax$RT)
    chrommmin<-data.frame()
    chrommmin<-data.frame(chrommin$int,RT=round(chrommin$RT,digits = 0))
    chrommmax<-data.frame()
    chrommmax<-data.frame(chromax$int,RT=round(chromax$RT,digits = 0))
    nmerged<-data.frame()
    merged<-data.frame()
    merged <- merge(chrommmin, chrommmax, by = "RT")

    score <- cor(merged[, "chrommin.int"], merged[, "chromax.int"])
    if (is.na(score)== TRUE){

      score<-0
      print("No aligned peaks")}else{
        plot(merged$RT,log(merged$chromax.int),type="l",col="red")
        lines(merged$RT,log(merged$chrommin.int),col="blue")}
    return(list(score=score,IntensityRatio=IntensityRatio, AssymetriRatio= AssymetriRatio, name=duplicated$name))
    }

