


plot_EIC<-function(fullms,peakID=333,ms=1,CE=0){

  if(ms==1){

  if (is.numeric(peakID) == FALSE) {
    for (i in fullms$ms$peakID){
      peak <- fullms$RawData1[which(fullms$RawData1[,7] == i), ]

      peak<-as.data.frame(peak)
      npeak<-peak[order(peak$RT),]

      p3 <- ggplot(data =npeak,aes(npeak$RT/60  ,y=npeak$intensity))+
        geom_area(data=npeak,aes(x=npeak$RT/60,y=npeak$intensity),fill="#00BFC4",alpha = 0.5,color="darkcyan")
        p3b<-p3+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
        plot<-p3b+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
        ggtitle(paste("Raw:",fullms$ms$name[fullms$ms$peakID==i]))
        ggsave(plot,filename=paste(i,"Raw.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')

      npeak<-as.data.frame(npeak)
    #  predict<-cbind(npeak,predict(fit))
      #fit.bs<-lm(npeak$intensity~ bs(npeak$RT/60))
      #lines(npeak$RT/60,predict(fit.bs,data.frame(x=npeak$RT/60)),col=2,lwd=2)
      npeakRT<-npeak$RT/60

      p3 <- ggplot(data =npeak,aes(x=npeakRT,y=npeak$intensity))+ geom_point()+
                     geom_smooth(method = "loess",se=F)+theme_light()+
                     geom_ribbon(aes(ymin=0,ymax=predict(loess(npeak$intensity~npeakRT))),fill="#00BFC4",alpha = 0.5)+
                  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
      plot<-p3+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
                 ggtitle(paste("Smooth:",fullms$ms$name[fullms$ms$peakID==i]))
      ggsave(plot,filename=paste(i,"Smooth.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')



      diff<-diff(npeak$`m/z`)*100

      tiff(filename=paste(i,"Scan.pdf",sep=""), res = 800, width = 6, height = 4, units = 'in')

      plot(diff,type="h",xlab="Scan",ylab="m/z difference",cex.lab=1.5,lwd=5,col="cyan3",panel.first=grid(lty=1))
      abline(h=mean(diff),col="#3366FF",lwd=3, lty=2)
      abline(h=0,col="black")
      title(main=paste("QC m/z scan:",fullms$ms$name[fullms$ms$peakID==i]))

      dev.off()

    }}else {

    i=peakID
    peak <- fullms$RawData[which(fullms$RawData[,7] == i), ]
    peak<-as.data.frame(peak)
    npeak<-peak[order(peak$RT),]
    npeakRT<-npeak$RT/60

    p3 <- ggplot(data =npeak,aes(npeak$RT/60  ,y=npeak$intensity))+
      geom_area(data=npeak,aes(x=npeak$RT/60,y=npeak$intensity),fill="#00BFC4",alpha = 0.5,color="darkcyan")
    p3b<-p3+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
    plot<-p3b+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
      ggtitle(paste("Raw:",fullms$ms$name[fullms$ms$peakID==i]))
    ggsave(plot,filename=paste(i,"Raw.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')



    p3 <- ggplot(data =npeak,aes(x=npeakRT,y=npeak$intensity))+ geom_point()+
      geom_smooth(method = "loess",se=F)+theme_light()+
      geom_ribbon(aes(ymin=0,ymax=predict(loess(npeak$intensity~npeakRT))),fill="#00BFC4",alpha = 0.5)+
      theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
    plot<-p3+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
      ggtitle(paste("Smooth:",fullms$ms$name[fullms$ms$peakID==i]))
    ggsave(plot,filename=paste(i,"Smooth.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')

    diff<-diff(npeak$`m/z`)*100

    tiff(filename=paste(i,"Scan.tiff",sep=""), res = 800, width = 6, height = 4, units = 'in')
    plot(diff,type="h",xlab="Scan",ylab="m/z difference",cex.lab=1.5,lwd=5,col="cyan3",panel.first=grid(lty=1))
    abline(h=mean(diff),col="#3366FF",lwd=3, lty=2)
    abline(h=0,col="black")
    title(main=paste("QC m/z scan:",fullms$ms$name[fullms$ms$peakID==i]))

    dev.off()
  }}else{

    if(CE==0){
    ## here is for ms1 from DIA experiment

    if (is.numeric(peakID) == FALSE) {
      for (i in fullms$annotation$PeakID){
        peak <- fullms$RawData1[which(fullms$RawData1[,7] == i), ]

        peak<-as.data.frame(peak)
        npeak<-peak[order(peak$RT),]

        p3 <- ggplot(data =npeak,aes(npeak$RT/60  ,y=npeak$intensity))+
          geom_area(data=npeak,aes(x=npeak$RT/60,y=npeak$intensity),fill="#00BFC4",alpha = 0.5,color="darkcyan")
        p3b<-p3+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
        plot<-p3b+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
          ggtitle(paste("Raw:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))
        ggsave(plot,filename=paste(i,"Raw.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')

        npeak<-as.data.frame(npeak)
       # predict<-cbind(npeak,predict(fit))
        #fit.bs<-lm(npeak$intensity~ bs(npeak$RT/60))
        #lines(npeak$RT/60,predict(fit.bs,data.frame(x=npeak$RT/60)),col=2,lwd=2)
        npeakRT<-npeak$RT/60

        p3 <- ggplot(data =npeak,aes(x=npeakRT,y=npeak$intensity))+ geom_point()+
          geom_smooth(method = "loess",se=F)+theme_light()+
          geom_ribbon(aes(ymin=0,ymax=predict(loess(npeak$intensity~npeakRT))),fill="#00BFC4",alpha = 0.5)+
          theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
        plot<-p3+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
          ggtitle(paste("Smooth:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))
        ggsave(plot,filename=paste(i,"Smooth.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')

        #

        diff<-diff(npeak$`m/z`)*100

        pdf(file=paste(i,"Scan.pdf",sep=""))
        plot(diff,type="h",xlab="Scan",ylab="m/z difference",cex.lab=1.5,lwd=5,col="cyan3",panel.first=grid(lty=1))
        abline(h=mean(diff),col="#3366FF",lwd=3, lty=2)
        abline(h=0,col="black")
        title(main=paste("QC m/z scan:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))

        dev.off()

      }}

    else {

      i=peakID
      peak <- fullms$RawData1[which(fullms$RawData1[,7] == i), ]
      peak<-as.data.frame(peak)
      npeak<-peak[order(peak$RT),]
      npeakRT<-npeak$RT/60

      p3 <- ggplot(data =npeak,aes(npeak$RT/60  ,y=npeak$intensity))+
        geom_area(data=npeak,aes(x=npeak$RT/60,y=npeak$intensity),fill="#00BFC4",alpha = 0.5,color="darkcyan")
      p3b<-p3+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
      plot<-p3b+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
        ggtitle(paste("Raw:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))
      ggsave(plot,filename=paste(i,"Raw.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')

      npeak<-as.data.frame(npeak)
      #predict<-cbind(npeak,predict(fit))
      #fit.bs<-lm(npeak$intensity~ bs(npeak$RT/60))
      #lines(npeak$RT/60,predict(fit.bs,data.frame(x=npeak$RT/60)),col=2,lwd=2)
      npeakRT<-npeak$RT/60

      p3 <- ggplot(data =npeak,aes(x=npeakRT,y=npeak$intensity))+ geom_point()+
        geom_smooth(method = "loess",se=F)+theme_light()+
        geom_ribbon(aes(ymin=0,ymax=predict(loess(npeak$intensity~npeakRT))),fill="#00BFC4",alpha = 0.5)+
        theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
      plot<-p3+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
        ggtitle(paste("Smooth:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))
      ggsave(plot,filename=paste(i,"Smooth.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')

      #

      diff<-diff(npeak$`m/z`)*100

      pdf(file=paste(i,"Scan.pdf",sep=""))
      plot(diff,type="h",xlab="Scan",ylab="m/z difference",cex.lab=1.5,lwd=5,col="cyan3",panel.first=grid(lty=1))
      abline(h=mean(diff),col="#3366FF",lwd=3, lty=2)
      abline(h=0,col="black")
      title(main=paste("QC m/z scan:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))

      dev.off()



  }

    }else{



      if (is.numeric(peakID) == FALSE) {
        for (i in fullms$annotation$PeakID){
          peak <- fullms$RawData2[which(fullms$RawData2[,7] == i), ]

          peak<-as.data.frame(peak)
          npeak<-peak[order(peak$RT),]


          p3 <- ggplot(data =npeak,aes(npeak$RT/60  ,y=npeak$intensity))+
            geom_area(data=npeak,aes(x=npeak$RT/60,y=npeak$intensity),fill="#00BFC4",alpha = 0.5,color="darkcyan")
          p3b<-p3+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
          plot<-p3b+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
            ggtitle(paste("Raw:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))
          ggsave(plot,filename=paste(i,"Raw.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')

          npeak<-as.data.frame(npeak)
          # predict<-cbind(npeak,predict(fit))
          #fit.bs<-lm(npeak$intensity~ bs(npeak$RT/60))
          #lines(npeak$RT/60,predict(fit.bs,data.frame(x=npeak$RT/60)),col=2,lwd=2)
          npeakRT<-npeak$RT/60

          p3 <- ggplot(data =npeak,aes(x=npeakRT,y=npeak$intensity))+ geom_point()+
            geom_smooth(method = "loess",se=F)+theme_light()+
            geom_ribbon(aes(ymin=0,ymax=predict(loess(npeak$intensity~npeakRT))),fill="#00BFC4",alpha = 0.5)+
            theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
          plot<-p3+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
            ggtitle(paste("Smooth:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))
          ggsave(plot,filename=paste(i,"Smooth.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')

          #

          diff<-diff(npeak$`m/z`)*100

          pdf(file=paste(i,"Scan.pdf",sep=""))
          plot(diff,type="h",xlab="Scan",ylab="m/z difference",cex.lab=1.5,lwd=5,col="cyan3",panel.first=grid(lty=1))
          abline(h=mean(diff),col="#3366FF",lwd=3, lty=2)
          abline(h=0,col="black")
          title(main=paste("QC m/z scan:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))

          dev.off()

        }}

      else {

        i=peakID
        peak <- fullms$RawData2[which(fullms$RawData2[,7] == i), ]
        peak<-as.data.frame(peak)
        npeak<-peak[order(peak$RT),]
        npeakRT<-npeak$RT/60


        p3 <- ggplot(data =npeak,aes(npeak$RT/60  ,y=npeak$intensity))+
          geom_area(data=npeak,aes(x=npeak$RT/60,y=npeak$intensity),fill="#00BFC4",alpha = 0.5,color="darkcyan")
        p3b<-p3+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
        plot<-p3b+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
          ggtitle(paste("Raw:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))
        ggsave(plot,filename=paste(i,"Raw.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')
        plot
        npeak<-as.data.frame(npeak)
        #predict<-cbind(npeak,predict(fit))
        #fit.bs<-lm(npeak$intensity~ bs(npeak$RT/60))
        #lines(npeak$RT/60,predict(fit.bs,data.frame(x=npeak$RT/60)),col=2,lwd=2)
        npeakRT<-npeak$RT/60

        p3 <- ggplot(data =npeak,aes(x=npeakRT,y=npeak$intensity))+ geom_point()+
          geom_smooth(method = "loess",se=F)+theme_light()+
          geom_ribbon(aes(ymin=0,ymax=predict(loess(npeak$intensity~npeakRT))),fill="#00BFC4",alpha = 0.5)+
          theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
        plot<-p3+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)+
          ggtitle(paste("Smooth:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))
        ggsave(plot,filename=paste(i,"Smooth.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')

        #

        diff<-diff(npeak$`m/z`)*100

        pdf(file=paste(i,"Scan.pdf",sep=""))
        plot(diff,type="h",xlab="Scan",ylab="m/z difference",cex.lab=1.5,lwd=5,col="cyan3",panel.first=grid(lty=1))
        abline(h=mean(diff),col="#3366FF",lwd=3, lty=2)
        abline(h=0,col="black")
        title(main=paste("QC m/z scan:",fullms$annotation$Metabolite[fullms$annotation$PeakID==i]))

        dev.off()




  }
    }}}


