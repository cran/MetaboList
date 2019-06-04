



ScoresDIA<-function(input,file,ID1,ID2,CE){

CE<-CE
  peak1 <- input[which(input[,11] == ID1), ]
  peak2 <- input[which(input[,11] == ID2), ]
  total<-rbind(peak1,peak2)

  i<-peak1$Metabolite
  IntensityRatio<-peak2$`Signal Intensity`/peak1$`Signal Intensity`
  AssymetriRatio<- peak2$`Peak Assymetry`/peak1$`Peak Assymetry`

  Peakmax<-total[which.max(total$Area),]
  Peakmin<-total[which.min(total$Area ),]

  peakIDmin<-Peakmin$`MS Level`
  peakIDmax<-Peakmax$`MS Level`


  if (peakIDmin == 2){

    peakmin <- file$RawData2[which(file$RawData2[,7] == Peakmin$PeakID), ]

  }else{
    peakmin <- file$RawData1[which(file$RawData1[,7] == Peakmin$PeakID), ]
    }

  if (peakIDmax == 2){

    peakmax <- file$RawData2[which(file$RawData2[,7] == Peakmax$PeakID), ]

  }else{
    peakmax <- file$RawData1[which(file$RawData1[,7] == Peakmax$PeakID), ]

  }


  peakmin <- as.data.frame(peakmin)
  peakmax <- as.data.frame(peakmax)
  npeakmin <- peakmin[order(peakmin$RT), ]
  npeakmax <- peakmax[order(peakmax$RT), ]

  chrommin <- data.frame(int = npeakmin$intensity, RT = (npeakmin$RT)/60)
  chromax <- data.frame(int = npeakmax$intensity,RT = (npeakmax$RT)/60)

  #### plot coelution

  p3 <- ggplot(data =chromax,aes(x=chromax[,2]  ,y=(chromax[,1])))+
    geom_area(data=chromax,aes(x=chromax[,2],y=(chromax[,1])),fill="#00BFC4",alpha = 0.5,color="darkcyan")+

    geom_area(data=chrommin,aes(x=chrommin[,2],y=(chrommin[,1])),fill="#F8766D",alpha = 0.75,color="tomato3")


  p3b<-p3+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
  plot<-p3b+ scale_x_continuous(name="Retention time (min)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)
  ggsave(plot,filename=paste(CE,"Coelution.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')


  chrommmin <- data.frame()
  chrommin <- data.frame(int = npeakmin$intensity, RT = (npeakmin$RT))
  chromax <- data.frame(int = npeakmax$intensity,RT = (npeakmax$RT))

  chrommmin <- data.frame(chrommin$int, RT = round(chrommin$RT,
                                                   digits = 0))
  chrommmax <- data.frame()
  chrommmax <- data.frame(chromax$int, RT = round(chromax$RT,
                                                  digits = 0))
  nmerged <- data.frame()
  merged <- data.frame()
  merged <- merge(chrommmin, chrommmax, by = "RT")

  ggplotRegression <- function (fit) {



    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +

      stat_smooth(method = "lm", col = "red") +
      labs(title = paste("Cor = ",signif(sqrt(summary(fit1)$`r.squared`), 3),
                         " P =",signif(summary(fit)$coef[2,4], 5)))
  }

    fit1 <- lm(chromax.int~ chrommin.int, data = merged)
    fitt<-ggplotRegression(fit1)
    fittt<-fitt+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
    plot<-fittt+ scale_x_continuous(name="Intensity Peak A",labels= scales::scientific)+ scale_y_continuous(name="Intensity Peak B",labels= scales::scientific)
    ggsave(plot,filename=paste(CE,"Correlation.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')
  plot

   score <- cor(merged[, "chrommin.int"], merged[, "chromax.int"])
  if (is.na(score)== TRUE){

    score<-0
    print("No aligned peaks")}else{


    #   tiff(file=paste(CE,"RTmerged.tiff",sep=""), res = 800, width = 6, height = 4, units = 'in')
      #  plot(merged$RT,merged$chromax.int,type="l",col="red",xlab="RT (s)",ylab="Intensity",cex.lab=1.5,lwd=3,panel.first=grid(lty=1))
      #lines(merged$RT,merged$chrommin.int,col=c("blue"),lwd=3)}


   p3 <- ggplot(data =merged,aes(x=merged$RT  ,y=merged$chromax.int))+
     geom_smooth(data=merged,aes(x=merged$RT  ,y=merged$chromax.int),color="darkcyan",size=1.4,method = "loess",se=F)+
     geom_smooth(data=merged,aes(x=merged$RT  ,y=merged$chrommin.int),color="tomato3",size=1.4,method = "loess",se=F)

   p3b<-p3+theme_light()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),axis.title.y = element_text(size=16,margin=margin(t=0,r=20,b=60,l=0)))
   plot<-p3b+ scale_x_continuous(name="Retention time (s)")+ scale_y_continuous(name="Intensity",labels= scales::scientific)
   ggsave(plot,filename =paste(CE,"RTmerged.tiff",sep=""), dpi = 800, width = 6, height = 4, units = 'in')




   # dev.off()
  return(list(score=score,IntensityRatio=IntensityRatio, AssymetriRatio= AssymetriRatio))
}

}
