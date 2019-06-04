FullMS<-function(file, database, rtw = 10, mzw = 0.0004, dmzgap = 50,
                       dmzdens = 20, drtgap = 25, drtsmall = 50, drtdens = 20, drtfill = 5,
                       drttotal = 100, minpeak = 5, recurs = 3, weight = 1,
                       SB = 2, SN = 1.5, minint = 1000, maxint = 9e+09, ion_mode = "positive",
                       ppm = TRUE,ended=6)
{


  if (is.numeric(mzw) == FALSE) {
    stop("Please, establish a numeric m/z window ")
  }


  ms1 <- enviPickwrap(file, MSlevel = c(1),dmzgap ,ppm=ppm ,dmzdens,
                      drtdens,drtgap,drtsmall,drtfill ,drttotal ,
                      minpeak,recurs,weight , SB, SN,minint , maxint ,ended,ion_mode = ion_mode)

   b <- as.data.frame(ms1[[8]])
  f<-mzw
  rtw<-rtw
  msbc=ms=name=X=fullms=rt=fullms2=rt2=to=fullmsset=msb=ions1bb= ppmdif=NULL
  ms <- data.frame()

  databasemzRT<-subset(database,!is.na(database$RT))


  if (is.numeric(rtw) == TRUE){
      for (i in 1:nrow(databasemzRT)){
      from<-5
      to<-as.numeric(length(databasemzRT[i,]))
      fullmsset <- databasemzRT[i, from:to]
      options(digits=9)
      fullms2 <-  as.numeric(fullmsset[!is.na(fullmsset)])
      rt2 <- as.numeric(as.character(databasemzRT[i, 2]))
      name<- databasemzRT[i, 4]
      b<-cbind(b,Isotope="Isotope")
      b<-cbind(b,ppmdif="ppmdif")

      msbc<-data.frame()
      for (i in 1:length(fullms2)){


          msb <- subset(b, b[, 1] > fullms2[i] - f & b[, 1] < fullms2[i] + f & b[,5]> rt2*60 - rtw & b[,5]< rt2*60 + rtw)

          if (nrow(msb) > 0){
          ppmdif<-((as.numeric(msb[1])-fullms2[i])/fullms2[i])*1000000
          msb[,13]<- ppmdif
          msb[,12]<-names(fullmsset)[i]
         }
          msbc<-rbind(msb,msbc)
      }

          if (nrow(msbc) == 0)  {
         next
             return(print(paste("No MS1 found for this metabolite, remove it from database or change algorithms parameters:",
                           name)))}


          msbc<-as.data.frame(msbc)


      ions1bb <- data.frame(name, msbc[, 1],msbc[, 3],msbc[,4],as.numeric(as.character(msbc[, 5]))/60
                            , as.numeric(as.character(msbc[, 7]))/60 -  as.numeric(as.character(msbc[, 6]))/60,
                            ((as.numeric(as.character(msbc[, 7]))/60) -  (as.numeric(as.character(msbc[, 5]))/60))/((as.numeric(as.character(msbc[, 5]))/60)-(as.numeric(as.character(msbc[, 6]))/60)),
                            msbc[,10],msbc[,12],msbc[,13])
      ms <- rbind(ions1bb, ms)
      }

    } else{

        for (i in 1:nrow(database)){
          from<-5
          to<-as.numeric(length(database[i,]))
          fullmsset <- database[i, from:to]
          options(digits=9)
          fullms2 <-  as.numeric(fullmsset[!is.na(fullmsset)])
          rt2 <- as.numeric(as.character(database[i, 2]))
          name<- database[i, 4]
          b<-cbind(b,Isotope="Isotope")
          b<-cbind(b,ppmdif="ppmdif")

          msbc<-data.frame()
          for (i in 1:length(fullms2)){
            msb <- subset(b, b[, 1] > fullms2[i] - f & b[, 1] < fullms2[i] + f )
            msb<-as.matrix(msb)
            if (nrow(msb) > 0){
              if (nrow(msb)>1){

                for (j in 1:nrow(msb)){
                  ppmdif<-((as.numeric(msb[j,1])-fullms2[i])/fullms2[i])*1000000
                  msb[j,13]<- ppmdif
                  msb[j,12]<-names(fullmsset)[i]
                  }
                } else{
              ppmdif<-((as.numeric(msb[1])-fullms2[i])/fullms2[i])*1000000
              msb[,13]<- ppmdif
              msb[,12]<-names(fullmsset)[i]} }

            #msb<-as.data.frame(msb)
            msbc<-rbind(msb,msbc)
            }



          if (nrow(msbc) == 0)  {
            next
            return(print(paste("No MS1 found for this metabolite, remove it from database or change algorithms parameters:",
                               name)))

          }
          msbc<-as.data.frame(msbc)


          ions1bb <- data.frame(name, as.numeric(as.character(msbc[,1])),as.numeric(as.character(msbc[,3])),as.numeric(as.character(msbc[,4])),as.numeric(as.character(msbc[, 5]))/60,

                                 as.numeric(as.character(msbc[, 7]))/60 -  as.numeric(as.character(msbc[, 6]))/60,

                                ((as.numeric(as.character(msbc[, 7]))/60) -  (as.numeric(as.character(msbc[, 5]))/60))/((as.numeric(as.character(msbc[, 5]))/60)-(as.numeric(as.character(msbc[, 6]))/60)),

                                as.numeric(as.character(msbc[,10])),(as.character(msbc[,12])),as.numeric(as.character(msbc[,13])))
          ms <- rbind(ions1bb, ms)
          }
      }



  if(nrow(ms)==0){

    return(print("No MS1 found for this library"))
  }else{

    colnames(ms) <- c("name","Ion Metabolite (m/z)","Maximum Intensity","Area","Rt(min)",
                      "Interval Peak (min)","Peak Assymetry","peakID","Isotope" ,"Error(ppm)")

  }
  rownames(ms) <- c()

  results<- list(RawData1=ms1$Scans[[2]],ms=ms,Peaklist=ms1$Peaklist,PP=ms1)

  write.csv(ms,paste0("FullMSannotated",ion_mode,".csv"))
  return(results)
}
