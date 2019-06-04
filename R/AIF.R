




AIF<-function(fileMS1,fileMS2,CE=0, database, rtw = 7, mzw = 0.005, dmzgap = 50,
              drtdens = 20, drtgap = 25, drtsmallMS1 = 100,drtsmallMS2=30, dmzdensMS1 = 15,dmzdensMS2=30, drtfill = 5,
              drttotal = 100, minpeakMS1 = 5,minpeakMS2=3, recurs = 2, weight = 2,
              SB = 3, SN = 2, minintMS1 = 1000,minintMS2=100, maxint = 9e+09, ion_mode = "positive",
              ppm = TRUE,ended=6)
{
  if (is.numeric(rtw) == FALSE) {
    stop("Please, establish a numeric retention time window ")
  }
  if (is.numeric(mzw) == FALSE) {
    stop("Please, establish a numeric m/z window ")
  }

  database<-as.data.frame(database)


  f <- mzw
  rtw <- rtw





  ms1 <- enviPickwrap(fileMS1, MSlevel = c(1), dmzgap, dmzdens=dmzdensMS1,
                      ppm=TRUE, drtgap, drtsmall=drtsmallMS1, drtdens, drtfill, drttotal, minpeak=minpeakMS1,
                      recurs, weight, SB, SN, minint=minintMS1, maxint, ion_mode = ion_mode,ended)

  ms2 <- enviPickwrap(fileMS2, MSlevel = c(2), dmzgap, dmzdens=dmzdensMS2,
                      ppm=TRUE, drtgap, drtsmall=drtsmallMS2, drtdens, drtfill, drttotal, minpeak=minpeakMS2,
                      recurs, weight, SB, SN, minint=minintMS2, maxint, ion_mode = ion_mode,ended)


  ions1bb=NULL
  dd=NULL
  dff=NULL
  xf=NULL
  MS1MS2=NULL
  MS2_f=NULL
  MS1_f=NULL
  MS1xx=NULL
  ionst =NULL
  x=NULL
  xx=NULL
  ions=NULL
  Xc=NULL
  jk_ms1=NULL
  jk_ms2=NULL
  jkr_ms2=NULL
  jkr=NULL
  relative=NULL
  rtmetab=NULL
  app_F=NULL
  cc=NULL
  ccr=NULL
  ccr2=NULL
  dif=NULL
  dif_ms2=NULL
  annotationb =NULL
  i=NULL
  ms=NULL
  ms_2=NULL
  ms1anot=NULL
  fragments=NULL
  trrr=NULL
  X=NULL
  MS2xxff=NULL
  MS1xxff=NULL
  results=NULL
  xxf=NULL
  xy=NULL
  metab<- data.frame()
  ionst_ms1_xxff<- data.frame()
  ionst_ms2_xxff<- data.frame()
  ions<- data.frame()
  ions1bb<- data.frame()
  ionst<- data.frame()
  aiflist <- data.frame()
  annotation <- data.frame()
  annotation_final <- data.frame()
  annotationb <- data.frame()
  prueba <- data.frame()
  database <- database
  dd <- data.frame()
  dff <- data.frame()
  xf_f_f<- data.frame()
  dff_f_f <- data.frame()
  dd_f_f <- data.frame()
  df2 <- data.frame()
  summary_ms2_f <- list()
  summary_ms2 <- list()
  summary2 <- list()
  summary <- list()
  prueba<-data.frame()

  b <- ms1[[8]]
  a <- ms2[[8]]
  database<-as.data.frame(database)

  for (i in 1:length(database[, 2])) {

    options(digits=9)
    fullms <- as.numeric(database[i, 2])
    fragments <- as.numeric(database[i, 3:length(database[i,
                                                          ])])
    fragments <-  fragments[!is.na(fragments)]
    name <- database[i, 1]
    options(digits=9)
    cc <- as.numeric(database[i, 2:length(database[i, ])])
    cc <-    cc[!is.na(cc)]

    ms <- subset(b, b[, 1] > fullms - f & b[, 1] < fullms + f)
    if (nrow(ms) == 0){next
      print(paste("No MS1 found for this metabolite, remove it from database or change algorithms parameters:",
                  name))
    }



    ms_2 <- cbind(ms, fullms)
    ms_2 <-cbind(ms_2,1)


    colnames(ms_2) <- c("m/z", "var_m/z", "max_int", "sum_int",
                        "RT", "minRT", "maxRT", "part_ID", "EIC_ID", "peak_ID",
                        "Score","fragments","MS")


    x <- data.frame()
    aiflist <- data.frame()
    for (i in 1:length(fragments)) {
      x<- subset(a, a[, 1] > fragments[i] - f & a[, 1] < fragments[i] + f)


      if (class(x) == "numeric" & nrow(x) != 0) {
        x <- as.data.frame(t(x))

        colnames(x) <- c("m/z", "var_m/z", "max_int",
                         "sum_int", "RT", "minRT", "maxRT", "part_ID",
                         "EIC_ID", "peak_ID","Score")
        aiflist <- rbind(aiflist, x)
        aiflist <- cbind(aiflist, fragments[i])
      }

      if (class(x) == "matrix" & nrow(x) != 0) {


        colnames(x) <- c("m/z", "var_m/z", "max_int",
                         "sum_int", "RT", "minRT", "maxRT", "part_ID",
                         "EIC_ID", "peak_ID","Score")


        x <- cbind(x, fragments=fragments[i])


        aiflist <- rbind(aiflist, x)

      }else {print(paste("Some MS2 ion not found for this metabolite, remove it from database or change algorithms parameters:",
                         name,fragments[i]))

        next

      }
      # next
    }

    #next
    if (nrow(aiflist) >0){
      xx <- data.frame()
      xx <- cbind(aiflist, 2)

      colnames(xx) <- c("m/z", "var_m/z", "max_int", "sum_int",
                        "RT", "minRT", "maxRT", "part_ID", "EIC_ID", "peak_ID","Score",
                        "fragments", "MS")
      }


    xy <- data.frame()

    xy <- rbind(ms_2, xx)
    xxf<- data.frame()
    xxf<-data.frame(xy)
    #xx <- data.frame(xy[, 1], xy[, 5], xy$minRT, xy$maxRT, xy$max_int, xy$sum_int, xy[, 5]/60, xy$MS)
    trrr <- xxf
    relative <- data.frame()
    X <- data.frame()
    ions <- data.frame()
    dd <- data.frame()
    dff <- data.frame()
    dF<- data.frame()
    ##aqui
    if (length(ms[, 1])>0){
      for (i in 1:length(ms[, 1])) {


        ions <- trrr[which(trrr[i, 5] - trrr[, 5] > -rtw &
                             trrr[i, 5] - trrr[, 5] < rtw), ]
        relative <- (ions$max_int)/(max(ions$max_int)) *
          100


        X <- data.frame(ions[, 1], ions[, 5]/60, ions[, 7]/60 - ions[, 6]/60, (ions[, 7] - ions[, 5])/(ions[,5] - ions$minRT), ions$sum_int, ions$max_int,
                        relative, ions$MS,ions$part_ID,ions$EIC_ID,ions$peak_ID,ions$fragments)
        is.na(X) <- sapply(X, is.infinite)
        colnames(X) <- c("Ion Metabolite (m/z)", "Rt(min)",
                         "Peak width (min)", "Peak Assymetry (B/A)", "Area",
                         "Signal Intensity", "Relative Intensity (%)",
                         "MS Level","PartID","EICID","PeakID","Fragments")


        ions <- data.frame()
        x <- data.frame()
        ions1bb <- data.frame()
        Xc <- data.frame()
        Xc <- subset(X, X$MS == "2")
        if (nrow(Xc) > 0) {
          for (i in 1:length(Xc[, 1])) {
            ions <- Xc[which(Xc[i, 1] - Xc[, 1] > -1 &
                               Xc[i, 1] - Xc[, 1] < 1), ]
            ions1bb <- data.frame(paste("", name), mean(ions$`Ion Metabolite (m/z)`),
                                  mean(ions$`Rt(min)`), mean(ions$`Peak width (min)`),
                                  sum(ions$Area), mean(ions$`Peak Assymetry (B/A)`), sum(ions$`Signal Intensity`),
                                  sum(ions$`Relative Intensity (%)`), ions$`MS Level`,ions$PartID, ions$EICID, ions$PeakID, as.numeric(nrow(Xc)),ions$Fragments)
            colnames(ions1bb) <- c("name", "Ion Metabolite (m/z)",
                                   "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity",
                                   "Relative intensity", "MS Level", "PartID","EICID","PeakID","Nions","Fragments")


            x <- rbind(ions1bb, x)
          }
          colnames(x) <- c("name", "Ion Metabolite (m/z)",
                           "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity",
                           "Relative intensity", "MS Level", "PartID","EICID","PeakID","Nions","Fragments")

          ionst <- aggregate(. ~ `Ion Metabolite (m/z)`,
                             data = x, FUN = mean,na.rm=TRUE, na.action=NULL)


          MS1xx <- subset(X, X$MS == "1")

          if (nrow(MS1xx)> 0 & nrow(ionst)> 0 ){

            MS1_f <- data.frame(MS1xx$`Ion Metabolite (m/z)`,
                                MS1xx$`Rt(min)`,MS1xx$`Peak width (min)`,
                                MS1xx$Area, MS1xx$`Peak Assymetry (B/A)`, MS1xx$`Signal Intensity`, MS1xx$`Relative Intensity (%)`, MS1xx$`MS Level`, MS1xx$PartID, MS1xx$EICID, MS1xx$PeakID,as.numeric(nrow(MS1xx)),MS1xx$Fragments)

            colnames(MS1_f) <- c("Ion Metabolite (m/z)",
                                 "Rt(min)", "Peak width (min)", "Area", "Peak Assymetry (B/A)","Signal Intensity", "Relative intensity",
                                 "MS Level","PartID","EICID","PeakID","Nions","Fragments")


            MS2_f <- data.frame(ionst$`Ion Metabolite (m/z)`,
                                ionst$`Rt(min)`, ionst$`Peak width (min)`,
                                ionst$Area,ionst$`Peak Assymetry (B/A)`, ionst$`Signal Intensity`, ionst$`Relative intensity`, ionst$`MS Level`, ionst$PartID, ionst$EICID,ionst$PeakID, ionst$Nions,ionst$Fragments)

            colnames(MS2_f) <- c("Ion Metabolite (m/z)",
                                 "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity", "Relative intensity",
                                 "MS Level","PartID","EICID","PeakID","Nions","Fragments")



            MS1MS2 <- rbind(MS1_f, MS2_f)
            xf <- data.frame(paste("", name), MS1MS2)
            dd <- data.frame(xf)
            colnames(dd) <- c("name","Ion Metabolite (m/z)",
                              "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity", "Relative intensity",
                              "MS Level","PartID","EICID","PeakID","Nions","Fragments")

            dff <- rbind(dd, dff)
          }

          if (nrow(MS1xx)== 0 & nrow(ionst)> 0 ){

            MS2_f <- data.frame(ionst$`Ion Metabolite (m/z)`,
                                ionst$`Rt(min)`, ionst$ `Peak width (min)`,
                                ionst$Area,ionst$`Peak Assymetry (B/A)`, ionst$`Signal Intensity`, ionst$`Relative intensity`,ionst$`MS Level`, ionst$PartID, ionst$EICID,ionst$PeakID, ionst$Nions,ionst$Fragments)

            colnames(MS2_f) <- c("Ion Metabolite (m/z)",
                                 "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity", "Relative intensity",
                                 "MS Level","PartID","EICID","PeakID","Nions","Fragments")

            xf <- data.frame(paste("", name), MS2_f)
            dd <- data.frame(xf)
            colnames(dd) <- c("name","Ion Metabolite (m/z)",
                              "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity", "Relative intensity",
                              "MS Level","PartID","EICID","PeakID","Nions","Fragments")

            dff <- rbind(dd, dff)

          }
          if (nrow(MS1xx)> 0 & nrow(ionst)== 0 ){

            MS1_f <- data.frame(MS1xx$`Ion Metabolite (m/z)`,
                                MS1xx$`Rt(min)`, MS1xx$`Peak width (min)`,
                                MS1xx$Area, MS1xx$`Peak Assymetry (B/A)`, MS1xx$`Signal Intensity`, MS1xx$`Relative Intensity (%)`, MS1xx$`MS Level`, MS1xx$PartID, MS1xx$EICID, MS1xx$PeakID,as.numeric(nrow(MS1xx)),MS1xx$Fragments)

            colnames(MS1_f) <- c("Ion Metabolite (m/z)",
                                 "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity", "Relative intensity",
                                 "MS Level","PartID","EICID","PeakID","Nions","Fragments")

            xf <- data.frame(paste("", name), MS1_f)
            dd <- data.frame(xf)
            colnames(dd) <- c("name","Ion Metabolite (m/z)",
                              "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity", "Relative intensity",
                              "MS Level","PartID","EICID","PeakID","Nions","Fragments")

            dff <- rbind(dd, dff)

          }

        } else {

          print(paste("No MS2 aligned for some MS1 for this metabolite, remove it from database or change algorithms parameters:", name))

          MS1xx <- subset(X, X$MS == "1")
          if (nrow(MS1xx) > 0) {
            for (i in 1:length(MS1xx[, 1])) {
              ions <- MS1xx[which(MS1xx[i, 1] - MS1xx[, 1] > -1 &
                                    MS1xx[i, 1] - MS1xx[, 1] < 1), ]
              ions1bb <- data.frame(paste("", name), mean(ions$`Ion Metabolite (m/z)`),
                                    mean(ions$`Rt(min)`), mean(ions$`Peak width (min)`),
                                    sum(ions$Area), mean(ions$`Peak Assymetry (B/A)`), sum(ions$`Signal Intensity`),
                                    sum(ions$`Relative Intensity (%)`), ions$`MS Level`,ions$PartID, ions$EICID, ions$PeakID, as.numeric(nrow(MS1xx)),ions$Fragments )
              x <- rbind(ions1bb, x)
            }

            colnames(x) <- c("name", "Ion Metabolite (m/z)",
                             "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity",
                             "Relative intensity", "MS Level", "PartID","EICID","PeakID","Nions","Fragments")

            dff<-x
          }

        }

        if (nrow(dff) == 0) {(return(print(paste("Remove this metabolite from database or change algorithms parameters:",
                                                 name))))}




        colnames(dff) <- c("name", "Ion Metabolite (m/z)",
                           "Rt(min)", "Peak width (min)", "Area","Peak Assymetry (B/A)", "Signal Intensity","Relative intensity",
                           "MS Level","PartID","EICID","PeakID","Nions","Fragments")


        dF<-rbind(dF,dff)


      }

      dF<-dF[!duplicated(dF[,12]),]
      ##aqui
      }else{

        NOms1<-data.frame()
        NOms1<- data.frame(paste("", name), trrr$m.z, trrr$max_int,trrr$sum_int,trrr$RT,
                           trrr$part_ID,trrr$EIC_ID,trrr$peak_ID,trrr$MS

        )


      }




    annotation <- dF


   # prueba <- rbind(annotation, prueba)


    if(nrow( annotation)>0){
      ms1anot <- annotation[which(annotation$`MS Level` == 1), ]


      metab <- ms1anot[which.max(ms1anot[, 5]), ]
      rtmetab <- metab[, 3]
      annotationb <- annotation[which(annotation[, 3] > rtmetab -
                                        rtw/60 & annotation[, 3] < rtmetab + rtw/60), ]
      MS1xxff <- subset(annotationb, annotationb$`MS Level` == "1")


      MS2xxff <- subset(annotationb, annotationb$`MS Level` == "2")


      if (nrow(MS2xxff) > 0 & nrow(MS1xxff) > 0) {
        ionst_ms1_xxff <- aggregate(. ~ `Ion Metabolite (m/z)`,
                                    data = MS1xxff, FUN = mean,na.rm=TRUE, na.action=NULL)
        ionst_ms2_xxff <- aggregate(. ~ `Ion Metabolite (m/z)`,
                                    data = MS2xxff, FUN = mean,na.rm=TRUE, na.action=NULL)

        annotation_final <- rbind(ionst_ms1_xxff, ionst_ms2_xxff)

      }



      if (nrow(MS2xxff) > 0 & nrow(MS1xxff) == 0) {

        ionst_ms2_xxff <- aggregate(. ~ `Ion Metabolite (m/z)`,
                                    data = MS2xxff, FUN = mean,na.rm=TRUE, na.action=NULL)


        annotation_final <- rbind(ionst_ms2_xxff)
      }

      if (nrow(MS2xxff) == 0 & nrow(MS1xxff) > 0) {

        ionst_ms1_xxff <- aggregate(. ~ `Ion Metabolite (m/z)`,
                                    data = MS1xxff, FUN = mean,na.rm=TRUE, na.action=NULL)

        annotation_final <-   rbind(ionst_ms1_xxff)

      }



      xf_f_f <- data.frame(paste("", name), annotation_final)
      xf_f_f <- xf_f_f[, -3]
      dd_f_f <- data.frame()
      dd_f_f <- data.frame(xf_f_f)

      names(dd_f_f)<-names(dff_f_f)
      dff_f_f <- rbind(dd_f_f, dff_f_f)
      jk_ms1 <- ionst_ms1_xxff$`Ion Metabolite (m/z)`
      jk_ms2 <- ionst_ms2_xxff$`Ion Metabolite (m/z)`
      app_F <- append(jk_ms1, jk_ms2)

      if (nrow(ionst_ms2_xxff) >0 ){
        jkr_ms2 <- round(jk_ms2, digits = 2)}


      jkr <- round(app_F, digits = 2)
      ccr <- round(cc, digits = 2)
      ccr2<- round(cc[-1], digits=2)
      dif_ms2 <- (setdiff(ccr2, jkr_ms2))
      dif <- (setdiff(ccr, jkr))
      summary2 <- list(NoDetected = c(dif))
      summary_ms2 <- list(NoDetected = c(dif_ms2))
      summary_ms2_f[[paste("", name)]] <- summary_ms2
      summary[[paste("", name)]] <- summary2





      colnames(dff_f_f) <- c("Metabolite", "Ion Metabolite (m/z)",
                             "Rt(min)", "Peak width (min)", "Area", "Peak Assymetry", "Signal Intensity", "Relative intensity",
                             "MS Level", "PartID","EICID","PeakID","Nions","Fragments")




      }

     prueba <- rbind(annotation, prueba)


    }

  dff_f_f$CE[dff_f_f$`MS Level`==1]<- 0
  dff_f_f$CE[dff_f_f$`MS Level`==2]<- CE
  results <- list(ms1 = ms1, ms2 = ms2, annotation = dff_f_f,
                  annotationfull = prueba, nodetected = summary, nodetectedMS2 = summary_ms2_f,RawData1=ms1$Scans[[2]], RawData2=ms2$Scans[[2]])
  write.csv2(dff_f_f,paste0("DIAannotated",CE,".csv"))
  write.csv2(prueba,paste0("Full_DIA",CE,".csv"))

  return(results)
}
