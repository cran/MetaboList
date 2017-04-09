annotation<- function (file, database, rtw = 3, mzw = 0.004, dmzgap = 20,
                 dmzdens = 8, drtgap = 1000, drtsmall = 5, drtdens = 5, drtfill = 10,
                 drttotal = 200, minpeak = 2, recurs = 200, weight = 0.25,
                 SB = 1.5, SN = 0.5, minint = 10, maxint = 9e+09, ion_mode = "positive",
                 ppm = TRUE)
{
  if (is.numeric(rtw) == FALSE) {
    stop("Please, establish a numeric retention time window ")
  }
  if (is.numeric(mzw) == FALSE) {
    stop("Please, establish a numeric retention time window ")
  }

  ms1 <- enviPickwrap(file, MSlevel = c(1), dmzgap, dmzdens,
                      ppm, drtgap, drtsmall, drtdens, drtfill, drttotal, minpeak,
                      recurs, weight, SB, SN, minint, maxint, ion_mode = ion_mode)
  ms2 <- enviPickwrap(file, MSlevel = c(2), dmzgap, dmzdens,
                      ppm, drtgap, drtsmall, drtdens, drtfill, drttotal, minpeak,
                      recurs, weight, SB, SN, minint, maxint, ion_mode = ion_mode)

  database <- database
  dd<-data.frame()
  dff <- data.frame()
  dff_f_f<- data.frame()
  dd_f_f<-data.frame()
  df2 <- data.frame()
  summary_ms2_f<-list()
  summary_ms2<-list()
  summary2 <- list()
  summary <- list()
  b <- ms1[[8]]
  a <- ms2[[8]]
  f <- mzw
  rtw<-rtw


  for (i in 1:length(database[,2])){

    fullms <- database[i, 2]
    fragments <- as.numeric(database[i, 3:length(database[i,])])
    name <- database[i, 1]
    cc <- as.numeric(database[i, 2:length(database[i,])])
    ms<-subset(b, b[,1]> fullms - f& b[, 1] < fullms + f )
    ms_2<-cbind( ms,1)
    colnames(  ms_2) <- c("m/z", "var_m/z", "max_int",
                          "sum_int", "RT", "minRT", "maxRT", "part_ID",
                          "EIC_ID", "peak_ID", "Score","MS")



    x <- data.frame()
    aiflist <- data.frame()
    for (i in 1:length(fragments)) {

      x <- a[which(a[, 1] > fragments[i] - f & a[,1] < fragments[i] + f), ]

      if(class(x)=="numeric"){

        x<-as.data.frame(t(x))


        colnames(x) <- c("m/z", "var_m/z", "max_int",
                         "sum_int", "RT", "minRT", "maxRT", "part_ID",
                         "EIC_ID", "peak_ID", "Score")
        aiflist <- rbind(aiflist, x)
      }

      if(class(x)=="matrix"){


        colnames( x) <- c("m/z", "var_m/z", "max_int",
                          "sum_int", "RT", "minRT", "maxRT", "part_ID",
                          "EIC_ID", "peak_ID", "Score")
        aiflist <- rbind(aiflist, x)
      }}


    xx<-cbind(aiflist,2)
    colnames( xx) <- c("m/z", "var_m/z", "max_int",
                       "sum_int", "RT", "minRT", "maxRT", "part_ID",
                       "EIC_ID", "peak_ID", "Score","MS")

    xy <- data.frame()
    xy <- rbind( ms_2, xx)


    xx <- data.frame(xy[, 1], xy[, 5], xy$minRT, xy$maxRT, xy$max_int, xy$sum_int, xy[, 5]/60,xy$MS)
    trrr <- xx

    relative <- data.frame()
    X <- data.frame()
    ions <- data.frame()
    dd<-data.frame()
    dff<-data.frame()

    for (i in 1:length(ms[, 1])) {
      ions <- trrr[which(trrr[i, 2] - trrr[, 2] > -rtw & trrr[i, 2] - trrr[, 2] < rtw), ]
      relative <- (ions$xy.max_int)/(max(ions$xy.max_int)) * 100

      X <- data.frame(ions[, 1], ions[, 2]/60, ions[, 4]/60 - ions$xy.minRT/60, (ions[, 4] - ions[,2])/(ions[, 2] - ions$xy.minRT), ions$xy.sum_int,
                      ions$xy.max_int, relative,ions$xy.MS)
      colnames(X) <- c("Ion Metabolite (m/z)", "Rt(min)",
                       "Interval Peak (min)", "Assymetry (B/A)", "Area",
                       "Maximum Intensity", "Relative Intensity (%)","MS")



      ions <- data.frame()
      x <- data.frame()
      ions1bb <- data.frame()
      Xc<-data.frame()

      Xc<-subset(X,X$MS=="2")

      if (nrow(Xc) >0) {


        for (i in 1:length(Xc[, 1])) {
          ions <- Xc[which(Xc[i, 1] - Xc[, 1] > -1 & Xc[i, 1] - Xc[, 1] < 1), ]


          ions1bb <- data.frame(paste("", name), mean(ions$`Ion Metabolite (m/z)`),
                                mean(ions$`Rt(min)`), mean(ions$`Interval Peak (min)`),
                                sum(ions$Area), sum(ions$`Maximum Intensity`),sum(ions$`Relative Intensity (%)`),ions$MS)
          x <- rbind(ions1bb, x)
        }



        colnames(x) <- c("name", "Ion Metabolite (m/z)",
                         "Rt(min)", "Interval Peak (min)", "Area", "Maximum Intensity","Relative intensity","MS")

        ionst <- aggregate(. ~ `Ion Metabolite (m/z)`,data = x, FUN = mean)
        MS1xx<-subset(X,X$MS=="1")


        MS1_f<- data.frame(MS1xx$`Ion Metabolite (m/z)`,MS1xx$`Rt(min)`,MS1xx$`Interval Peak (min)`,MS1xx$Area,MS1xx$`Maximum Intensity`,MS1xx$MS)
        MS2_f<- data.frame(ionst$`Ion Metabolite (m/z)`,ionst$`Rt(min)`,ionst$`Interval Peak (min)`,ionst$Area,ionst$`Maximum Intensity`,ionst$MS)
        colnames(MS1_f) <- c("Ion Metabolite (m/z)", "Rt(min)",
                             "Interval Peak (min)", "Area",
                             "Maximum Intensity","MS")
        colnames(MS2_f) <- c("Ion Metabolite (m/z)", "Rt(min)",
                             "Interval Peak (min)", "Area",
                             "Maximum Intensity","MS")


        MS1MS2<-rbind(MS1_f,MS2_f)

        xf <- data.frame(paste("", name), MS1MS2)


        dd <- data.frame(xf)
        dff <- rbind(dd, dff)

        if (nrow(dff)==0){print("No metabolites were annotated")}
      }

    }

    colnames(dff) <- c("Metabolite", "Ion Metabolite (m/z)",
                       "Rt(min)", "Interval Peak (min)", "Area", "Maximum Intensity","MS")
    annotation <- dff

    ann<-metab<-annotation[which.max(annotation[,6]),]
    rtmetab<-metab[,3]

    annotationb<- annotation[which(annotation[, 3]  >  rtmetab  -rtw/60 & annotation[, 3] <   rtmetab  +rtw/60), ]

    MS1xxff<-subset( annotationb, annotationb$MS=="1")
    MS2xxff<-subset( annotationb, annotationb$MS=="2")

    #  if (nrow(MS2xxff)==0){ print(paste("",name),

    #                               "No MS2 fragments found in this retention time, please change retention time windows or remove metabolite from database")}

    ionst_ms1_xxff <- aggregate(. ~ `Ion Metabolite (m/z)`,data = MS1xxff, FUN = mean)
    ionst_ms2_xxff <- aggregate(. ~ `Ion Metabolite (m/z)`,data = MS2xxff, FUN = mean)



    annotation_final<-rbind(ionst_ms1_xxff,ionst_ms2_xxff)
    xf_f_f <- data.frame(paste("", name), annotation_final)
    xf_f_f<-xf_f_f[,-3]
    dd_f_f <- data.frame(xf_f_f)
    dff_f_f <- rbind(dd_f_f, dff_f_f)

    jk_ms1 <-  ionst_ms1_xxff$`Ion Metabolite (m/z)`
    jk_ms2 <-  ionst_ms2_xxff$`Ion Metabolite (m/z)`
    app_F<-append(  jk_ms1 ,  jk_ms2 )

    jkr_ms2<-round(   jk_ms2,digits=2)
    jkr <- round( app_F, digits = 2)
    ccr <- round(cc, digits = 2)
    dif_ms2<-(setdiff(ccr,jkr_ms2))
    dif <- (setdiff(ccr, jkr))

    summary2 <- list(NoDetected = c(dif))
    summary_ms2<-list(NoDetected=c(dif_ms2))
    summary_ms2_f[[paste("", name)]] <- summary_ms2
    summary[[paste("", name)]] <- summary2

  }

  colnames(dff_f_f) <- c("Metabolite", "Ion Metabolite (m/z)",
                         "Rt(min)", "Interval Peak (min)", "Area", "Maximum Intensity","MS")
  results <- list(ms1 = ms1[[8]], ms2 = ms2[[8]], annotation=dff_f_f,annotationfull = annotation,
                  nodetected = summary,nodetectedMS2 = summary_ms2_f)
  return(results)
}
