
AIF<-function(fileMS1,fileMS2,CE=0, database,rtw=7,
              ppm_tol=10,dmzgap=50,drtdens =20,drtgap=25,
              drtsmallMS1=100,drtsmallMS2=30,dmzdensMS1=15,
              dmzdensMS2=30,drtfill=5,drttotal=100,minpeakMS1=5,
              minpeakMS2=3,recurs=2,weight=2,SB=3, SN=2,
              minintMS1=1000,minintMS2=100,maxint=9e+09,
              ion_mode="positive",ppm=TRUE,ended=6){


  database<-as.data.frame(database)



  ms1 <- enviPickwrap(fileMS1, MSlevel = c(1), dmzgap, dmzdens=dmzdensMS1,
                      ppm=TRUE, drtgap, drtsmall=drtsmallMS1, drtdens, drtfill, drttotal, minpeak=minpeakMS1,
                      recurs, weight, SB, SN, minint=minintMS1, maxint, ion_mode = ion_mode,ended)

  ms2 <- enviPickwrap(fileMS2, MSlevel = c(2), dmzgap, dmzdens=dmzdensMS2,
                      ppm=TRUE, drtgap, drtsmall=drtsmallMS2, drtdens, drtfill, drttotal, minpeak=minpeakMS2,
                      recurs, weight, SB, SN, minint=minintMS2, maxint, ion_mode = ion_mode,ended)

  b <- ms1[[8]]
  a <- ms2[[8]]

  xxy<-data.frame()

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

    sel <- b[which(abs((b[,1]-fullms)/fullms)*1e+06 <= ppm_tol),]
    if (class(sel) == "numeric") {
      sel <- as.data.frame(t(sel))
    }
    if(nrow(sel)!=0){
    ppm_error<-(abs((sel[,1]-fullms)/fullms))*1e+06
    ms_1_list <- cbind(sel, ppm_error)
    ms_1_list_def <-cbind(ms_1_list,1)

    colnames(ms_1_list_def) <- c("m/z", "var_m/z", "max_int", "sum_int",
                        "RT", "minRT", "maxRT", "part_ID", "EIC_ID", "peak_ID",
                        "Score","ppm_error","MS")

    x <- data.frame()
    aiflist <- data.frame()
    ms_2_list <-data.frame()
    ms_2_list_def <-data.frame()
    xy <- data.frame()
    for (i in 1:length(fragments)) {
      x <- a[which(abs((a[,1]-fragments[i])/fragments[i])*1e+06 <= ppm_tol),]

      if (class(x) == "numeric") {
        x <- as.data.frame(t(x))

        ppm_error<-(abs((x[,1]-fragments[i])/fragments[i]))*1e+06
        ms_2_list <- cbind(x, ppm_error)
        ms_2_list_def <-cbind(ms_2_list,2)
        colnames(ms_2_list_def) <- c("m/z", "var_m/z", "max_int", "sum_int",
                                     "RT", "minRT", "maxRT", "part_ID", "EIC_ID", "peak_ID",
                                     "Score","ppm_error","MS")
        aiflist <- rbind(aiflist, ms_2_list_def)
      }

      if (class(x) == "matrix" & nrow(x) != 0) {

        ppm_error<-(abs((x[,1]-fragments[i])/fragments[i]))*1e+06
        ms_2_list <- cbind(x, ppm_error)
        ms_2_list_def <-cbind(ms_2_list,2)
        colnames(ms_2_list_def) <- c("m/z", "var_m/z", "max_int", "sum_int",
                                     "RT", "minRT", "maxRT", "part_ID", "EIC_ID", "peak_ID",
                                     "Score","ppm_error","MS")
        aiflist <- rbind(aiflist, ms_2_list_def)

      }else
        next

      }


    xy <- rbind(aiflist, ms_1_list_def)

    B<-(xy[,7]-xy[,5])
    A<-(xy[,5]-xy[,6])

    for (i in seq_along(B)) {
      if (B[i]==0){
        B[i]<-0.001
      }}
    for (i in seq_along(A)) {
      if (A[i]==0){
        A[i]<-0.001
      }
      }
    xy <- data.frame(paste("", name), xy[,1],xy[,5]/60,xy[,7]/60-xy[,6]/60,xy[,4],xy[,3],xy[,10],xy[,12],xy[,13],B/A

                     )


    colnames(xy) <- c("name", "Ion Metabolite (m/z)","Rt(min)", "Peak width (min)", "Area", "Signal Intensity",
                      "PeakID","Error (ppm)","MS Level","Peak Assymetry")

    xxy <- rbind(xxy, xy)

     }


  }

  write.csv2(xxy,paste0("Full_DIA",CE,".csv"))
  results <- list(ms1 = ms1, ms2 = ms2,
                  annotationfull =  xxy,RawData1=ms1$Scans[[2]], RawData2=ms2$Scans[[2]])

  return(  results )
  }



