
Isotopic<-function(library,name=c("M+H"),adducts) {

  getpeaksM <- matrix(nrow = nrow(library), ncol = 15, 0)
  spectrumlist<-list()

  for (i in 1:nrow(library)){

  formulaB<-library[i,1]
  RT<-library[i,2]
  Mass<-library[i,3]
  CompoundName<-library[i,4]
  formula <- elemental.formula(formulaB)



  adductsR<-adducts[which(adducts$Name==name),]
  charge<-adductsR$Charge
  mult<-adductsR$Mult
  inputFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0,
                       Br = 0, Cl = 0, F = 0, Si = 0,Mg=0,D=0)
  inputFormula[names(formula)] <- formula

  simulation <- function(inputFormula) {
    massCarbon <- sum(sample(c(12, 13.0033548378), size = inputFormula$C,
                             replace = TRUE, prob = c(0.9893, 0.0107)))
    massHydrogen <- sum(sample(c(1.0078250321, 2.014101778),
                               size = inputFormula$H, replace = TRUE, prob = c(0.999885,
                                                                               0.000115)))

    massDeuterium<- inputFormula$D * 2.014102

    massNitrogen <- sum(sample(c(14.0030740052, 15.0001088984),
                               size = inputFormula$N, replace = TRUE, prob = c(0.99632,
                                                                               0.00368)))

    massMagnesium <- sum(sample(c(23.985041700, 24.98583692,25.982592929),
                                size = inputFormula$Mg, replace = TRUE, prob = c(0.7899,
                                                                                 0.10000,0.1101)))

    massOxygen <- sum(sample(c(15.9949146221, 16.9991315,
                               17.9991604), size = inputFormula$O, replace = TRUE,
                             prob = c(0.99757, 0.00038, 0.00205)))
    massSulfer <- sum(sample(c(31.97207069, 32.9714585, 33.96786683,
                               35.96708088), size = inputFormula$S, replace = TRUE,
                             prob = c(0.9493, 0.0076, 0.0429, 2e-04)))
    massPhosphorus <- inputFormula$P * 30.97376151
    massBromine <- sum(sample(c(78.9183376, 80.916291), size = inputFormula$Br,
                              replace = TRUE, prob = c(0.5069, 0.4931)))
    massChlorine <- sum(sample(c(34.96885271, 36.9659026),
                               size = inputFormula$Cl, replace = TRUE, prob = c(0.7578,
                                                                                0.2422)))
    massFluorine <- inputFormula$F * 18.9984032
    massArsenic  <- inputFormula$As * 74.9215965

    massSilicon <- sum(sample(c(27.9769265327, 28.97649472,
                                29.97377022), size = inputFormula$Si, replace = TRUE,
                              prob = c(0.922297, 0.046832, 0.030872)))


    massMolecule <- sum(massCarbon, massHydrogen, massNitrogen,
                        massOxygen, massSulfer, massPhosphorus, massBromine,
                        massChlorine, massFluorine, massSilicon,massMagnesium,massArsenic,massDeuterium)
    mz <- massMolecule/abs(charge)
    return(mz)
  }
  sim <- replicate(10000, expr = simulation(inputFormula))

  mz<-simulation(inputFormula)
  b <- seq(from = min(sim) - (1/(2 * abs(charge))), to = max(sim) +
             1, by = 1/abs(charge))
  bins <- cut(sim, breaks = b)
  massIso <- round(tapply(sim, bins, mean), digits = 6)
  intensity <- as.vector(table(bins))

  massIsoAd<- massIso+adductsR$Mass
  mz<- mult*mz+adductsR$Mass
  spectrum <- data.frame(massIsoAd, intensity)
  spectrum <- spectrum[spectrum$intensity != 0, ]
  spectrum$percent <- with(spectrum, round(intensity/max(intensity) *
                                             100, digits = 2))
  row.names(spectrum) <- 1:(nrow(spectrum))

  percent<-spectrum$percent


  from <- 5

  to <- (from + as.numeric(length(massIsoAd)) -1)

  getpeaksM[i,1] <- as.character(formulaB)
  getpeaksM[i,2] <- RT
  getpeaksM[i,3] <-  Mass
  getpeaksM[i,4] <-  as.character(CompoundName)
  getpeaksM[i,from:to] <- massIsoAd

  spectrumlist[[as.character(CompoundName)]]<-list(spectrum)

  b<-spectrumlist
  }

  getpeaksM<-as.data.frame(getpeaksM)

  colnames(getpeaksM)[1] <- c("Formula")
  colnames(getpeaksM)[2] <- c("RT")
  colnames(getpeaksM)[3] <- c("Neutral Mass")
  colnames(getpeaksM)[4] <- c("CompoundName")
  colnames(getpeaksM)[from:15]<-paste("",name)

  getpeaksM[getpeaksM == 0] <- NA

  getpeaksM<- getpeaksM[,colSums(is.na(getpeaksM))<nrow(getpeaksM)]

  write.csv(getpeaksM,paste0("",name,".csv"))
  return(list(getpeaksM=getpeaksM,b=b,spectrumlist=spectrumlist,monoisotopic=mz))


}

