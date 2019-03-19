
CE.isolation<-function (file, output)
{
  if (file.exists(file)) {
    lines <- readLines(file)
    MS1 <- grep("          msLevel=\"1\"", lines)
    if (length(MS1) > 0) {
      fullMS <- TRUE
    }else {fullMS <- FALSE}
    MS2 <- grep("          msLevel=\"2\"", lines)
    lCe <- grep("          collisionEnergy=.+\"", lines)
    if (length(MS2) > 0 || length(lCe) > 0) {
      if (length(lCe) > 0) {
        Ce <- grep("          collisionEnergy=.+\"",
                   lines, value = T)
        Ce <- unlist(strsplit(Ce, "          collisionEnergy="))
        Ce <- unlist(strsplit(Ce[Ce != ""], "\""))
        Ce <- as.numeric(Ce[Ce != ""])
        CE <- unique(Ce)
        #lines <- gsub("          msLevel=\"2\"", "          msLevel=\"1\"",
         #             lines)
        scans <- grep("<scan", lines)
        runs <- grep("</msRun>", lines)
        v <- vector()
        if (fullMS == TRUE) {
          v <- append(v, c(CE[1:length(which(MS2 < MS1[1]))],
                           "fullMS"))
          v <- append(v, CE[length(which(MS2 < MS1[1])) +
                              1:length(CE)])
          v <- v[which(v != "NA")]
        }else {v <- as.character(CE)}

        for (i in 1:length(v)) {
          pos <- seq(i, length(scans), length(v))
          lines2write <- c(1:(scans[1] - 1))
          for (x in pos[1:length(pos) - 1]) {
            if (x == scans[length(scans)]) {
              lines2write <- append(lines2write, c(scans[x]:c(runs -
                                                                1)))
            }
            else {
              lines2write <- append(lines2write, c(scans[x]:(scans[x +
                                                         1] - 1)))
            }
          }
          lines2write <- append(lines2write, c(runs:length(lines)))
          write(lines[lines2write], file = paste(c(output,
                                                   v[i], "pp.mzXML"), collapse = ""))
        }
      }

    }

  }

}
