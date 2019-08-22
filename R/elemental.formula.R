elemental.formula<-function (elemental.formula){
  chr <- gregexpr("[[:upper:]][[:lower:]]{0,1}", elemental.formula)
for (i in 1:length(chr[[1]])) {
  y <- attr(chr[[1]], which = "match.length")[i]
  z <- substr(elemental.formula, chr[[1]][i], chr[[1]][i] + 
                y - 1)
  if (!(z == "C" | z == "H" | z == "N" | z == "O" | z == 
        "S" | z == "P" | z == "Br" | z == "Cl" | z == "F" | 
        z == "I" | z == "Si" | z == "Sn" | z == "Mg"| z == "As"| z == "D")) 
    warning(paste("Elemental formula", elemental.formula, 
                  "contains element not of C,H,N,O,S,P,Br,Cl,F,I,Si,Mg,As,D,Sn."))
}
GetAtoms <- function(elemental.formula, element) {
  reg.exp <- paste(element, "[[:digit:]]*(?![[:lower:]])", 
                   sep = "")
  x <- gregexpr(reg.exp, elemental.formula, perl = TRUE)
  if (x[[1]][1] != -1) {
    n <- vector(mode = "numeric", length = length(x[[1]]))
    for (i in 1:length(x[[1]])) {
      y <- attr(x[[1]], which = "match.length")[i]
      z <- substr(elemental.formula, x[[1]][i], x[[1]][i] + 
                    y - 1)
      number <- as.numeric(strsplit(z, split = element)[[1]][2])
      if (is.na(number)) {
        n[i] <- 1
      }
      else {
        n[i] <- number
      }
      atoms <- sum(n)
    }
  }
  else {
    atoms <- 0
  }
  return(atoms)
}
elements <- c("C", "H", "N", "O", "S", "P", "Br", "Cl", "F", 
              "I", "Si", "Sn","Mg","D")
result <- as.list(sapply(elements, function(x) {
  GetAtoms(elemental.formula, x)
}))
return(result)
}