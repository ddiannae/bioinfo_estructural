setwd("/home/diana/Workspace/PhD/semestre2/bioinfo_estructural/")

adiedros <- read.delim("1ayi-adiedros.tsv", header = F, 
                       sep = "\t", col.names = c("residuos", "amino", "chain", "phi", "psi", "omega"),
                       colClasses = c("integer", "character", "character", "numeric", "numeric", "numeric"))
phipsi <- c("phi", "psi")
beta <- matrix(c(-133, 171, -133, 92, -55, 92, -55, 171), nrow = 4, byrow = T, dimnames = list(NULL, phipsi))
alpha1 <- matrix(c(-53, -40, -53, -66, -159, -66, -159, -40), nrow = 4, byrow = T,  dimnames = list(NULL, phipsi))
alpha2 <- matrix(c(64, 96, 49, 96, 49, 20, 64, 20), nrow = 4, byrow = T, dimnames = list(NULL, phipsi))

isInSquare <- function(p, cordref) {
  return(p[, "phi"] <= max(cordref[, "phi"]) & p[, "phi"] >= min(cordref[, "phi"]) & 
           p[, "psi"] <= max(cordref[, "psi"]) & p[, "psi"] >= min(cordref[, "psi"]) )
}

setHEC <- function(p){
  print(p)
  coordref <- list(beta, alpha1, alpha2)
  asignacion <- sapply(coordref, isInSquare, p = p)
  print(asignacion)
  if(asignacion[1]) {
      return("E") 
  } else if(asignacion[2] | asignacion[3]) {
      return("H")
  } else {
    return("C")
  }
}

HECs <- lapply(1:nrow(adiedros), function(i){
  asig <- setHEC(adiedros[i, ])
  asig
})

asignaciones <- do.call(rbind, HECs)
adiedros[, "asig"] <- asignaciones
write.table(adiedros, "1ayi-adiedros-asignaciones.tsv", col.names = T, row.names = F, quote = F)



