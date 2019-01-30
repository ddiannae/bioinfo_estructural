library(data.table)


NNparams <- list(c("H" = -7.9, "S" = -22.2),
                 c("H" = -7.2, "S" = -20.4),
                 c("H" = -7.2, "S" = -21.3),
                 c("H" = -8.5, "S" = -22.7),
                 c("H" = -8.4, "S" = -22.4),
                 c("H" = -7.8, "S" = -21.0),
                 c("H" = -8.2, "S" = -22.2),
                 c("H" = -10.6,"S" = -27.2),
                 c("H" = -9.8, "S" = -24.4),
                 c("H" = -8.0, "S" = -19.9),
                 # initiation costs
                 c("H" = 0.1, "S" = -2.8),
                 c("H" = 2.3, "S" = 4.1),
                 # symmetry correction
                 c("H" = 0, "S" = -1.4) )
names(NNparams) = c('AA/TT','AT/TA','TA/AT','CA/GT','GT/CA','CT/GA','GA/CT','CG/GC','GC/CG','GG/CC','G','A','sym')
complement <- list("T", "A", "C", "G")
names(complement) <- c("A", "T", "G", "C")
Dc <- 3.3
E1c <- -17.55
winsize <- 15

DNA <- "TCGCGCCCTGCCAGCATTTCAACAGGAGGATGCCAGACTAAAGCATAATCAGCAGAGTCATTATCTCCGCTTTTCCATGCTCTGACTCTTGCCTGAGGAATAGCTTTGCGCAGTGCCTCAATCCACCATTGGGTATCGAACGTTGGGTGATAAAAGATGATATCCATACTGACTCCCGAAAAGCGTTGTGCGAAATTTATTCGCACTTATCGTTATGATCTACAAAGGCCACCAGCATAACAAATCCGTGGTCGGTGGCAAAAAAAGCAGATTTCGCTTATTAAAACCACACATTGATTGAAGTTTGAATAAACGCGCGATTTTTTCAAAAAGTTTGTTGACCTCAGGTCATGATTTCCCTAAATTAGCGCCCGTTCCAGCAAGACAGGAACGACAATTTGGTGAGGTGTCCGAGTGGCTGAAGGAGCACGCCTGGAAAGTGTGTATACGG"
nDNA <- nchar(DNA)

getDNAWindow <- function(text, r) {
  substr(text, r, r+ (winsize-1))
}

prom.params <- lapply(1:(nDNA-(199+(winsize-1))), function (k) {
  DNAseqs50 <- sapply(k:(k+49), getDNAWindow, text = DNA)
  deltaG1s <- sapply(DNAseqs50, duplex_deltaG, t = 37)
  E1 <- mean(deltaG1s)
  
  DNAseqs100 <- sapply((k+100):(k+199), getDNAWindow, text = DNA)
  deltaG2s <- sapply(DNAseqs100, duplex_deltaG, t = 37)
  E2 <- mean(deltaG2s)
  list(bp = k, D = E1 - E2, E1 = E1)
})

prom.params <- rbindlist(prom.params)
positive.signals <- prom.params[D > Dc, ]


DNAlist <- strsplit(toupper(DNA), "")[[1]]
nDNA <- length(DNA)

for(k in 1: (nDNA-213)) {

  print(paste("i: ", k, k+49))
  sumGE1 <- 0
  for(i in k:(k+49)) {
    print(i)
    sumGE1 <- sumGE1 + duplex_deltaG(paste(DNAlist[i:(i+14)], collapse =""), 37)
  }
  promG <- sumGE1 / 50
  
  sumGE2 <- 0
  print(paste("j: ", (k+99), (k+199)))
  for(j in (k + 99):(k+199)) {
    sumGE2 <- sumGE2 + duplex_deltaG(DNA[j:(j+14)], 37)
  }
  promG2 <- sumGE2 /100

  D <- promG - promG2
}


duplex_deltaG <- function (DNAseq, t) {
  total_dG = 0
  tK = 273.15 + t
  
  sequence <- unlist(strsplit(toupper(DNAseq), ""))
  nseq <- length(sequence)
  
  for(n in 1:(nseq-1)) {
    DNAstep <- paste(sequence[n], sequence[n+1], "/", 
                     complement[sequence[n]], complement[sequence[n+1]], sep="")
    if(is.null(NNparams[[DNAstep]])) {
        DNAstep <- paste(rev(strsplit(DNAstep, NULL)[[1]]), collapse="" )
    }
    dg <- (1000 * NNparams[[DNAstep]]["H"] - tK * NNparams[[DNAstep]]["S"])/1000
    names(dg) <- NULL
    total_dG = total_dG + dg
  }
  
  nt <- sequence[1]
  if(is.null(NNparams[[nt]])) {
    nt <- complement[[nt]]
  }
  
  total_dG = total_dG + (1000 * NNparams[[nt]]["H"] - tK * NNparams[[nt]]["S"]) /1000
  
  nt <- sequence[nseq]
  if(is.null(NNparams[[nt]])) {
    nt <- complement[[nt]]
  }
  
  total_dG = total_dG + (1000 * NNparams[[nt]]["H"] - tK * NNparams[[nt]]["S"]) /1000
  
  mid <- round(nseq/2)
  
  sym <- TRUE
  for(i in 1:(mid-1)) {
    if(sequence[i] != complement[sequence[nseq - (i-1)]]) {
      sym <- FALSE
      break
    }
  }
  
  if(sym) {
    total_dG = total_dG + (1000 * NNparams[["sym"]]["H"] - tK * NNparams[["sym"]]["S"]) /1000
  }
  names(total_dG) <- NULL
  return(total_dG)  
}


duplex_deltaG("CGTTGA", 37)

  
