library(data.table)
library(ggplot2)

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
w <- 15
t <- 37
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

getDNAWindow <- function(text, r) {
  substr(text, r, r+ (winsize-1))
}


# DNA.df : Data frame with columns: name, sequence
calculateNNdG <- function(DNA.df, winsize, temp) {
  
  all.seq.params <- parallel::mclapply(X = 1:nrow(DNA.df), mc.cores = 3, FUN = function(i) {
    DNAseq <- DNA.df[i, "sequence"]
    name <- DNA.df[i, "name"]
    nDNA <- nchar(DNAseq)
    print(paste(name, DNAseq))
    prom.params <- lapply(1:(nDNA-(199+(winsize-1))), function (k) {
      DNAseqs50 <- sapply(k:(k+49), getDNAWindow, text = DNAseq)
      deltaG1s <- sapply(DNAseqs50, duplex_deltaG, t = temp)
      E1 <- mean(deltaG1s)
      
      DNAseqs100 <- sapply((k+100):(k+199), getDNAWindow, text = DNAseq)
      deltaG2s <- sapply(DNAseqs100, duplex_deltaG, t = temp)
      E2 <- mean(deltaG2s)
      return(list(bp =-400+(k-1), D = E1 - E2, E1 = E1))
    })
    prom.params <- as.data.frame(rbindlist(prom.params))
    prom.params["name"] <- name
    return(prom.params)
  })
  all.seq.params <- as.data.frame(rbindlist(all.seq.params))
}

promotores <- fread("K12_400_50_sites", sep = "\\", col.names = c("name", "sequence"), 
                    colClasses = c("char", "char", "NULL"))

nndgs <- calculateNNdG(promotores, w, t)
fwrite(nndgs, file = "NNdG_results.tsv", sep = "\t", col.names = T, row.names = F)

#ggplot(prom.params, aes(x = bp, y = D)) + geom_line() + geom_point()
#  scale_y_continuous(breaks=seq(-2,2,0.5))


