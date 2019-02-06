library(data.table)
library(ggplot2)
library(seqinr)
library(dplyr)
library(reshape2)
setwd("/home/diana/Workspace/PhD/semestre2/bioinfo_estructural/")

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
names(NNparams) = c('AA/TT','AT/TA','TA/AT','CA/GT','GT/CA',
                    'CT/GA','GA/CT','CG/GC','GC/CG','GG/CC','G','A','sym')
complement <- list("T" = "A", "A" = "T", "C" = "G", "G" = "C")

getDNAWindow <- function(text, r, w) {
  substr(text, r, r+(w-1))
}

getRevComplement <- function(text) {
  comp <- lapply(strsplit(text, "")[[1]], function(c){
    return(complement[[c]])
  })
  return(paste(rev(comp), collapse="" ))
}

sumExtremos <- function(seq, tK) {
  nseq <- nchar(seq)
  exts <- c(substr(seq, 1, 1), substr(seq, nseq, nseq))
  dG_ext <- sapply(exts, function(nt){
    if(is.null(NNparams[[nt]])) {
      nt <- complement[[nt]]
    }
    return((1000 * NNparams[[nt]]["H"] - tK * NNparams[[nt]]["S"]) /1000)
  })
  return(sum(dG_ext))
}

sumSym <- function(seq, tK) {
  nseq <- nchar(seq)
  mid <- trunc(nseq/2)
  a <- substr(seq, 1, mid)
  b <- getRevComplement(substr(seq, nseq-mid + 1, nseq))
  
  if(a == b) {
    return((1000 * NNparams[["sym"]]["H"] - tK * NNparams[["sym"]]["S"]) /1000)
  }
  return(0)
}

sumStep <- function(step, tK) {
  if(is.null(NNparams[[step]])) {
    step <- paste(rev(strsplit(step, NULL)[[1]]), collapse="" )
  }
  return((1000 * NNparams[[step]]["H"] - tK * NNparams[[step]]["S"])/1000)
}

duplex_deltaG <- function (DNAseqs, t) {
  total_dG = c()
  tK = 273.15 + t
  DNAseqs <- sapply(DNAseqs, toupper)
  firstseq <- unlist(strsplit(DNAseqs[[1]], ""))
  nseq <- length(firstseq)
  
  firstdG <- sapply(1:(nseq-1), function(n) {
    DNAstep <- paste(firstseq[n], firstseq[n+1], "/", 
                     complement[firstseq[n]], complement[firstseq[n+1]], sep="")
    return(sumStep(DNAstep, tK))
  })
  
  total_dG <- append(total_dG, sum(firstdG))
  
  for(i in 2:length(DNAseqs)) {
    prevVal <- total_dG[(i-1)]
    prevfirst <- paste(substr(DNAseqs[i-1], 1, 2), "/", 
                       complement[substr(DNAseqs[i-1], 1, 1)], 
                       complement[substr(DNAseqs[i-1], 2, 2)], sep="")
    newlast <- paste(substr(DNAseqs[i], nseq-1 , nseq), "/", 
                     complement[substr(DNAseqs[i], nseq-1, nseq-1)], 
                     complement[substr(DNAseqs[i], nseq, nseq)], sep="")
    total_dG <- append(total_dG, prevVal - sumStep(prevfirst, tK) + sumStep(newlast, tK))
  }
  
  ext_sym <- sapply(DNAseqs, function(dna) {
    sumExtremos(dna, tK) + sumSym(dna, tK)
  })
  
  total_dG = total_dG + ext_sym
  names(total_dG) <- DNAseqs
  return(total_dG)  
}


# DNA.df : Data frame with columns: name, sequence
calculateNNdG <- function(DNA_df, winsize, temp) {
  
  all.seq.params <- parallel::mclapply(X = 1:nrow(DNA_df), mc.cores = 3, FUN = function(i) {
    DNAseq <- DNA_df[i, "sequence"]
    name <- DNA_df[i, "name"]
    nDNA <- nchar(DNAseq)
    print(paste(name, DNAseq))
    prom.params <- lapply(1:(nDNA-(199+(winsize-1))), function (k) {
       DNAseqs50 <- sapply(k:(k+49), getDNAWindow, text = DNAseq, w = winsize)
       deltaG1s <- duplex_deltaG(DNAseqs50, temp)
       E1 <- mean(deltaG1s)
       
       DNAseqs100 <- sapply((k+100):(k+199), getDNAWindow, text = DNAseq, w = winsize)
       deltaG2s <- duplex_deltaG(DNAseqs100, temp)
       E2 <- mean(deltaG2s)
       D = E1 - E2
       return(list(bp= k, D = E1 - E2, E1 = E1))
     })
    prom.params <- rbindlist(prom.params)
    prom.params[, "name" := name]
    return(prom.params)
  })
  all.seq.params <- rbindlist(all.seq.params)
}

promotores <- fread("K12_400_50_sites", sep = "\\", col.names = c("name", "sequence"), 
                    colClasses = c("char", "char", "NULL"))

nndgs <- calculateNNdG(promotores, 15, 37)

ggplot(nndgs, aes(x = bp - 400, y = D, color = name)) + geom_line()  +
  theme(legend.position="none") + xlab("bp")

mindgs <- nndgs %>% group_by(name) %>% summarize(minD = min(D), bp = bp[which.max(D)]) 

ggplot(mindgs, aes(x = minD)) + geom_histogram(bins = 25, fill = "tan1")  +
  theme(legend.position="none") + xlab("min D value") + theme_bw()
ggplot(mindgs, aes(x = bp-400)) + geom_histogram(bins = 25, fill = "tan1")  +
  theme(legend.position="none")  + xlab("bp") + theme_bw()

promotores <- read.fasta("promotor_sequences.faa")
promotores <- promotores[1:100]
lista.promotores <- lapply(promotores, paste, collapse = "")
lista.promotores <- lapply(lista.promotores, toupper)
promotores <- data.table(name = unlist(lapply(strsplit(names(promotores), split ="\t"), first)), sequence = unlist(lista.promotores))

nndgs <- calculateNNdG(promotores, 15, 37)

ggplot(nndgs, aes(x = bp - 500, y = D, color = name)) + geom_line()  +
  theme(legend.position="none") + xlab("bp")
ggplot(nndgs, aes(x = bp - 500, y = E1, color = name)) + geom_line()  +
  theme(legend.position="none") + xlab("bp")

E1cutoffs <- c(seq(-15, -20, -0.25))
Dcutoffs <- c(seq(0.5, 3.5, 0.25))
sensit_matrix <- matrix(rep(0, length(E1cutoffs)*length(Dcutoffs)),
                     nrow = length(E1cutoffs), byrow = T)
precis_matrix <- matrix(rep(0, length(E1cutoffs)*length(Dcutoffs)),
                        nrow = length(E1cutoffs), byrow = T)
rownames(sensit_matrix) <- E1cutoffs
colnames(sensit_matrix) <- Dcutoffs
rownames(precis_matrix) <- E1cutoffs
colnames(precis_matrix) <- Dcutoffs

for (e in E1cutoffs) {
  for (d in Dcutoffs) {
    
    nndgs <- nndgs %>% mutate(passc1 = D > d, passc2 = E1 > e)
    pass <- nndgs %>% filter(passc1 == TRUE, passc2 == TRUE) %>% select(bp, name)
    
    p <- lapply(unique(pass[,"name"]), function(name){
      signals <- sort(pass[pass$name == name, "bp"])
      truesig <- data.frame("start" = signals[1], "end" = signals[1], "name" = name)
      for(s in signals[-1]) {
        nl <- nrow(truesig)
        if(abs(truesig[nl, "end"] - s) < 25) {
          truesig[nl, "end"] = s
        } else {
          truesig <- rbind(truesig, list("start" = s, "end" = s, "name" = name))
        }
      }
      return(truesig)
    })
    
    p <- rbindlist(p)
    p <- p %>% mutate(start = start- 500, end = end - 500) %>%
      mutate(truepos = (start >= -150 & start <= 50)  | 
               (end >= -150 & end <= 50))
    truepositives <- p %>% filter(truepos == TRUE)
    falsepositives <- p %>% filter(truepos == FALSE)
    truepositives <- truepositives %>% group_by(name) %>% 
      summarise(cTSS = min(abs(start)), start = start[which.min(abs(start))])
    
    falsenegatives <- sum(!(promotores$name %in% truepositives$name))
    sensitivity = nrow(truepositives)/( nrow(truepositives) + falsenegatives)
    precision =  nrow(truepositives)/( nrow(truepositives) + nrow(falsepositives))
    precis_matrix[as.character(e), as.character(d)] = precision
    sensit_matrix[as.character(e), as.character(d)] = sensitivity
  }
}

sensit_melt <- melt(sensit_matrix, id.vars = c("E1", "D"),
                    measure.vars = "sensitivity")
names(sensit_melt) <- c("E1", "D", "sensitivity")
ggplot() + stat_contour(data = sensit_melt, 
                        aes(x = E1, y = D, z = sensitivity,
                        color = ..level..), 
                        breaks = seq(0, 1, 0.1)) +
  scale_color_continuous(name = "sensitivity") 


precis_melt <- melt(precis_matrix, id.vars = c("E1", "D"),
                    measure.vars = "precision")
names(precis_melt) <- c("E1", "D", "precision")
ggplot() + stat_contour(data = precis_melt, 
                        aes(x = E1, y = D, z = precision,
                            color = ..level..), 
                        breaks = seq(0, 1, 0.1)) +
          scale_color_continuous(name = "precision") 

