setwd("/home/diana/Workspace/PhD/semestre2/bioinfo_estructural/")

code <- list("Ala"="A", "Cys"="C","Asp"="D", "Glu"="E", "Phe"="F", 
      "Gly"="G","His"="H", "Ile"="I","Lys"="K", "Leu"="L","Met"="M",
      "Asn"="N","Pro"="P", "Gln"="Q","Arg"="R", "Ser"="S","Thr"="T", 
      "Val"="V", "Trp"="W", "Tyr"="Y")
names(code) <- toupper(names(code))

lines <- readLines("data/1ayi.pdb")
ids <- grep("^ATOM*", lines)
atom.lines <- lines[ids]
atom.lines <- lapply(atom.lines, function(line){
      spl <- strsplit(line, split= " ")[[1]]
      spl <- spl[which(spl != "")]
      return(c(spl[4], spl[6]))
  })
atom.lines <- as.data.frame(do.call(rbind, atom.lines))
colnames(atom.lines) <- c("amino", "pos")
atom.lines <- atom.lines[!duplicated(atom.lines$pos), ]
fasta <- apply(atom.lines, 1, function(al){
  code[al["amino"]][[1]]
})
fasta <- paste(fasta, collapse =  "")
fasta <- paste(">", fasta, sep = "\n")
write(fasta, file = "P1.faa")
