setwd("/home/diana/Workspace/PhD/semestre2/bioinfo_estructural/Tarea_2/")

lines <- readLines("P1_multifasta.faa")
rama <- lines[4]
psiblast <- lines[6]

rama <- paste("*", rama, "*", sep = "")
nchar(rama) == nchar(psiblast)
len <- nchar(rama)
rama <- strsplit(rama, "")[[1]]
psiblast <- strsplit(psiblast, "")[[1]]

porc.ident <- (sum(rama == psiblast)/len) * 100
Hs <- (sum(psiblast[rama == psiblast] == "H")/len)* 100
Cs <- (sum(psiblast[rama == psiblast] == "C")/len)* 100
Es <- (sum(psiblast[rama == psiblast] == "E")/len)* 100

library(ComplexHeatmap)
library(circlize)


ba = HeatmapAnnotation(
  text = anno_text(fpsiblast, rot = 90, offset = unit(1, "npc"), just = "right"),
  annotation_height = max_text_width(fpsiblast)
)

ta <- HeatmapAnnotation(
  text = anno_text(rama, rot = 90, offset = unit(1, "npc"), just = "right"),
  annotation_height = max_text_width(rama)
)

Heatmap(cbind(rama, psiblast), name = "seq_comparison", col = c("*" = "tan1", "C" = "skyblue3", "E" = "violetred1", "H" = "springgreen3"), 
        cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = T, column_title = "Secuencias",  
        show_column_names = F, show_row_names = F,rect_gp = gpar(col= "white"), heatmap_legend_param = list(title = "char"))

