library(ComplexHeatmap)
library(circlize)

setwd("/home/diana/Workspace/PhD/semestre2/bioinfo_estructural/Tarea_2")

pssm <- read.table("P1_adap.pssm", nrows = 86, header = T)
m1 <- pssm[, 2:22]
rownames(m1) <- paste(rownames(m1), m1[,1], sep = "-")
m1 <- data.matrix(m1[, -1])
min(m1)

positivos <- which(m1 > 0, arr.ind = T)
positivos[ order(as.numeric(row.names(positivos))),] 

max(m)

Heatmap(m1, col = colorRamp2(c(min(m1), 10), c("tan1", "blue")), column_title = "Matriz 1",   
        cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = T, show_column_names = T, show_row_names = T)

m2 <-  pssm[, 23:42]
rownames(m2) <- rownames(m1)
colnames(m2) <- colnames(m1)
min(m2)
max(m2)
Heatmap(m2, col = colorRamp2(c(min(m2), 60), c("tan1", "blue")), column_title = "Matriz 2",   
        cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = T, show_column_names = T, show_row_names = T)
