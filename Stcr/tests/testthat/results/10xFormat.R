setwd("D:/elee4472/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.MAIT.RData") # seurat v2 file
load("../testdata/ident_MAIT.RData") # seurat v2 file
load("../testdata/seurat.gdT.RData") # seurat v2 file
load("../testdata/ident_gdT.RData") # seurat v2 file


library(Seurat)
library(Matrix)


#MAIT
MAIT_F3_barcodes = paste0(substr(colnames(seurat.MAIT@data)[grepl("^F3",colnames(seurat.MAIT@data))],4,20),"-1")
MAIT_F5_barcodes = paste0(substr(colnames(seurat.MAIT@data)[grepl("^F5",colnames(seurat.MAIT@data))],4,20),"-1")

matrix_dir = "D:/elee4472/1.Projects/TCELL_eunmin/2. Single cell RNAseq/R_publication/data/F3/outs/raw_gene_bc_matrices/GRCm38.91/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
matF3 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(matF3) = barcode.names$V1
rownames(matF3) = feature.names$V1

MAIT_matF3 = matF3[rownames(seurat.MAIT@data),MAIT_F3_barcodes]
rm(matF3)
rm(feature.names)
rm(barcode.names)

matrix_dir = "D:/elee4472/1.Projects/TCELL_eunmin/2. Single cell RNAseq/R_publication/data/F3/outs/raw_gene_bc_matrices/GRCm38.91/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
matF5 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(matF5) = barcode.names$V1
rownames(matF5) = feature.names$V1

MAIT_matF5 = matF5[rownames(seurat.MAIT@data),MAIT_F5_barcodes]
rm(matF5)
rm(feature.names)
rm(barcode.names)

colnames(MAIT_matF3) = paste0("F3-",substr(colnames(MAIT_matF3),1,16))
colnames(MAIT_matF5) = paste0("F5-",substr(colnames(MAIT_matF5),1,16))

MAIT_mat= cbind(MAIT_matF3,MAIT_matF5)

writeMM(obj = MAIT_mat, file="MAIT_manual_filtered_outs/matrix.mtx")
write(x = rownames(MAIT_mat), file = "MAIT_manual_filtered_outs/genes.tsv")
write(x = colnames(MAIT_mat), file = "MAIT_manual_filtered_outs/barcodes.tsv")

projection = seurat.MAIT@dr$umap@cell.embeddings[colnames(MAIT_mat),]
colnames(projection) = c("UMAP-1", "UMAP-2")
write.csv(projection,row.names = TRUE, file = "MAIT_manual_filtered_outs/projection.csv", quote = FALSE)

clusters = as.data.frame(ident_MAIT[colnames(MAIT_mat)])
colnames(clusters) = c("Cluster")
write.csv(clusters,row.names = TRUE, file = "MAIT_manual_filtered_outs/clusters.csv", quote = FALSE)



#gdT
gdT_F3_barcodes = paste0(substr(colnames(seurat.gdT@data)[grepl("^F3",colnames(seurat.gdT@data))],4,20),"-1")
gdT_F5_barcodes = paste0(substr(colnames(seurat.gdT@data)[grepl("^F5",colnames(seurat.gdT@data))],4,20),"-1")


matrix_dir = "D:/elee4472/1.Projects/TCELL_eunmin/2. Single cell RNAseq/R_publication/data/F3/outs/raw_gene_bc_matrices/GRCm38.91/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
matF3 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(matF3) = barcode.names$V1
rownames(matF3) = feature.names$V1

gdT_matF3 = matF3[rownames(seurat.gdT@data),gdT_F3_barcodes]
rm(matF3)
rm(feature.names)
rm(barcode.names)

matrix_dir = "D:/elee4472/1.Projects/TCELL_eunmin/2. Single cell RNAseq/R_publication/data/F3/outs/raw_gene_bc_matrices/GRCm38.91/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
matF5 <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(matF5) = barcode.names$V1
rownames(matF5) = feature.names$V1

gdT_matF5 = matF5[rownames(seurat.gdT@data),gdT_F5_barcodes]
rm(matF5)
rm(feature.names)
rm(barcode.names)

colnames(gdT_matF3) = paste0("F3-",substr(colnames(gdT_matF3),1,16))
colnames(gdT_matF5) = paste0("F5-",substr(colnames(gdT_matF5),1,16))

gdT_mat= cbind(gdT_matF3,gdT_matF5)

writeMM(obj = gdT_mat, file="gdT_manual_filtered_outs/matrix.mtx")
write(x = rownames(gdT_mat), file = "gdT_manual_filtered_outs/genes.tsv")
write(x = colnames(gdT_mat), file = "gdT_manual_filtered_outs/barcodes.tsv")

projection = seurat.gdT@dr$umap@cell.embeddings[colnames(gdT_mat),]
colnames(projection) = c("UMAP-1", "UMAP-2")
write.csv(projection,row.names = TRUE, file = "gdT_manual_filtered_outs/projection.csv", quote = FALSE)

clusters = as.data.frame(ident_gdT[colnames(gdT_mat)])
colnames(clusters) = c("Cluster")
write.csv(clusters,row.names = TRUE, file = "gdT_manual_filtered_outs/clusters.csv", quote = FALSE)

