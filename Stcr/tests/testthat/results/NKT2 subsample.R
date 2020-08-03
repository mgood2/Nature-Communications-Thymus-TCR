# init
library(Seurat)
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.NKT.RData") # seurat v2 file
load("../testdata/ident_NKT.RData")
load("../testdata/ensemblGenes2018-10-15.RData")





DimPlot(seurat.NKT, reduction.use="umap") + theme(axis.ticks = element_blank(),
                                                  axis.title = element_blank(),
                                                  axis.text = element_blank()) + ggsave("subsampled_NKT2_only_100.pdf", height=7, width= 7)


set.seed(10)
sampled_NKT2_cells = sample(NKT2_cells, as.integer(length(NKT2_cells)*0.75))

# NKT_cells = c(sampled_NKT2_cells)
NKT_cells = c(colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N1", "N2","N7")], sampled_NKT2_cells)
sub_seurat.NKT75 = SubsetData(seurat.NKT, cells.use = NKT_cells )

PCAuse = 10
set.seed(42)
sub_seurat.NKT75 <- FindClusters(sub_seurat.NKT75,
                                 resolution = 0.4,
                                 reduction.type="pca", dims.use = 1:PCAuse, save.SNN = TRUE, force.recalc = TRUE)

DimPlot(sub_seurat.NKT75, reduction.use = "umap")+ theme(axis.ticks = element_blank(),
                                                             axis.title = element_blank(),
                                                             axis.text = element_blank())  + ggsave("subsampled_NKT2_only_75.pdf", height=7, width= 7)




# NKT75cells = c(colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N1", "N2","N7")], NKT_cells)



# newcluster = as.vector(seurat.NKT@ident)
# names(newcluster) = names(seurat.NKT@ident)
# newcluster = newcluster[NKT75cells]
# newcluster[NKT_cells] = as.vector(sub_seurat.NKT75@ident)
# 
# sub_seurat.NKT75_all = SubsetData(seurat.NKT, cells.use = NKT75cells )
# sub_seurat.NKT75_all@ident = factor(newcluster)
# 
# DimPlot(sub_seurat.NKT75_all, reduction.use = "umap")+ theme(axis.ticks = element_blank(),
#                                                              axis.title = element_blank(),
#                                                              axis.text = element_blank())  + ggsave("subsampled_NKT2_only_75.pdf", height=7, width= 7)


# NKT2 cell only cluster

set.seed(10)
sampled_NKT2_cells = sample(NKT2_cells, as.integer(length(NKT2_cells)*0.50))

# NKT_cells = c(sampled_NKT2_cells)
NKT_cells = c(colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N1", "N2","N7")], sampled_NKT2_cells)

sub_seurat.NKT50 = SubsetData(seurat.NKT, cells.use = NKT_cells )

PCAuse = 10
set.seed(42)
sub_seurat.NKT50 <- FindClusters(sub_seurat.NKT50,
                                 resolution = 0.4,
                               reduction.type="pca", dims.use = 1:PCAuse, save.SNN = TRUE, force.recalc = TRUE)

DimPlot(sub_seurat.NKT50, reduction.use = "umap")+ theme(axis.ticks = element_blank(),
                                                             axis.title = element_blank(),
                                                             axis.text = element_blank())  + ggsave("subsampled_NKT2_only_50.pdf", height=7, width= 7)

# 
# 
# NKT50cells = c(colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N1", "N2","N7")], NKT_cells)
# 
# 
# newcluster = as.vector(seurat.NKT@ident)
# names(newcluster) = names(seurat.NKT@ident)
# newcluster = newcluster[NKT50cells]
# newcluster[NKT_cells] = as.vector(sub_seurat.NKT50@ident)
# 
# sub_seurat.NKT50_all = SubsetData(seurat.NKT, cells.use = NKT50cells )
# sub_seurat.NKT50_all@ident = factor(newcluster)
# 
# DimPlot(sub_seurat.NKT50_all, reduction.use = "umap")+ theme(axis.ticks = element_blank(),
#                                                              axis.title = element_blank(),
#                                                              axis.text = element_blank())  + ggsave("subsampled_NKT2_only_50.pdf", height=7, width= 7)
# 


set.seed(10)
sampled_NKT2_cells = sample(NKT2_cells, as.integer(length(NKT2_cells)*0.25))

# NKT_cells = c(sampled_NKT2_cells)
NKT_cells = c(colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N1", "N2","N7")], sampled_NKT2_cells)
sub_seurat.NKT25 = SubsetData(seurat.NKT, cells.use = NKT_cells )

PCAuse = 10
set.seed(42)
sub_seurat.NKT25 <- FindClusters(sub_seurat.NKT25,
                                 resolution = 0.4,
                                 reduction.type="pca", dims.use = 1:PCAuse, save.SNN = TRUE, force.recalc = TRUE)
DimPlot(sub_seurat.NKT25, reduction.use = "umap") + theme(axis.ticks = element_blank(),
                                                          axis.title = element_blank(),
                                                          axis.text = element_blank()) + ggsave("subsampled_NKT2_only_25.pdf", height=7, width= 7)
# 
# 
# NKT25cells = c(colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N1", "N2","N7")], NKT_cells)
# 
# 
# newcluster = as.vector(seurat.NKT@ident)
# names(newcluster) = names(seurat.NKT@ident)
# newcluster = newcluster[NKT25cells]
# newcluster[NKT_cells] = as.vector(sub_seurat.NKT25@ident)
# 
# sub_seurat.NKT25_all = SubsetData(seurat.NKT, cells.use = NKT25cells )
# sub_seurat.NKT25_all@ident = factor(newcluster)
# 
# DimPlot(sub_seurat.NKT25_all, reduction.use = "umap") + theme(axis.ticks = element_blank(),
#                                                               axis.title = element_blank(),
#                                                               axis.text = element_blank()) + ggsave("subsampled_NKT2_only_25.pdf", height=7, width= 7)



set.seed(10)
sampled_NKT2_cells = sample(NKT2_cells, as.integer(length(NKT2_cells)*0.10))

# NKT_cells = c(sampled_NKT2_cells)
NKT_cells = c(colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N1", "N2","N7")], sampled_NKT2_cells)
sub_seurat.NKT10 = SubsetData(seurat.NKT, cells.use = NKT_cells )

PCAuse = 10
set.seed(42)
sub_seurat.NKT10 <- FindClusters(sub_seurat.NKT10,
                                 resolution = 0.4,
                                 reduction.type="pca", dims.use = 1:PCAuse, save.SNN = TRUE, force.recalc = TRUE)

NKT10cells = c(colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N1", "N2","N7")], NKT_cells)


DimPlot(sub_seurat.NKT10, reduction.use = "umap") + theme(axis.ticks = element_blank(),
                                                              axis.title = element_blank(),
                                                              axis.text = element_blank()) + ggsave("subsampled_NKT2_only_10.pdf", height=7, width= 7)

# newcluster = as.vector(seurat.NKT@ident)
# names(newcluster) = names(seurat.NKT@ident)
# newcluster = newcluster[NKT25cells]
# newcluster[NKT_cells] = as.vector(sub_seurat.NKT10@ident)
# 
# sub_seurat.NKT10_all = SubsetData(seurat.NKT, cells.use = NKT10cells )
# sub_seurat.NKT10_all@ident = factor(newcluster)
# 
# DimPlot(sub_seurat.NKT10_all, reduction.use = "umap") + theme(axis.ticks = element_blank(),
#                                                               axis.title = element_blank(),
#                                                               axis.text = element_blank()) + ggsave("subsampled_NKT2_only_10.pdf", height=7, width= 7)

# 
# 
# NKT2_cells = colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N3", "N4","N5", "N6")]
# table(seurat.NKT@ident)
# set.seed(10)
# sampled_NKT2_cells = sample(NKT2_cells, 240)
# 
# NKT_cells = c(colnames(seurat.NKT@scale.data)[seurat.NKT@ident %in% c("N1", "N2","N7")], sampled_NKT2_cells)
# 
# sub_seurat.NKT = SubsetData(seurat.NKT, cells.use = NKT_cells )
# 
# PCAuse = 10
# set.seed(42)
# sub_seurat.NKT <- FindClusters(sub_seurat.NKT,
#                                resolution = 0.6,
#                                reduction.type="pca", dims.use = 1:PCAuse, save.SNN = TRUE, force.recalc = TRUE)
# DimPlot(sub_seurat.NKT, reduction.use = "umap") + ggsave("subsampled_NKT_240cells.pdf", height=7, width= 7)

