# init
setwd("D:/Dropbox/R projects/Stcr/tests/testthat")
load("../testdata/seurat.NKT.RData") # seurat v2 file
load("../testdata/ident_NKT.RData")


ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")

sset = MergeSamples(ssetF3,ssetF5) 

sset_NKT = ImportSeurat(sset,seurat.NKT)

NKT_cdr3_mat = SetCDR3Mat(sset_NKT, cell.ind = colnames(seurat.NKT@scale.data), chain1 = c("TRA"), chain2 = c("TRB"))

NKT_cdr3_pair_mat = SetCDR3PairMat(sset_NKT, NKT_cdr3_mat, chain1 = c("TRA"), chain2 = c("TRB"))

mat_cdr3_pair_NKT2 = NKT_cdr3_pair_mat[[1]]
df_clonotype = NKT_cdr3_pair_mat[[2]]

NKTp_ind = names(seurat.NKT@ident[seurat.NKT@ident %in% c("N1")])
NKT1_ind = names(seurat.NKT@ident[seurat.NKT@ident %in% c("N2")])
NKT2_ind = names(seurat.NKT@ident[seurat.NKT@ident %in% c("N3","N4","N5","N6")])
NKT17_ind = names(seurat.NKT@ident[seurat.NKT@ident %in% c("N7")])

df_cdr3_pair_NKT2 = data.frame(mat_cdr3_pair_NKT2)
df_cdr3_pair_NKT2$celltype = ""
df_cdr3_pair_NKT2$celltype[rownames(df_cdr3_pair_NKT2) %in% NKTp_ind] = "NKTp"
df_cdr3_pair_NKT2$celltype[rownames(df_cdr3_pair_NKT2) %in% NKT1_ind] = "NKT1"
df_cdr3_pair_NKT2$celltype[rownames(df_cdr3_pair_NKT2) %in% NKT2_ind] = "NKT2"
df_cdr3_pair_NKT2$celltype[rownames(df_cdr3_pair_NKT2) %in% NKT17_ind] = "NKT17"

for(i in 1:nrow(df_clonotype)){
  df_clonotype$celltype[i] =  paste( unique(df_cdr3_pair_NKT2[mat_cdr3_pair_NKT2[,i] > 0,'celltype']), collapse = ",")
}

df_cdr3_pair_NKT2$cluster = ident_NKT[rownames(df_cdr3_pair_NKT2)]
for(i in 1:nrow(df_clonotype)){
  df_clonotype$cluster[i] =  paste( unique(df_cdr3_pair_NKT2[mat_cdr3_pair_NKT2[,i] > 0,'cluster']), collapse = ",")
}

df_clonotype2 = subset(df_clonotype, n > 0)
df_clonotype2$canonical = "noncanonical"
df_clonotype2[grepl("TRAV11",df_clonotype2$chain1_v_gene) & grepl("TRAJ18",df_clonotype2$chain1_j_gene) &
                grepl("(TRBV13|TRBV29|TRBV1\\*)",df_clonotype2$chain2_v_gene),]$canonical = "canonical"

write.csv(df_clonotype2, file="supple_clonotype_NKT.csv")


ident = ident_NKT
new_ident = as.vector(ident)
names(new_ident) = names(ident)
new_ident[new_ident %in% c("N1")] = "NKTp"
new_ident[new_ident %in% c("N2")] = "NKT1"
new_ident[new_ident %in% c("N3","N4","N5","N6")] = "NKT2"
new_ident[new_ident %in% c("N7")] = "NKT17"
new_ident = factor(new_ident)


mat_cluster_clonotypes = matrix(data=0, ncol=nlevels(new_ident),
                                nrow =ncol(mat_cdr3_pair_NKT2))
colnames(mat_cluster_clonotypes) = levels(new_ident)
rownames(mat_cluster_clonotypes) = colnames(mat_cdr3_pair_NKT2)
for(i in 1:nlevels(new_ident)){
  cluster_cell_ind = names(new_ident)[new_ident == levels(new_ident)[i]]
  mat_cdr3_pair_NKT2_cluster = mat_cdr3_pair_NKT2[rownames(mat_cdr3_pair_NKT2) %in% cluster_cell_ind,]
  mat_cluster_clonotypes[,levels(new_ident)[i]] = colSums(mat_cdr3_pair_NKT2_cluster)
}
rownames(mat_cluster_clonotypes) = df_clonotype$clonotype

colnames(mat_cluster_clonotypes) = levels(new_ident)


geo_mean = matrix(data = 0, ncol = nlevels(new_ident),
                  nrow = nlevels(new_ident))
rownames(geo_mean) = levels(new_ident)
colnames(geo_mean) = levels(new_ident)


for(i in 1:nlevels(new_ident)){
  for(j in 1:nlevels(new_ident)){
    if (i != j){
      sumA = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),i]/sum(mat_cluster_clonotypes[,i]))
      sumB = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),j]/sum(mat_cluster_clonotypes[,j]))
      geo_mean[i,j] = sqrt(sumA*sumB)
    }
  }
}

geo_mean = geo_mean[c("NKTp","NKT1","NKT2","NKT17"),c("NKTp","NKT1","NKT2","NKT17")]


mat_breaks <- seq(0, 0.38, length.out = 50)

write.csv(geo_mean, file="NKT_Geomean_value.csv")

pheatmap::pheatmap(geo_mean, color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name =
                                                             "RdBu")))(100)[51:100], breaks = mat_breaks,
         cluster_rows = F, cluster_cols = F, angle_col = c("0"), filename = "Figure7A_clonotypic_overlap_NKT.pdf" )

