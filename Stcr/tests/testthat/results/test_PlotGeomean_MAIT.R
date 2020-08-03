# init
setwd("D:/Dropbox/R projects/Stcr/tests/testthat")
load("../testdata/seurat.MAIT.RData") # seurat v2 file
load("../testdata/ident_MAIT.RData")


ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_MAIT = ImportSeurat(sset,seurat.MAIT)

MAIT_cdr3_mat = SetCDR3Mat(sset_MAIT, cell.ind = colnames(seurat.MAIT@scale.data), chain1 = c("TRA"), chain2 = c("TRB"))

MAIT_cdr3_pair_mat = SetCDR3PairMat(sset_MAIT, MAIT_cdr3_mat, chain1 = c("TRA"), chain2 = c("TRB"))

mat_cdr3_pair_MAIT2 = MAIT_cdr3_pair_mat[[1]]
df_clonotype = MAIT_cdr3_pair_mat[[2]]

S1_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M1")])
S2_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M2","M3","M4")])
MAIT1_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M5")])
MAIT17_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M6","M7","M8")])

df_cdr3_pair_MAIT2 = data.frame(mat_cdr3_pair_MAIT2)
df_cdr3_pair_MAIT2$celltype = ""
df_cdr3_pair_MAIT2$celltype[rownames(df_cdr3_pair_MAIT2) %in% S1_ind] = "S1"
df_cdr3_pair_MAIT2$celltype[rownames(df_cdr3_pair_MAIT2) %in% S2_ind] = "S2"
df_cdr3_pair_MAIT2$celltype[rownames(df_cdr3_pair_MAIT2) %in% MAIT1_ind] = "MAIT1"
df_cdr3_pair_MAIT2$celltype[rownames(df_cdr3_pair_MAIT2) %in% MAIT17_ind] = "MAIT17"


for(i in 1:nrow(df_clonotype)){
  df_clonotype$celltype[i] =  paste( unique(df_cdr3_pair_MAIT2[mat_cdr3_pair_MAIT2[,i] > 0,'celltype']), collapse = ",")
}


df_cdr3_pair_MAIT2$cluster = ident_MAIT[rownames(df_cdr3_pair_MAIT2)]
for(i in 1:nrow(df_clonotype)){
  df_clonotype$cluster[i] =  paste( unique(df_cdr3_pair_MAIT2[mat_cdr3_pair_MAIT2[,i] > 0,'cluster']), collapse = ",")
}

df_clonotype2 = subset(df_clonotype, n > 0)
df_clonotype2$canonical = "noncanonical"
df_clonotype2[grepl("TRAV1\\*",df_clonotype2$chain1_v_gene) & grepl("TRAJ33",df_clonotype2$chain1_j_gene) &
                grepl("(TRBV13|TRBV19)",df_clonotype2$chain2_v_gene),]$canonical = "canonical"

write.csv(df_clonotype2, file="supple_clonotype_MAIT.csv")




mat_cluster_clonotypes = matrix(data=0, ncol=nlevels(ident_MAIT),
                                nrow =ncol(mat_cdr3_pair_MAIT2))
colnames(mat_cluster_clonotypes) = levels(ident_MAIT)
rownames(mat_cluster_clonotypes) = colnames(mat_cdr3_pair_MAIT2)
for(i in 1:nlevels(ident_MAIT)){
  cluster_cell_ind = names(ident_MAIT)[ident_MAIT == levels(ident_MAIT)[i]]
  mat_cdr3_pair_MAIT2_cluster = mat_cdr3_pair_MAIT2[rownames(mat_cdr3_pair_MAIT2) %in% cluster_cell_ind,]
  mat_cluster_clonotypes[,levels(ident_MAIT)[i]] = colSums(mat_cdr3_pair_MAIT2_cluster)
}
rownames(mat_cluster_clonotypes) = df_clonotype$clonotype

colnames(mat_cluster_clonotypes) = levels(ident_MAIT)


geo_mean = matrix(data = 0, ncol = nlevels(ident_MAIT),
                  nrow = nlevels(ident_MAIT))
rownames(geo_mean) = levels(ident_MAIT)
colnames(geo_mean) = levels(ident_MAIT)


for(i in 1:nlevels(ident_MAIT)){
  for(j in 1:nlevels(ident_MAIT)){
    if (i != j){
      sumA = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),i]/sum(mat_cluster_clonotypes[,i]))
      sumB = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),j]/sum(mat_cluster_clonotypes[,j]))
      geo_mean[i,j] = sqrt(sumA*sumB)
    }
  }
}

geo_mean = geo_mean[c("M1","M2","M3","M4","M5","M6","M7","M8"),c("M1","M2","M3","M4","M5","M6","M7","M8")]

write.csv(geo_mean, file="MAIT_Geomean_value.csv")
mat_breaks <- seq(0, 0.38, length.out = 50)

pheatmap::pheatmap(geo_mean, color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name =
                                                             "RdBu")))(100)[51:100], breaks = mat_breaks,
         cluster_rows = F, cluster_cols = F, angle_col = c("0"), filename = "Figure7A_clonotypic_overlap_MAIT.pdf" )

