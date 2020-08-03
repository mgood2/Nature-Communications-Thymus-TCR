# init
setwd("D:/elee4472/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.gdT.RData") # seurat v2 file
load("../testdata/ident_gdT.RData")
ssetF3 = Load10xVDJ('../testdata/mTCR-F12-gdT/out_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-G1-gdT/out_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5)

sset_gdT = ImportSeurat(sset,seurat.gdT)

gdT_cdr3_mat = SetCDR3Mat(sset_gdT, cell.ind = colnames(seurat.gdT@scale.data), chain1 = c("TRG"), chain2 = c("TRD","Multi"), v_gene2 = "DV",
                          j_gene2="DJ")

gdT_cdr3_pair_mat = SetCDR3PairMat(sset_gdT, gdT_cdr3_mat, chain1 = c("TRG"), chain2 = c("TRD","Multi"))

plt = PlotClusterHeatmap(gdT_cdr3_pair_mat, ident_gdT, compute = "geomean")
print(plt)
plt2 = PlotClusterHeatmap(gdT_cdr3_pair_mat, ident_gdT, compute = "logFC")
print(plt2)

mat_cdr3_pair_gdT2 = gdT_cdr3_pair_mat[[1]]
df_clonotype = gdT_cdr3_pair_mat[[2]]

gdTp_ind = names(ident_gdT[ident_gdT %in% c("G1","G2","G3")])
gdT17i_ind = names(ident_gdT[ident_gdT %in% c("G4","G5")])
gdT17_ind = names(ident_gdT[ident_gdT %in% c("G6-1","G6-2")])
gdT1i_ind = names(ident_gdT[ident_gdT %in% c("G7-1")])
gdT1_ind = names(ident_gdT[ident_gdT %in% c("G7-2")])

df_cdr3_pair_gdT2 = data.frame(mat_cdr3_pair_gdT2)
df_cdr3_pair_gdT2$celltype = ""
df_cdr3_pair_gdT2[rownames(df_cdr3_pair_gdT2) %in% gdTp_ind,]$celltype = "gdTp"
df_cdr3_pair_gdT2[rownames(df_cdr3_pair_gdT2) %in% gdT1i_ind,]$celltype = "gdT1i"
df_cdr3_pair_gdT2[rownames(df_cdr3_pair_gdT2) %in% gdT1_ind,]$celltype = "gdT1"
df_cdr3_pair_gdT2[rownames(df_cdr3_pair_gdT2) %in% gdT17i_ind,]$celltype = "gdT17i"
df_cdr3_pair_gdT2[rownames(df_cdr3_pair_gdT2) %in% gdT17_ind,]$celltype = "gdT17"

for(i in 1:nrow(df_clonotype)){
  df_clonotype$celltype[i] =  paste( unique(df_cdr3_pair_gdT2[mat_cdr3_pair_gdT2[,i] > 0,'celltype']), collapse = ",")
}


df_cdr3_pair_gdT2$cluster = ident_gdT[rownames(df_cdr3_pair_gdT2)]
for(i in 1:nrow(df_clonotype)){
  df_clonotype$cluster[i] =  paste( unique(df_cdr3_pair_gdT2[mat_cdr3_pair_gdT2[,i] > 0,'cluster']), collapse = ",")
}

df_clonotype2 = subset(df_clonotype, n > 0)

write.csv(df_clonotype2, file="supple_clonotype_gdT.csv")




mat_cluster_clonotypes = matrix(data=0, ncol=nlevels(ident_gdT),
                                nrow =ncol(mat_cdr3_pair_gdT2))
colnames(mat_cluster_clonotypes) = levels(ident_gdT)
rownames(mat_cluster_clonotypes) = colnames(mat_cdr3_pair_gdT2)
for(i in 1:nlevels(ident_gdT)){
  cluster_cell_ind = names(ident_gdT)[ident_gdT == levels(ident_gdT)[i]]
  mat_cdr3_pair_gdT2_cluster = mat_cdr3_pair_gdT2[rownames(mat_cdr3_pair_gdT2) %in% cluster_cell_ind,]
  mat_cluster_clonotypes[,levels(ident_gdT)[i]] = colSums(mat_cdr3_pair_gdT2_cluster)
}
rownames(mat_cluster_clonotypes) = df_clonotype$clonotype

colnames(mat_cluster_clonotypes) = levels(ident_gdT)


geo_mean = matrix(data = 0, ncol = nlevels(ident_gdT),
                  nrow = nlevels(ident_gdT))
rownames(geo_mean) = levels(ident_gdT)
colnames(geo_mean) = levels(ident_gdT)


for(i in 1:nlevels(ident_gdT)){
  for(j in 1:nlevels(ident_gdT)){
    if (i != j){
      sumA = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),i]/sum(mat_cluster_clonotypes[,i]))
      sumB = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),j]/sum(mat_cluster_clonotypes[,j]))
      geo_mean[i,j] = sqrt(sumA*sumB)
    }
  }
}

geo_mean = geo_mean[c("G1","G2","G3","G4","G5","G6-1","G6-2","G7-1","G7-2"),c("G1","G2","G3","G4","G5","G6-1","G6-2","G7-1","G7-2")]

write.csv(geo_mean, file="gdT_Geomean_value.csv")
mat_breaks <- seq(0, 0.38, length.out = 50)

pheatmap::pheatmap(geo_mean, color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name =
                                                             "RdBu")))(100)[51:100], breaks = mat_breaks,
         cluster_rows = F, cluster_cols = F, angle_col = c("0"), filename = "Figure7A_clonotypic_overlap_gdT.pdf" )

