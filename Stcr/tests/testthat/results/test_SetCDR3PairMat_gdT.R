# init
setwd("D:/Dropbox/R projects/Stcr/tests/testthat")
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


