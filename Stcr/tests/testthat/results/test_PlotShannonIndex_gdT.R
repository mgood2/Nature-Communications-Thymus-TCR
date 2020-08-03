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


clonotype_used = colnames(mat_cdr3_pair_gdT2)
mat_pseudotime_clonotype = matrix(data=0, ncol=nlevels(ident_gdT),
                                  nrow = length(clonotype_used))
rownames(mat_pseudotime_clonotype) = clonotype_used
colnames(mat_pseudotime_clonotype) = levels(ident_gdT)

for(i in 1:nlevels(ident_gdT)){
  cluster_cell_ind = names(ident_gdT)[ident_gdT == levels(ident_gdT)[i]]
  
  c_mat_cdr3_pair_gdT2 = mat_cdr3_pair_gdT2[rownames(mat_cdr3_pair_gdT2) %in% cluster_cell_ind,]
  c_mat_cdr3_pair_gdT2 = c_mat_cdr3_pair_gdT2[,clonotype_used]
  mat_pseudotime_clonotype[,levels(ident_gdT)[i]] = colSums(c_mat_cdr3_pair_gdT2)
}

cluster_shannon_indices = apply(mat_pseudotime_clonotype, 2,ComputeShannonIndex)

df_si = data.frame(cluster_shannon_indices)
df_si$cluster = rownames(df_si)
df_si = df_si[c("G1","G2","G3","G4","G5","G6-1","G6-2","G7-1","G7-2"),]
df_si$cluster = factor(df_si$cluster, levels = c("G1","G2","G3","G4","G5","G6-1","G6-2","G7-1","G7-2"))
plt = ggplot(df_si, aes(x=cluster,
                   y=cluster_shannon_indices,
                   group=1)) + geom_line() + geom_point(size=3)  +
  xlab("") +
  ylab("Normalized diversity Index") +
  ylim(low=0.3, high=1) +
  # scale_color_manual(values=c("#F8766D","#00BE67","#00A9FF")) +
  theme_bw() +
  theme(
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=20),
    plot.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x =element_blank()
  ) + ggsave(paste0("Figure7C_shannon_gdT.pdf"),width=7,height=7)

plt
