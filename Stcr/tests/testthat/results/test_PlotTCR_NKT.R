# init
setwd("D:/elee4472/1.Projects/R projects/Stcr/tests/testthat")
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


clonotypes=as.vector(df_clonotype2$clonotype)[df_clonotype2$canonical == "canonical"]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_canonical.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[df_clonotype2$canonical == "noncanonical"]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_noncanonical.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRAV11",df_clonotype2$chain1_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_TRAV11_TRAV11D.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRAJ18",df_clonotype2$chain1_j_gene)]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_TRAJ18.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRAJ50",df_clonotype2$chain1_j_gene)]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_TRAJ50.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRAV13",df_clonotype2$chain1_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_TRAV13.pdf", width=7, height=7)

clonotypes=intersect(as.vector(df_clonotype2$clonotype)[grepl("TRAJ50",df_clonotype2$chain1_j_gene)],as.vector(df_clonotype2$clonotype)[grepl("TRAV13",df_clonotype2$chain1_v_gene)])
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_TRAV13_TRAJ50.pdf", width=7, height=7)



clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRBV1\\*",df_clonotype2$chain2_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_TRBV1.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRBV13-1",df_clonotype2$chain2_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_TRBV13-1.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRBV13-2",df_clonotype2$chain2_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_TRBV13-2.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRBV13-3",df_clonotype2$chain2_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_NKT2)[rowSums(mat_cdr3_pair_NKT2[,clonotypes]) > 0]

plt = PlotTCR(sset_NKT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_NKT_TRBV13-3.pdf", width=7, height=7)