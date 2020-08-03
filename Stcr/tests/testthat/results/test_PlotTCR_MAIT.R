# init
setwd("D:/elee4472/1.Projects/R projects/Stcr/tests/testthat")
library(Seurat)
library(Stcr)
library(ggplot2)
load("../testdata/seurat.MAIT.RData") # seurat v2 file
load("../testdata/ident_MAIT.RData")


ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_MAIT = ImportSeurat(sset,seurat.MAIT, version = "v2")

MAIT_cdr3_mat = SetCDR3Mat(sset_MAIT, cell_ind = colnames(seurat.MAIT@scale.data), chain1 = c("TRA"), chain2 = c("TRB"))

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
# df_clonotype2[grepl("TRAV1\\*",df_clonotype2$chain1_v_gene) & grepl("TRAJ33",df_clonotype2$chain1_j_gene) &
#                 grepl("(TRBV13|TRBV19)",df_clonotype2$chain2_v_gene),]$canonical = "canonical"

df_clonotype2[grepl("TRAV1\\*",df_clonotype2$chain1_v_gene) | grepl("TRAJ33",df_clonotype2$chain1_j_gene),]$canonical = "canonical"

write.csv(df_clonotype2, file="supple_clonotype_MAIT.csv")


clonotypes=as.vector(df_clonotype2$clonotype)[df_clonotype2$canonical == "canonical"]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_canonical.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[df_clonotype2$canonical == "noncanonical"]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_noncanonical.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRAV1\\*",df_clonotype2$chain1_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_TRAV1.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRAJ33",df_clonotype2$chain1_j_gene)]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_TRAJ33.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRBV19",df_clonotype2$chain2_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_TRBV19.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRBV13-1",df_clonotype2$chain2_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_TRBV13-1.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRBV13-2",df_clonotype2$chain2_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_TRBV13-2.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRBV13-3",df_clonotype2$chain2_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_TRBV13-3.pdf", width=7, height=7)



# major non canonical

clonotypes=as.vector(df_clonotype2$clonotype)[df_clonotype2$canonical == "noncanonical"]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

sort(table(subset(MAIT_cdr3_mat[[2]], cellName %in% clonotype_ind)$v_gene),decreasing = T)
sort(table(subset(MAIT_cdr3_mat[[2]], cellName %in% clonotype_ind)$d_gene),decreasing = T)
sort(table(subset(MAIT_cdr3_mat[[2]], cellName %in% clonotype_ind)$j_gene),decreasing = T)


clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRAV12D-2",df_clonotype2$chain1_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_TRAV12D-2.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRAV16",df_clonotype2$chain1_v_gene)]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_TRAV16.pdf", width=7, height=7)

clonotypes=as.vector(df_clonotype2$clonotype)[grepl("TRAJ17",df_clonotype2$chain1_j_gene)]
clonotype_ind = rownames(mat_cdr3_pair_MAIT2)[rowSums(mat_cdr3_pair_MAIT2[,clonotypes]) > 0]

plt = PlotTCR(sset_MAIT, cell_ind=clonotype_ind)
plt + ggsave("TCRplot_MAIT_TRAJ17.pdf", width=7, height=7)
