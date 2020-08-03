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


ident = ident_NKT
new_ident = as.vector(ident)
names(new_ident) = names(ident)
new_ident[new_ident %in% c("N1")] = "NKTp"
new_ident[new_ident %in% c("N2")] = "NKT1"
new_ident[new_ident %in% c("N3","N4","N5","N6")] = "NKT2"
new_ident[new_ident %in% c("N7")] = "NKT17"
new_ident = factor(new_ident)


plt = PlotDiversity(NKT_cdr3_pair_mat, new_ident)
print(plt)

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

clonotype_used = colnames(mat_cdr3_pair_NKT2)
mat_pseudotime_clonotype = matrix(data=0, ncol=nlevels(new_ident),
                                  nrow = length(clonotype_used))
rownames(mat_pseudotime_clonotype) = clonotype_used
colnames(mat_pseudotime_clonotype) = levels(new_ident)

for(i in 1:nlevels(new_ident)){
  cluster_cell_ind = names(new_ident)[new_ident == levels(new_ident)[i]]
  
  c_mat_cdr3_pair_NKT2 = mat_cdr3_pair_NKT2[rownames(mat_cdr3_pair_NKT2) %in% cluster_cell_ind,]
  c_mat_cdr3_pair_NKT2 = c_mat_cdr3_pair_NKT2[,clonotype_used]
  mat_pseudotime_clonotype[,levels(new_ident)[i]] = colSums(c_mat_cdr3_pair_NKT2)
}

cluster_shannon_indices = apply(mat_pseudotime_clonotype, 2,ComputeShannonIndex)

df_si = data.frame(cluster_shannon_indices)
df_si$cluster = rownames(df_si)
df_si = df_si[c("NKTp","NKT1","NKT2","NKT17"),]
df_si$cluster = factor(df_si$cluster, levels = c("NKTp","NKT1","NKT2","NKT17"))
plt = ggplot(df_si, aes(x=cluster,
                   y=cluster_shannon_indices,
                   group=1)) + geom_line() + geom_point(size=3)  +
  xlab("") +
  ylab("Normalized diversity Index") +
  ylim(low=0.4, high=1) +
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
  ) + ggsave(paste0("Figure7C_shannon_NKT.pdf"),width=7,height=7)

plt
