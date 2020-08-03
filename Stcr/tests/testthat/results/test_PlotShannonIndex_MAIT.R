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


clonotype_used = colnames(mat_cdr3_pair_MAIT2)
mat_pseudotime_clonotype = matrix(data=0, ncol=nlevels(ident_MAIT),
                                  nrow = length(clonotype_used))
rownames(mat_pseudotime_clonotype) = clonotype_used
colnames(mat_pseudotime_clonotype) = levels(ident_MAIT)

for(i in 1:nlevels(ident_MAIT)){
  cluster_cell_ind = names(ident_MAIT)[ident_MAIT == levels(ident_MAIT)[i]]
  
  c_mat_cdr3_pair_MAIT2 = mat_cdr3_pair_MAIT2[rownames(mat_cdr3_pair_MAIT2) %in% cluster_cell_ind,]
  c_mat_cdr3_pair_MAIT2 = c_mat_cdr3_pair_MAIT2[,clonotype_used]
  mat_pseudotime_clonotype[,levels(ident_MAIT)[i]] = colSums(c_mat_cdr3_pair_MAIT2)
}

cluster_shannon_indices = apply(mat_pseudotime_clonotype, 2,ComputeShannonIndex)

df_si = data.frame(cluster_shannon_indices)
df_si$cluster = rownames(df_si)
df_si = df_si[c("M1","M2","M3","M4","M5","M6","M7","M8"),]
df_si$cluster = factor(df_si$cluster, levels = c("M1","M2","M3","M4","M5","M6","M7","M8"))
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
  ) + ggsave(paste0("Figure7C_shannon_MAIT.pdf"),width=7,height=7)

plt
