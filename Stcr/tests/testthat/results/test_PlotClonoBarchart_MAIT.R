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


cellcount_clonotypes = colSums(mat_cdr3_pair_MAIT2)

df = data.frame(cellcount = cellcount_clonotypes, clonotypes = names(cellcount_clonotypes), stringsAsFactors = F)
sub_df = subset(df, cellcount > 1)
sub_df = sub_df[rev(order(sub_df$cellcount)),]

long_cellcount_clonotypes =  reshape2::melt(mat_cdr3_pair_MAIT2[,rownames(sub_df)])
colnames(long_cellcount_clonotypes) = c("cell","clonotype","cellcounts")
long_cellcount_clonotypes = long_cellcount_clonotypes[long_cellcount_clonotypes$cellcounts > 0,]

long_cellcount_clonotypes$cluster = ""
clust = as.vector(ident_MAIT)
names(clust) = names(ident_MAIT)


for(i in 1:length(long_cellcount_clonotypes$cell)){
  long_cellcount_clonotypes$cluster[i] = clust[as.character(long_cellcount_clonotypes$cell[i])]
}
library(dplyr)
summary_df = long_cellcount_clonotypes %>% group_by(clonotype,cluster) %>% count(cellcounts) %>% as.data.frame()

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ggcolors = gg_color_hue(nlevels(ident_MAIT))
ggcolors = ggcolors[c(2,8,1,5,6,7,4,3)]
ggplot(summary_df, aes(x = factor(clonotype, levels=unique(summary_df$clonotype)
), y =n, fill = factor(cluster))) +
  scale_fill_manual(values=ggcolors) +
  geom_bar(stat="identity") + 
  ylab("cell counts") +
  xlab("clonotype") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # legend.position = "none"
        ) +
  ggsave("Figure7B_bargraph_clonotype_MAIT.pdf", width = 20)
dev.off()
