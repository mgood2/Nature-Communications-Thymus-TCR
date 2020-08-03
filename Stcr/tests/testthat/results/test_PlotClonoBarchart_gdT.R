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


cellcount_clonotypes = colSums(mat_cdr3_pair_gdT2)

df = data.frame(cellcount = cellcount_clonotypes, clonotypes = names(cellcount_clonotypes), stringsAsFactors = F)
sub_df = subset(df, cellcount > 1)
sub_df = sub_df[rev(order(sub_df$cellcount)),]

long_cellcount_clonotypes =  reshape2::melt(mat_cdr3_pair_gdT2[,rownames(sub_df)])
colnames(long_cellcount_clonotypes) = c("cell","clonotype","cellcounts")
long_cellcount_clonotypes = long_cellcount_clonotypes[long_cellcount_clonotypes$cellcounts > 0,]

clust = as.vector(ident_gdT)
names(clust) = names(ident_gdT)

long_cellcount_clonotypes$cluster = ""
for(i in 1:length(long_cellcount_clonotypes$cell)){
  long_cellcount_clonotypes$cluster[i] = clust[as.character(long_cellcount_clonotypes$cell[i])]
}
library(dplyr)
summary_df = long_cellcount_clonotypes %>% group_by(clonotype,cluster) %>% count(cellcounts) %>% as.data.frame()
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ggcolors = gg_color_hue(7)


ggcolors=c(ggcolors[1],"#C04A00","#E88A00",ggcolors[3:6],"#f59bff","#B748C4")
ggcolors=ggcolors[c(5,4,6,1,7,8,9,3,2)]
ggplot(summary_df, aes(x = factor(clonotype, levels=unique(summary_df$clonotype)
), y =n, fill = factor(cluster))) +
  scale_fill_manual(values=ggcolors) +
  geom_bar(stat="identity") + 
  ylab("cell counts") +
  xlab("clonotype") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)) +
  ggsave("Figure7B_bargraph_clonotype_gdT.pdf", width = 20)
dev.off()


