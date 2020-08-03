# init
setwd("D:/elee4472/1.Projects/R projects/Stcr/tests/testthat")
library(Stcr)
library(Seurat)
library(ggplot2)
load("../testdata/seurat.NKT.RData") # seurat v2 file
load("../testdata/ident_NKT.RData")


ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_NKT = ImportSeurat(sset,seurat.NKT)

NKT_cdr3_mat = SetCDR3Mat(sset_NKT, cell_ind = colnames(seurat.NKT@scale.data), chain1 = c("TRA"), chain2 = c("TRB"))

NKT_cdr3_pair_mat = SetCDR3PairMat(sset_NKT, NKT_cdr3_mat, chain1 = c("TRA"), chain2 = c("TRB"))

mat_cdr3_pair_NKT2 = NKT_cdr3_pair_mat[[1]]
df_clonotype = NKT_cdr3_pair_mat[[2]]

NKTp_ind = names(seurat.NKT@ident[seurat.NKT@ident %in% c("N1")])
NKT1_ind = names(seurat.NKT@ident[seurat.NKT@ident %in% c("N2")])
NKT2_ind = names(seurat.NKT@ident[seurat.NKT@ident %in% c("N3","N4","N5","N6")])
NKT17_ind = names(seurat.NKT@ident[seurat.NKT@ident %in% c("N7")])

NKT0_ind = c( "F3-TCGGGACCAAAGTGCG", "F3-TGTATTCTCTGGTTCC", "F5-CCTCAGTAGTGGGATC" ,"F5-GGAGCAAAGGTTCCTA",
              "F5-TTATGCTAGGCCCGTT","F5-ATCCGAAGTTCTGTTT", "F5-GGCTCGACAGGTCGTC", "F5-AGCTCTCAGAGGTTGC",
              "F3-CCTCTGAGTAGGCTGA", "F3-TCCACACGTTCGAATC" ,"F3-TCCACACAGCTGTCTA")

NKT0_cano = c("F3-TGTATTCTCTGGTTCC","F3-TCCACACGTTCGAATC","F3-TCCACACAGCTGTCTA")
NKT0_noncano = c("F5-CCTCAGTAGTGGGATC","F5-TTATGCTAGGCCCGTT")
NKT0_TCR_minus = c("F3-TCGGGACCAAAGTGCG","F3-CCTCTGAGTAGGCTGA","F5-AGCTCTCAGAGGTTGC","F5-GGCTCGACAGGTCGTC","F5-GGAGCAAAGGTTCCTA","F5-ATCCGAAGTTCTGTTT")


seurat.NKT@meta.data[NKT0_cano,"Yexprs"]
seurat.NKT@meta.data[NKT0_noncano,"Yexprs"]
seurat.NKT@meta.data[NKT0_TCR_minus,"Yexprs"]

df_cdr3_pair_NKT2 = data.frame(mat_cdr3_pair_NKT2)
df_cdr3_pair_NKT2$celltype = ""
df_cdr3_pair_NKT2$celltype[rownames(df_cdr3_pair_NKT2) %in% NKTp_ind] = "NKTp"
df_cdr3_pair_NKT2$celltype[rownames(df_cdr3_pair_NKT2) %in% NKT1_ind] = "NKT1"
df_cdr3_pair_NKT2$celltype[rownames(df_cdr3_pair_NKT2) %in% NKT2_ind] = "NKT2"
df_cdr3_pair_NKT2$celltype[rownames(df_cdr3_pair_NKT2) %in% NKT17_ind] = "NKT17"
df_cdr3_pair_NKT2$celltype[rownames(df_cdr3_pair_NKT2) %in% NKT0_ind] = "NKT0"



celltypes = df_cdr3_pair_NKT2$celltype
names(celltypes) = rownames(df_cdr3_pair_NKT2)
celltypes = factor(celltypes)

plt = PlotFrequency(NKT_cdr3_pair_mat, ident = celltypes, reverse = TRUE, selected_ind = NULL,
              chain1_v_gene = "TRAV11", chain1_d_gene = NULL, chain1_j_gene = "TRAJ18", 
              chain2_v_gene = "(TRBV13|TRBV29|TRBV1\\*)", chain2_d_gene = NULL, chain2_j_gene = NULL)

plt$data$cluster = factor(plt$data$cluster, levels=c("NKT0","NKTp","NKT1", "NKT2", "NKT17"))

plt +  theme_bw() + ylab("") +
  geom_point(size=4) +
  theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    plot.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), 
    axis.line = element_line(),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20),
    axis.text.x =element_text(colour=c("black","black", "lightgreen","darkred","blue"), angle = 90, vjust = 0.5),
    legend.position = "none") + ggsave("NKT_noncano_freq_added_NKT0.pdf",width=7,height=7)
dev.off()



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

canonicalList = list(AV=c("TRAV11*01", "TRAV11D*01"), AJ =c("TRAJ18*01"),BV = c("TRBV1*01","TRBV29*01","TRBV13-3*01","TRBV13-2*01","TRBV13-1*02"))
df_clonotype_labeled = df_clonotype2
ident = ident_NKT
mat_cdr3_pair = mat_cdr3_pair_NKT2



df_freq = ComputeFreq(mat_cdr3_pair,
                      df_clonotype_labeled,
                      ident,
                      canonicalList)
df_freq = df_freq[c("N1","N2","N3","N4","N5","N6","N7"),]
plt=PlotProportionsByGene(df_freq)
plt + ggsave("canonical_bar_NKT.pdf", width = 14, height= 7)
print(plt)
