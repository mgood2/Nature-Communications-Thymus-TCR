# init
library(Stcr)
setwd("D:/elee4472/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.MAIT.RData") # seurat v2 file
load("../testdata/ident_MAIT.RData")
load("../testdata/ensemblGenes2018-10-15.RData")

ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_MAIT = ImportSeurat(sset,seurat.MAIT, version="v2")

MAIT_cdr3_mat = SetCDR3Mat(sset_MAIT, cell_ind = colnames(seurat.MAIT@scale.data),
                           chain1 = c("TRA"), chain2 = c("TRB"))

MAIT_cdr3_pair_mat = SetCDR3PairMat(sset_MAIT, MAIT_cdr3_mat, chain1 = c("TRA"), chain2 = c("TRB"))


plt = PlotFrequency(MAIT_cdr3_pair_mat, sset_MAIT@ident, canonical = FALSE,
                    chain1_v_gene = "TRAV1\\*",
                    chain1_j_gene = "TRAJ33", chain2_v_gene = "(TRBV13|TRBV19)" )
print(plt)

# mat_cdr3_pair_MAIT2 = MAIT_cdr3_pair_mat[[1]]
# df_clonotype = MAIT_cdr3_pair_mat[[2]]
# 
# S1_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M1")])
# S2_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M2","M3","M4")])
# MAIT1_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M5")])
# MAIT17_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M6","M7","M8")])
# 
# df_cdr3_pair_MAIT2 = data.frame(mat_cdr3_pair_MAIT2)
# df_cdr3_pair_MAIT2$celltype = ""
# df_cdr3_pair_MAIT2$celltype[rownames(df_cdr3_pair_MAIT2) %in% S1_ind] = "S1"
# df_cdr3_pair_MAIT2$celltype[rownames(df_cdr3_pair_MAIT2) %in% S2_ind] = "S2"
# df_cdr3_pair_MAIT2$celltype[rownames(df_cdr3_pair_MAIT2) %in% MAIT1_ind] = "MAIT1"
# df_cdr3_pair_MAIT2$celltype[rownames(df_cdr3_pair_MAIT2) %in% MAIT17_ind] = "MAIT17"
# 
# 
# for(i in 1:nrow(df_clonotype)){
#   df_clonotype$celltype[i] =  paste( unique(df_cdr3_pair_MAIT2[mat_cdr3_pair_MAIT2[,i] > 0,'celltype']), collapse = ",")
# }
# 
# 
# df_cdr3_pair_MAIT2$cluster = ident_MAIT[rownames(df_cdr3_pair_MAIT2)]
# for(i in 1:nrow(df_clonotype)){
#   df_clonotype$cluster[i] =  paste( unique(df_cdr3_pair_MAIT2[mat_cdr3_pair_MAIT2[,i] > 0,'cluster']), collapse = ",")
# }
# 
# df_clonotype2 = subset(df_clonotype, n > 0)
# df_clonotype2$canonical = "noncanonical"
# df_clonotype2[grepl("TRAV1\\*",df_clonotype2$chain1_v_gene) & grepl("TRAJ33",df_clonotype2$chain1_j_gene) &
#                 grepl("(TRBV13|TRBV19)",df_clonotype2$chain2_v_gene),]$canonical = "canonical"
# #write.csv(df_clonotype2, file="supple_clonotype_MAIT.csv")
# 
# canonicalList = list(AV=c("TRAV1*01"), AJ =c("TRAJ33*01"),BV = c("TRBV19*01", "TRBV13-1*02","TRBV13-2*01","TRBV13-3*01"))
# df_clonotype_labeled = df_clonotype2
# ident = ident_MAIT
# mat_cdr3_pair = mat_cdr3_pair_MAIT2
# 
# df_freq = ComputeFreq(mat_cdr3_pair,
#                       df_clonotype_labeled,
#                       ident,
#                       canonicalList)
# 
# df_freq[,c("cano_pct", "noncano_pct")]
# df_freq$cluster = rownames(df_freq)
# 
# 
# ggplot(df_freq, aes(x=cluster,y=noncano_pct*100)) +
#   geom_point(size=3) +
#   geom_line(group = 0) +
#   xlab("") +
#   ylab("% of non-canonical TCRs") +
#   ylim(0,50) # +
#   # ggsave("MAIT_noncanonical_frequency.pdf",width=7,height=7)
# 
# dev.off()


