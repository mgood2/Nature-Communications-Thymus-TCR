# init
library(ggplot2)
setwd("D:/elee4472/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.MAIT.RData") # seurat v2 file
load("../testdata/ident_MAIT.RData")


ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_MAIT = ImportSeurat(sset,seurat.MAIT, version="v2")

MAIT_cdr3_mat = SetCDR3Mat(sset_MAIT, cell_ind = colnames(seurat.MAIT@scale.data),
                           chain1 = c("TRA"), chain2 = c("TRB"))

MAIT_cdr3_pair_mat = SetCDR3PairMat(sset_MAIT, MAIT_cdr3_mat, chain1 = c("TRA"), chain2 = c("TRB"))

plt = PlotFrequency(MAIT_cdr3_pair_mat, ident = seurat.MAIT@ident,
                    reverse = FALSE,
                    chain1_v_gene = "TRAV1\\*",
                    chain1_j_gene = "TRAJ33")
plt + ggplot2::ylab("% of canonical TCRs") +
  theme_bw() +
  theme(line = element_blank()) +
  ggplot2::ggsave("canonical_v_j_alpha_only_MAIT.pdf")

plt = PlotFrequencyByGene(MAIT_cdr3_pair_mat, ident = sset_MAIT@ident, canonical = TRUE,
                    chain1_v_gene = "TRAV1\\*", chain1_j_gene = "TRAJ33", chain2_v_gene = "(TRBV13|TRBV19)")
print(plt)

S1_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M1")])
S2_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M2","M3","M4")])
MAIT1_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M5")])
MAIT17_ind = names(seurat.MAIT@ident[seurat.MAIT@ident %in% c("M6","M7","M8")])

newcluster = as.vector(sset_MAIT@ident)
names(newcluster) = names(sset_MAIT@ident)
newcluster[names(newcluster) %in% S1_ind] = "S1"
newcluster[names(newcluster) %in% S2_ind] = "S2"
newcluster[names(newcluster) %in% MAIT1_ind] = "MAIT1"
newcluster[names(newcluster) %in% MAIT17_ind] = "MAIT17"

plt = PlotFrequencyByGene(MAIT_cdr3_pair_mat, ident = newcluster, canonical = TRUE,
                          chain1_v_gene = "TRAV1\\*", chain1_j_gene = "TRAJ33", chain2_v_gene = "(TRBV13|TRBV19)")
print(plt)
