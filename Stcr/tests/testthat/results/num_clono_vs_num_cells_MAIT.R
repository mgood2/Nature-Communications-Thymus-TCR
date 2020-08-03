ComputeShannonIndex <- function(vector){
  # k = length(vector)
  k = sum(vector > 0)
  return(-sum((vector/sum(vector))*log(vector/sum(vector)), na.rm=TRUE)/log(k))
}


# MAIT
library(Stcr)
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
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


mat_cdr3_pair = MAIT_cdr3_pair_mat[[1]]
df_clonotype = MAIT_cdr3_pair_mat[[2]]
ident = as.factor(seurat.MAIT@ident)
df_mat_cdr3_pair = data.frame(mat_cdr3_pair)
MAITNN = table(colSums(df_mat_cdr3_pair))