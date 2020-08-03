ComputeShannonIndex <- function(vector){
  # k = length(vector)
  k = sum(vector > 0)
  return(-sum((vector/sum(vector))*log(vector/sum(vector)), na.rm=TRUE)/log(k))
}


# nkT
library(Stcr)
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.NKT.RData") # seurat v2 file
load("../testdata/ident_NKT.RData")
load("../testdata/ensemblGenes2018-10-15.RData")

ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_NKT = ImportSeurat(sset,seurat.NKT, version="v2")

NKT_cdr3_mat = SetCDR3Mat(sset_NKT, cell_ind = colnames(seurat.NKT@scale.data),
                           chain1 = c("TRA"), chain2 = c("TRB"))

NKT_cdr3_pair_mat = SetCDR3PairMat(sset_NKT, NKT_cdr3_mat, chain1 = c("TRA"), chain2 = c("TRB"))


mat_cdr3_pair = NKT_cdr3_pair_mat[[1]]
df_clonotype = NKT_cdr3_pair_mat[[2]]
ident = as.factor(seurat.NKT@ident)
df_mat_cdr3_pair = data.frame(mat_cdr3_pair)


F3_NKT = mat_cdr3_pair[grepl("F3",rownames(mat_cdr3_pair)),]
F5_NKT = mat_cdr3_pair[grepl("F5",rownames(mat_cdr3_pair)),]
ComputeShannonIndex(colSums(F3_NKT)[colSums(F3_NKT)>1]) # 0.9681774
ComputeShannonIndex(colSums(F5_NKT)[colSums(F5_NKT)>1]) # 0.8856559
