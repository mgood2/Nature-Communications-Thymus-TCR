ComputeShannonIndex <- function(vector){
  # k = length(vector)
  k = sum(vector > 0)
  return(-sum((vector/sum(vector))*log(vector/sum(vector)), na.rm=TRUE)/log(k))
}


# gdT
library(Stcr)
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.gdT.RData") # seurat v2 file
load("../testdata/ident_gdT.RData")
load("../testdata/ensemblGenes2018-10-15.RData")

ssetF3 = Load10xVDJ('../testdata/mTCR-F12-gdT/out_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-G1-gdT/out_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")

sset = MergeSamples(ssetF3,ssetF5) 

sset_gdT = ImportSeurat(sset,seurat.gdT, version="v2")

gdT_cdr3_mat = SetCDR3Mat(sset_gdT, cell_ind = colnames(seurat.gdT@scale.data), chain1 = c("TRG"), chain2 = c("TRD","Multi"), v_gene2 = "DV",
                          j_gene2="DJ")



gdT_cdr3_pair_mat = SetCDR3PairMat(sset_gdT, gdT_cdr3_mat, chain1 = c("TRG"), chain2 = c("TRD","Multi"))




mat_cdr3_pair = gdT_cdr3_pair_mat[[1]]
df_clonotype = gdT_cdr3_pair_mat[[2]]
ident = as.factor(seurat.gdT@ident)
df_mat_cdr3_pair = data.frame(mat_cdr3_pair)

gdtNN = table(colSums(df_mat_cdr3_pair))