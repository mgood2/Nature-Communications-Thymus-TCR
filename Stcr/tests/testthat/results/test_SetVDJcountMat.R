# init
load("../testdata/seurat.NKT.RData") # seurat v2 file

ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")
sset = MergeSamples(ssetF3,ssetF5)
sset_NKT = ImportSeurat(sset,seurat.NKT)

test_that('check SetVDJcountMat on NKT data total TRA v_gene counts',{

  mat_AV = SetVDJcountMat(sset_NKT, chain = "TRA", v_gene = "AV")
  expect_equal(sum(rowSums(mat_AV) > 0),2833)
})
