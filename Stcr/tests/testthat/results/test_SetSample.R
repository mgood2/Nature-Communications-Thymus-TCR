# init
dataDir = '../testdata/mTCR-S1/outs_IMGT/'
sset = Load10xVDJ(dataDir, filtered = FALSE)

# 1. single sampleID
test_that('set sampleID using single input',{
  sampleID = "F3"
  sset = SetSample(sset,sampleID)
  print(length(sset@contig$cellName))
  expect_equal(length(unique(gsub('(.*)-[ATCG]*','\\1',sset@contig$cellName))),1)
})

# 2. multiple sampleID
test_that('set sampleID using vector input',{
  nF3 = (nrow(sset@contig)-10)
  sampleID = c(rep("F3", nF3),rep("F2",10))
  sset = SetSample(sset,sampleID)
  expect_gt(length(unique(gsub('(.*)-[ATCG]*','\\1',sset@contig$cellName))),1)
  expect_match(sset@contig$cellName[1],"F3")
  expect_match(sset@contig$cellName[length(sset@contig$cellName)],"F2")
})

