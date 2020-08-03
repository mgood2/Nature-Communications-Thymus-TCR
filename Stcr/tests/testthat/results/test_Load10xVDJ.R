# init
dataDir = '../testdata/mTCR-S1/outs_IMGT/'
sset = Load10xVDJ(dataDir, filtered = FALSE)

# 1. Type checking
test_that('check type of sset',{
  expect_type(sset,'S4')
  expect_match(class(sset)[1],"StcrSet")
})

