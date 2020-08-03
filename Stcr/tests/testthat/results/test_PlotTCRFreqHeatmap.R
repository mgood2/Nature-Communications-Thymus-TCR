
# init
load("../testdata/seurat.gdT.RData") # seurat v2 file

ssetF3 = Load10xVDJ('../testdata/mTCR-F12-gdT/out_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-G1-gdT/out_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5)

sset_gdT = ImportSeurat(sset,seurat.gdT)


test_that('check ploting on gdT total GV vs DV freq log2FC',{
  
  mat_GV = SetVDJcountMat(sset_gdT, chain = "TRG", v_gene = "GV")
  
  mat_DV = SetVDJcountMat(sset_gdT, chain = "TRD", v_gene = "DV")
  
  row_cells = rownames(mat_GV)
  
  mat_GVDV = ((t(mat_GV) %*% mat_DV))/length(unique(row_cells))
  
  i = "G1"
  cluster_cell_ind = names(seurat.gdT@ident)[seurat.gdT@ident == i]
  num_cluster_cell = length(cluster_cell_ind)
  c_mat_GV = mat_GV[rownames(mat_GV) %in% cluster_cell_ind,]
  c_mat_DV = mat_DV[rownames(mat_DV) %in% cluster_cell_ind,]
  c_mat_GVDV = ((t(c_mat_GV) %*% c_mat_DV))/num_cluster_cell
  
  ph_test = PlotTCRFreqHeatmap(c_mat_GVDV,mat_GVDV)

  expect_type(ph_test, "list")
})



