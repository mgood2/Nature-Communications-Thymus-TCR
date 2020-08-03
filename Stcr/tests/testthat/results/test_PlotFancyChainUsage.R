# init
setwd("D:/elee4472/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.gdT.RData") # seurat v2 file
load("../testdata/ident_gdT.RData")
ssetF3 = Load10xVDJ('../testdata/mTCR-F12-gdT/out_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-G1-gdT/out_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5)

sset_gdT = ImportSeurat(sset,seurat.gdT)

gdT_cdr3_mat = SetCDR3Mat(sset_gdT, cell.ind = colnames(seurat.gdT@scale.data), chain1 = c("TRG"), chain2 = c("TRD","Multi"), v_gene2 = "DV",
                          j_gene2="DJ")

gdT_cdr3_pair_mat = SetCDR3PairMat(sset_gdT, gdT_cdr3_mat, chain1 = c("TRG"), chain2 = c("TRD","Multi"))

plt = PlotFancyChainUsage(gdT_cdr3_mat, cellLimit = 4)
print(plt)
# viz using circlize
library(circlize)
library(RColorBrewer)
mat = gdT_cdr3_mat[[1]]
mat = mat[rowSums(mat) > 2,colSums(mat) > 2]

pdf("test.pdf")

circos.par(gap.degree = c(rep(1, nrow(mat)-1), 10, rep(1, ncol(mat)-1), 15), start.degree = 5)

rcols <- rep(brewer.pal(12, "Paired"), nrow(mat)/12 + 1)[1:nrow(mat)]
ccols <- rep(brewer.pal(12, "Paired"), ncol(mat)/12 + 1)[1:ncol(mat)]

names(rcols) <- sort(rownames(mat))
names(ccols) <- sort(colnames(mat))

chordDiagram(mat, annotationTrack = "grid",
             grid.col = c(rcols, ccols),
             preAllocateTracks = list(track.height = 0.2), transparency = 0.5)

circos.trackPlotRegion(track.index = 1, bg.border = NA,
                       panel.fun = function(x, y) {
                         sector.name = get.cell.meta.data("sector.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                       }
)

circos.clear()

dev.off()

