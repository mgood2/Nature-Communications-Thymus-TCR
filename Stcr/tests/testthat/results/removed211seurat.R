library(Seurat) #(v2)
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
load(file="cellind_removed_221cells.RData")
load("../testdata/ensemblGenes2018-10-15.RData")

load("../testdata/seurat.MAIT.RData") # seurat v2 file
# load("../testdata/ident_MAIT.RData")
load("../testdata/ensemblGenes2018-10-15.RData")

seurat.MAIT.removed221 = SubsetData(seurat.MAIT, cells =  colnames(seurat.MAIT@scale.data)[(colnames(seurat.MAIT@scale.data) %in% cellind_removed_221cells)])

library(scran)

sce_MAIT = as.SingleCellExperiment(seurat.MAIT.removed221)

var.fit <- trendVar(sce_MAIT, method="spline", parametric=TRUE, 
                    use.spikes=FALSE)#, design = batchDesign)
var.out <- decomposeVar(sce_MAIT, var.fit)
hvg_MAIT <- var.out[which(var.out$FDR < 0.05 & var.out$bio > .2),]


seurat.MAIT.removed221@var.genes = rownames(hvg_MAIT)
PCA = 50
seurat.MAIT.removed221 <- RunPCA(seurat.MAIT.removed221, pcs.compute = PCA, weight.by.var = FALSE)
qplot(x = seq(1:PCA),y = seurat.MAIT.removed221@dr$pca@sdev,
      xlab = "PC", ylab = "Eigenvalue")


PCAuse = 25
seurat.MAIT.removed221 <- RunTSNE(seurat.MAIT.removed221, dims.use = 1:PCAuse, do.fast = T, seed.use = 42, perplexity=25)
seurat.MAIT.removed221 <- FindClusters(seurat.MAIT.removed221, reduction.type="pca", dims.use = 1:PCAuse, save.SNN = TRUE, force.recalc = TRUE)
seurat.MAIT.removed221 <- RunUMAP(seurat.MAIT.removed221, dims.use = 1:PCAuse)
save(seurat.MAIT.removed221, file="seurat.MAIT.removed221.RData")
load(file="seurat.MAIT.removed221.RData")

DimPlot(seurat.MAIT.removed221, reduction.use = "umap") +
  ylab("UMAP2") +
  xlab("UMAP1") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), 
            axis.line = element_blank(),
            axis.ticks=element_blank()
  ) + ggsave("UMAP_seurat.MAIT.removed221.pdf", height=7, width=7)
  
genes = c("Cd24a", "Cd44", "Plzf", "Gata3", "Cd4", "Cd8a", "Rorc", "Ccr6", "Tbx21", "Il2rb", "Ifng", "Cxcr3")
genes_symbols = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% genes]
for( i in genes_symbols){
  df = data.frame(x=seurat.MAIT.removed221@dr$umap@cell.embeddings[,1], 
                  y=seurat.MAIT.removed221@dr$umap@cell.embeddings[,2], 
                  expression = seurat.MAIT.removed221@scale.data[i,])
  df = df[order(df$expression),]
  ggplot(df,aes(x=x, y=y, colour=expression)) +
    geom_point() + 
    scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                           breaks=c(0,max(df$expression)),
                           labels=c(0,round(as.numeric(max(df$expression)), digits = 2))) +
    ylab("UMAP2") +
    xlab("UMAP1") +
    theme_bw() +
    theme(    plot.title = element_text(hjust = 0.5,size = 40),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_blank(),
              panel.border=element_blank(),
              plot.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(), 
              axis.line = element_blank(),
              axis.ticks=element_blank()
    ) + ggsave(paste0("UMAP_seurat.MAIT.removed221_",ensemblGenes[i,"external_gene_name"],".pdf"), width = 7, height=7)
}

