#2020-02-10
library(Seurat)
library(xlsx)
library(ggplot2)
# library(Stcr)
load("./tests/testdata/seurat.MAIT.RData") # seurat v2 file
load("./tests/testdata/ensemblGenes2018-10-15.RData")
N1_geneset = read.xlsx("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                       startRow = 5, colIndex = 2)
test = N1_geneset
test = unique(as.vector(test[,1]))
test = test[!is.na(test)]
N1_geneset = test[test != ""]


N2_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 3)
test = N2_geneset
test = unique(as.vector(test[,1]))
N2_geneset = test[test != ""]

N3_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 4)
test = N3_geneset
test = unique(as.vector(test[,1]))
N3_geneset = test[test != ""]

N4_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 5)
test = N4_geneset
test = unique(as.vector(test[,1]))
N4_geneset = test[test != ""]

N5_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 6)
test = N5_geneset
test = unique(as.vector(test[,1]))
N5_geneset = test[test != ""]

N6_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 7)
test = N6_geneset
test = unique(as.vector(test[,1]))
N6_geneset = test[test != ""]

N7_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 8)
test = N7_geneset
test = unique(as.vector(test[,1]))
N7_geneset = test[test != ""]

N1_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N1_geneset]
N2_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N2_geneset]
N3_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N3_geneset]
N4_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N4_geneset]
N5_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N5_geneset]
N6_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N6_geneset]
N7_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N7_geneset]



N1_score =  colMeans(as.matrix(seurat.MAIT@scale.data)[rownames(seurat.MAIT@scale.data) %in% N1_geneset_ensem,])
N2_score =  colMeans(as.matrix(seurat.MAIT@scale.data)[rownames(seurat.MAIT@scale.data) %in% N2_geneset_ensem,])
N3_score =  colMeans(as.matrix(seurat.MAIT@scale.data)[rownames(seurat.MAIT@scale.data) %in% N3_geneset_ensem,])
N4_score =  colMeans(as.matrix(seurat.MAIT@scale.data)[rownames(seurat.MAIT@scale.data) %in% N4_geneset_ensem,])
N5_score =  colMeans(as.matrix(seurat.MAIT@scale.data)[rownames(seurat.MAIT@scale.data) %in% N5_geneset_ensem,])
N6_score =  colMeans(as.matrix(seurat.MAIT@scale.data)[rownames(seurat.MAIT@scale.data) %in% N6_geneset_ensem,])
N7_score =  colMeans(as.matrix(seurat.MAIT@scale.data)[rownames(seurat.MAIT@scale.data) %in% N7_geneset_ensem,])

df = data.frame(x = seurat.MAIT@dr$umap@cell.embeddings[,1],
                y = seurat.MAIT@dr$umap@cell.embeddings[,2],
                N1_score,N2_score,N3_score,N4_score,N5_score,N6_score,N7_score,
                clusters =factor(seurat.MAIT@ident, levels=paste0("M",seq(1,9))))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors= c(gg_color_hue(8)[2],gg_color_hue(8)[8],gg_color_hue(8)[1],gg_color_hue(8)[5],
          gg_color_hue(8)[6],gg_color_hue(8)[7],gg_color_hue(8)[4],gg_color_hue(8)[3])

ggplot(df, aes(x=clusters, y = N1_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  ylim(0,max(df$N1_score)) +
  ylab("N1 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N1_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N2_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  ylim(0,max(df$N2_score)) +
  ylab("N2 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N2_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N3_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  ylim(0,max(df$N3_score)) +
  ylab("N3 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N3_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N4_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  ylim(0,max(df$N4_score)) +
  ylab("N4 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N4_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N5_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  ylim(0,max(df$N5_score)) +
  ylab("N5 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N5_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N6_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  ylim(0,max(df$N6_score)) +
  ylab("N6 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N6_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N7_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  ylim(0,max(df$N7_score)) +
  ylab("N7 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N7_scaled_average_score.pdf", width=9, height=7)



df = data.frame(x=seurat.MAIT@dr$umap@cell.embeddings[,1], 
                y=seurat.MAIT@dr$umap@cell.embeddings[,2], 
                N1_score,N2_score,N3_score,N4_score,N5_score,N6_score,N7_score)
df = df[order(df$N1_score),]
ggplot(df,aes(x=x, y=y, colour=N1_score)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$N1_score)),
                         labels=c(0,round(as.numeric(max(df$N1_score)), digits = 2))) +
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
  ) + ggsave(paste0("MAIT_UMAP_N1_average_score.pdf"), width = 7, height=7)
df = df[order(df$N2_score),]
ggplot(df,aes(x=x, y=y, colour=N2_score)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$N2_score)),
                         labels=c(0,round(as.numeric(max(df$N2_score)), digits = 2))) +
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
  ) + ggsave(paste0("MAIT_UMAP_N2_average_score.pdf"), width = 7, height=7)
df = df[order(df$N3_score),]
ggplot(df,aes(x=x, y=y, colour=N3_score)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$N3_score)),
                         labels=c(0,round(as.numeric(max(df$N3_score)), digits = 2))) +
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
  ) + ggsave(paste0("MAIT_UMAP_N3_average_score.pdf"), width = 7, height=7)
df = df[order(df$N4_score),]
ggplot(df,aes(x=x, y=y, colour=N4_score)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$N4_score)),
                         labels=c(0,round(as.numeric(max(df$N4_score)), digits = 2))) +
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
  ) + ggsave(paste0("MAIT_UMAP_N4_average_score.pdf"), width = 7, height=7)
df = df[order(df$N5_score),]
ggplot(df,aes(x=x, y=y, colour=N5_score)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$N5_score)),
                         labels=c(0,round(as.numeric(max(df$N5_score)), digits = 2))) +
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
  ) + ggsave(paste0("MAIT_UMAP_N5_average_score.pdf"), width = 7, height=7)
df = df[order(df$N6_score),]
ggplot(df,aes(x=x, y=y, colour=N6_score)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$N6_score)),
                         labels=c(0,round(as.numeric(max(df$N6_score)), digits = 2))) +
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
  ) + ggsave(paste0("MAIT_UMAP_N6_average_score.pdf"), width = 7, height=7)
df = df[order(df$N7_score),]
ggplot(df,aes(x=x, y=y, colour=N7_score)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$N7_score)),
                         labels=c(0,round(as.numeric(max(df$N7_score)), digits = 2))) +
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
  ) + ggsave(paste0("MAIT_UMAP_N7_average_score.pdf"), width = 7, height=7)



# scaled
cd = as.matrix(seurat.MAIT@scale.data)
cd = cd[rowSums(cd) > 0,]
zscore = t(scale(t(cd)))
zscore[zscore > 2] = 2
zscore[zscore < -2] = -2
N1_score =  colMeans(zscore[rownames(zscore) %in% N1_geneset_ensem,])
N2_score =  colMeans(zscore[rownames(zscore) %in% N2_geneset_ensem,])
N3_score =  colMeans(zscore[rownames(zscore) %in% N3_geneset_ensem,])
N4_score =  colMeans(zscore[rownames(zscore) %in% N4_geneset_ensem,])
N5_score =  colMeans(zscore[rownames(zscore) %in% N5_geneset_ensem,])
N6_score =  colMeans(zscore[rownames(zscore) %in% N6_geneset_ensem,])
N7_score =  colMeans(zscore[rownames(zscore) %in% N7_geneset_ensem,])

df = data.frame(x = seurat.MAIT@dr$umap@cell.embeddings[,1],
                y = seurat.MAIT@dr$umap@cell.embeddings[,2],
                N1_score,N2_score,N3_score,N4_score,N5_score,N6_score,N7_score,
                clusters =factor(seurat.MAIT@ident, levels=paste0("M",seq(1,9))))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors= c(gg_color_hue(8)[2],gg_color_hue(8)[8],gg_color_hue(8)[1],gg_color_hue(8)[5],
          gg_color_hue(8)[6],gg_color_hue(8)[7],gg_color_hue(8)[4],gg_color_hue(8)[3])

ggplot(df, aes(x=clusters, y = N1_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N1_score)) +
  ylim(-2,2) +
  ylab("N1 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N1_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N2_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N2_score)) +
  ylim(-2,2) +
  ylab("N2 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N2_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N3_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N3_score)) +
  ylim(-2,2) +
  ylab("N3 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N3_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N4_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N4_score)) +
  ylim(-2,2) +
  ylab("N4 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N4_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N5_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N5_score)) +
  ylim(-2,2) +
  ylab("N5 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N5_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N6_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N6_score)) +
  ylim(-2,2) +
  ylab("N6 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N6_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = N7_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N7_score)) +
  ylim(-2,2) +
  ylab("N7 signature score") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
            axis.text.x=element_text(size=20, face="bold"),
            axis.text.y=element_text(size=20, face="bold"),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # axis.line = element_blank(),
            legend.position = "none",
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT_violin_N7_scaled_average_score.pdf", width=9, height=7)

