#2020-02-10
library(Seurat)
library(xlsx)
library(ggplot2)
# library(Stcr)
load("./tests/testdata/seurat.NKT.RData") # seurat v2 file
load("./tests/testdata/ensemblGenes2018-10-15.RData")
M1_geneset = read.xlsx("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                       startRow = 5, colIndex = 9)
test = M1_geneset
test = unique(as.vector(test[,1]))
test = test[!is.na(test)]
M1_geneset = test[test != ""]


M2_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 10)
test = M2_geneset
test = unique(as.vector(test[,1]))
M2_geneset = test[test != ""]

M3_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 11)
test = M3_geneset
test = unique(as.vector(test[,1]))
M3_geneset = test[test != ""]

M4_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 12)
test = M4_geneset
test = unique(as.vector(test[,1]))
M4_geneset = test[test != ""]

M5_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 13)
test = M5_geneset
test = unique(as.vector(test[,1]))
M5_geneset = test[test != ""]

M6_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 14)
test = M6_geneset
test = unique(as.vector(test[,1]))
M6_geneset = test[test != ""]

M7_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 15)
test = M7_geneset
test = unique(as.vector(test[,1]))
M7_geneset = test[test != ""]

M8_geneset = read.xlsx2("./tests/testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 16)
test = M8_geneset
test = unique(as.vector(test[,1]))
M8_geneset = test[test != ""]

M1_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% M1_geneset]
M2_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% M2_geneset]
M3_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% M3_geneset]
M4_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% M4_geneset]
M5_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% M5_geneset]
M6_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% M6_geneset]
M7_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% M7_geneset]
M8_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% M8_geneset]




# scaled
cd = as.matrix(seurat.NKT@scale.data)
cd = cd[rowSums(cd) > 0,]
zscore = t(scale(t(cd)))
zscore[zscore > 2] = 2
zscore[zscore < -2] = -2
M1_score =  colMeans(zscore[rownames(zscore) %in% M1_geneset_ensem,])
M2_score =  colMeans(zscore[rownames(zscore) %in% M2_geneset_ensem,])
M3_score =  colMeans(zscore[rownames(zscore) %in% M3_geneset_ensem,])
M4_score =  colMeans(zscore[rownames(zscore) %in% M4_geneset_ensem,])
M5_score =  colMeans(zscore[rownames(zscore) %in% M5_geneset_ensem,])
M6_score =  colMeans(zscore[rownames(zscore) %in% M6_geneset_ensem,])
M7_score =  colMeans(zscore[rownames(zscore) %in% M7_geneset_ensem,])
M8_score =  colMeans(zscore[rownames(zscore) %in% M8_geneset_ensem,])

df = data.frame(x = seurat.NKT@dr$umap@cell.embeddings[,1],
                y = seurat.NKT@dr$umap@cell.embeddings[,2],
                M1_score,M2_score,M3_score,M4_score,M5_score,M6_score,M7_score,M8_score,
                clusters =factor(seurat.NKT@ident, levels=paste0("N",seq(1,7))))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors= c(gg_color_hue(7)[5],gg_color_hue(7)[6],gg_color_hue(7)[2],gg_color_hue(7)[3],
          gg_color_hue(7)[1],gg_color_hue(7)[7],gg_color_hue(7)[4])

ggplot(df, aes(x=clusters, y = M1_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$M1_score)) +
  ylim(-2,2) +
  ylab("M1 signature score") +
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
  ) + ggsave("NKT_violin_M1_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = M2_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$M2_score)) +
  ylim(-2,2) +
  ylab("M2 signature score") +
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
  ) + ggsave("NKT_violin_M2_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = M3_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N3_score)) +
  ylim(-2,2) +
  ylab("M3 signature score") +
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
  ) + ggsave("NKT_violin_M3_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = M4_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N4_score)) +
  ylim(-2,2) +
  ylab("M4 signature score") +
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
  ) + ggsave("NKT_violin_M4_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = M5_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N5_score)) +
  ylim(-2,2) +
  ylab("M5 signature score") +
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
  ) + ggsave("NKT_violin_M5_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = M6_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N6_score)) +
  ylim(-2,2) +
  ylab("M6 signature score") +
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
  ) + ggsave("NKT_violin_M6_scaled_average_score.pdf", width=9, height=7)
ggplot(df, aes(x=clusters, y = M7_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N7_score)) +
  ylim(-2,2) +
  ylab("M7 signature score") +
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
  ) + ggsave("NKT_violin_M7_scaled_average_score.pdf", width=9, height=7)

ggplot(df, aes(x=clusters, y = M8_score, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  # ylim(0,max(df$N7_score)) +
  ylim(-2,2) +
  ylab("M8 signature score") +
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
  ) + ggsave("NKT_violin_M8_scaled_average_score.pdf", width=9, height=7)
