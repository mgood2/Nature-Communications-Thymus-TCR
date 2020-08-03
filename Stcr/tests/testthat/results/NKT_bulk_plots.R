setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.NKT.RData") # seurat v2 file
load("../testdata/ident_NKT.RData")
load("../testdata/ensemblGenes2018-10-15.RData")
#TODO: figure1_C adjust contrast grid plot
#bulk markers
logcd = as.matrix(seurat.NKT@scale.data)
zscore_logcd= t(scale(t(logcd)))
zscore_logcd = zscore_logcd[rowSums(zscore_logcd)>0,] 

load(file="../testdata/bulk_NKTp_gdTP.RData")
colnames(NKTp_gdTP)[7] = "Pval"
NKTp_gdTP_FC2_p0.001 = subset(NKTp_gdTP,log2FC > log2(2) & Pval < 0.001)[,"GeneID"]
NKTp_zscore = colMeans(zscore_logcd[rownames(zscore_logcd) %in% NKTp_gdTP_FC2_p0.001,])
load(file="../testdata/bulk_NKT1_gdT1.RData")
colnames(NKT1_gdT1)[7] = "Pval"
NKT1_gdT1_FC2_p0.001 = subset(NKT1_gdT1,log2FC > log2(2) & Pval < 0.001)[,"GeneID"]
NKT1_zscore = colMeans(zscore_logcd[rownames(zscore_logcd) %in% NKT1_gdT1_FC2_p0.001,])
load(file="../testdata/bulk_NKT2_gdT2.RData")
colnames(NKT2_gdT2)[7] = "Pval"
NKT2_gdT2_FC2_p0.001 = subset(NKT2_gdT2,log2FC > log2(2) & Pval < 0.001)[,"GeneID"]
NKT2_zscore = colMeans(zscore_logcd[rownames(zscore_logcd) %in% NKT2_gdT2_FC2_p0.001,])
load(file="../testdata/bulk_NKT17_gdT17.RData")
colnames(NKT17_gdT17)[7] = "Pval"
NKT17_gdT17_FC2_p0.001 = subset(NKT17_gdT17,log2FC > log2(2) & Pval < 0.001)[,"GeneID"]
NKT17_zscore = colMeans(zscore_logcd[rownames(zscore_logcd) %in% NKT17_gdT17_FC2_p0.001,])


load(file="../testdata/bulk_gdTp.RData")
colnames(gdTp)[7] = "Pval"
gdTp_FC2_p0.001 = subset(gdTp,log2FC > log2(2) & Pval < 0.001)[,"GeneID"]
gdTp_zscore = colMeans(zscore_logcd[rownames(zscore_logcd) %in% gdTp_FC2_p0.001,])
load(file="../testdata/bulk_gdT1.RData")
colnames(gdT1)[7] = "Pval"
gdT1_FC2_p0.001 = subset(gdT1,log2FC > log2(2) & Pval < 0.001)[,"GeneID"]
gdT1_zscore = colMeans(zscore_logcd[rownames(zscore_logcd) %in% gdT1_FC2_p0.001,])
load(file="../testdata/bulk_gdT2.RData")
colnames(gdT2)[7] = "Pval"
gdT2_FC2_p0.001 = subset(gdT2,log2FC > log2(2) & Pval < 0.001)[,"GeneID"]
gdT2_zscore = colMeans(zscore_logcd[rownames(zscore_logcd) %in% gdT2_FC2_p0.001,])
load(file="../testdata/bulk_gdT17.RData")
colnames(gdT17)[7] = "Pval"
gdT17_FC2_p0.001 = subset(gdT17,log2FC > log2(2) & Pval < 0.001)[,"GeneID"]
gdT17_zscore = colMeans(zscore_logcd[rownames(zscore_logcd) %in% gdT17_FC2_p0.001,])
gdT17i_markers = read.table(file="../testdata/DEGs_symbols_17i.txt", stringsAsFactors = F)[,1]
gdT17i_markers_ensembl = rownames(ensemblGenes)[ensemblGenes[,"external_gene_name"] %in% gdT17i_markers]
gdT17i_zscore = colMeans(zscore_logcd[rownames(zscore_logcd) %in% gdT17i_markers_ensembl,])

df = data.frame(
  x = seurat.NKT@dr$umap@cell.embeddings[,1],
  y = seurat.NKT@dr$umap@cell.embeddings[,2],
  NKTp_zscore,
  NKT1_zscore,
  NKT2_zscore,
  NKT17_zscore,
  gdTp_zscore,
  gdT1_zscore,
  gdT2_zscore,
  gdT17_zscore,
  gdT17i_zscore)


# not_NKT_ind = rownames(seurat.NKT@meta.data)[seurat.NKT@meta.data$final_celltype != "NKT"]
# df$NKTp_zscore[rownames(df) %in% not_NKT_ind] = NA
df = df[order(df$NKTp_zscore),]
# df = df[rev(rownames(df)),]
p_bulk_NKTp = ggplot(df) + geom_point(aes(x=x,y=y,colour=NKTp_zscore)) +
  scale_colour_gradientn(colors=rev(c("#300000", "red","#eeeeee")))+
  # ylab("iNKT") +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 20),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), 
            axis.line = element_blank(),
            axis.ticks=element_blank(),
            legend.position = "None") + ggsave("NKT_bulk_marker_NKTp.pdf", width=7, height=7)

# not_NKT_ind = rownames(seurat.NKT@meta.data)[seurat.NKT@meta.data$final_celltype != "NKT"]
# df$NKT1_zscore[rownames(df) %in% not_NKT_ind] = NA
df = df[order(df$NKT1_zscore),]
# df = df[rev(rownames(df)),]
p_bulk_NKT1 = ggplot(df) + geom_point(aes(x=x,y=y,colour=NKT1_zscore)) +
  scale_colour_gradientn(colors=rev(c("#300000", "red","#eeeeee")))+
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 20),
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
            axis.ticks=element_blank(),
            legend.position = "None")+ ggsave("NKT_bulk_marker_NKT1.pdf", width=7, height=7)

# not_NKT_ind = rownames(seurat.NKT@meta.data)[seurat.NKT@meta.data$final_celltype != "NKT"]
# df$NKT2_zscore[rownames(df) %in% not_NKT_ind] = NA
df = df[order(df$NKT2_zscore),]
# df = df[rev(rownames(df)),]
p_bulk_NKT2 = ggplot(df) + geom_point(aes(x=x,y=y,colour=NKT2_zscore)) +
  scale_colour_gradientn(colors=rev(c("#300000", "red","#eeeeee")))+
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 20),
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
            axis.ticks=element_blank(),
            legend.position = "None")+ ggsave("NKT_bulk_marker_NKT2.pdf", width=7, height=7)

# not_NKT_ind = rownames(seurat.NKT@meta.data)[seurat.NKT@meta.data$final_celltype != "NKT"]
# df$NKT17_zscore[rownames(df) %in% not_NKT_ind] = NA
df = df[order(df$NKT17_zscore),]
# df = df[rev(rownames(df)),]
p_bulk_NKT17 = ggplot(df) + geom_point(aes(x=x,y=y,colour=NKT17_zscore)) +
  scale_colour_gradientn(colors=rev(c("#300000", "red","#eeeeee")))+
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 20),
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
            axis.ticks=element_blank(),
            legend.position = "None") + ggsave("NKT_bulk_marker_NKT17.pdf", width=7, height=7)