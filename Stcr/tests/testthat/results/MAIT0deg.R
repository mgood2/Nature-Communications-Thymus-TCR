# init
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.MAIT.RData") # seurat v2 file
library(Seurat)
seurat.MAIT@data = as.matrix(seurat.MAIT@data)
seurat.MAIT@scale.data = as.matrix(seurat.MAIT@scale.data)
seurat.MAIT@raw.data = as.matrix(seurat.MAIT@raw.data)
seuratv3.MAIT = UpdateSeuratObject(seurat.MAIT)
as.loom(seuratv3.MAIT,filename= "MAIT.loom")
load("../testdata/ident_MAIT.RData")
load("../testdata/ensemblGenes2018-10-15.RData")
library(ggplot2)
#nature immu 2019 FRancois
MAIT0DEG=c("Itm2a",
           "Nrgn",
           "Lef1",
           "Bcl2",
           "Slamf6",
           "Izumo1r",
           "Gimap6",
           "Malat1",
           "Id3",
           "Cd28",
           "Ccr9",
           "Cd2",
           "Slc29a1",
           "Gimap5",
           "Satb1",
           "Ccr7",
           "Il21r",
           "Tox",
           "Marcksl1",
           "Adk",
           "Cd27",
           "Cldn10",
           "Cd8a",
           "Tuba1a",
           "Cd4",
           "Actn1",
           "Cd81",
           "Atp1b1",
           "Myb",
           "Cd5",
           "Kcnn4",
           "Ubac2",
           "Sox4",
           "Rasgrp1",
           "Ptprc",
           "Hivep3",
           "Ssbp2",
           "Tspan13",
           "Ldhb",
           "Cd8b1",
           "Tnfrsf9",
           "Ly6d",
           "Cytip",
           "Mrps6",
           "Gsn",
           "Myl10",
           "Rhoh",
           "Themis",
           "Egr2",
           "Nab2",
           "Dntt",
           "Cd69",
           "Cd24a",
           "Gpr83",
           "Asap1",
           "Tmem108",
           "Frat2",
           "Rnf167",
           "Chst2",
           "Tespa1",
           "Spint2",
           "H2-Q2",
           "Fam169b",
           "Egr1",
           "Lztfl1",
           "Tnfsf8",
           "Zfp330",
           "Tubb2b",
           "Cux1",
           "St6gal1",
           "Prkca",
           "Slc25a17",
           "Arpp21",
           "Zcchc12",
           "Mpp4",
           "Ubash3a",
           "Bach2",
           "Traf4",
           "Pip4k2a",
           "Basp1",
           "Patz1",
           "Aqp11",
           "Nr4a1",
           "Lad1",
           "Zfp281",
           "Siae",
           "Mir142hg",
           "Tsc22d1",
           "Vamp1")

MAIT1DEG =c("Bcl2",
            "Cd7",
            "AW112010",
            "Ctla2a",
            "Gimap6",
            "H2-Q7",
            "Nkg7",
            "Ccl5",
            "Ly6c2",
            "Ms4a6b",
            "Klrd1",
            "Gimap4",
            "Klra9",
            "Gzma",
            "Ms4a4b",
            "Cxcr3",
            "Gimap3",
            "Klrk1",
            "Klra3",
            "Xcl1",
            "Klra1",
            "H2-K1",
            "Inpp4b",
            "Hcst",
            "Gm8369",
            "H2-D1",
            "Fgl2",
            "Gimap5",
            "Ctsw",
            "Dusp2",
            "Ctsd",
            "Klre1",
            "Pglyrp1",
            "Gimap1",
            "H2-Q6",
            "Zfp36l2",
            "Fcer1g",
            "Ifit3",
            "Ifitm10",
            "Cd160",
            "Gramd3",
            "Slamf7",
            "Klrc2",
            "Klrb1c",
            "Gm4208",
            "Fasl",
            "Il2rb",
            "Samd3",
            "Itga1",
            "Efhd2",
            "Gimap8",
            "Serpina3g",
            "Gimap7",
            "Klrc1",
            "H2-Q4",
            "Txnip",
            "Tsc22d3",
            "Tbx21",
            "Hsd11b1",
            "Ccl4",
            "Adgre5",
            "Lrrk1",
            "Arsb",
            "Dennd4a",
            "Coro2a",
            "Il10rb",
            "Gm19585")

MAIT17bDEG =c("Lgals3",
              "Tmem176a",
              "Serpinb1a",
              "Pxdc1",
              "Rorc",
              "Tmem176b",
              "Il17re",
              "Fos",
              "Rgs1",
              "Furin",
              "Lingo4",
              "Sema4b",
              "Cd7",
              "S100a4",
              "Capg",
              "Lgals3",
              "Il7r",
              "Il18r1",
              "Tmem176a",
              "Icos",
              "Serpinb1a",
              "Lmo4",
              "S100a11",
              "Pxdc1",
              "Rorc",
              "Tmem176b",
              "Sptssa",
              "Blk",
              "Fam83a",
              "Actn2",
              "Sdc1",
              "Il17re",
              "Selenop",
              "Emb",
              "Fos",
              "S100a6",
              "Ly6g5b",
              "Cxcr6",
              "Ckb",
              "Aqp3",
              "Prr13",
              "Rgs1",
              "Furin",
              "Gm2a",
              "Lgmn",
              "Ramp3",
              "Ccr6",
              "Ccr2",
              "Jaml",
              "Krt83",
              "F2r",
              "Socs3",
              "St3gal6",
              "Id2",
              "Junb",
              "Smox",
              "Lingo4",
              "Maf",
              "Znrf1",
              "Psap",
              "Ypel3",
              "Tnfrsf25",
              "Plin3",
              "Ltb4r1",
              "Ly6a",
              "Ltb",
              "B3galt2",
              "Camk2d",
              "Ero1l",
              "Il1r1",
              "Itgae",
              "Ctss",
              "Prelid2",
              "St6galnac3",
              "Igf1r",
              "Avpi1",
              "Il23r",
              "Cited4",
              "Clstn3",
              "Cd72",
              "Smpdl3a",
              "Cysltr1",
              "5830411N06Rik",
              "Fam129a",
              "Kcnc1",
              "Rnase4",
              "Ahr",
              "Ifi203")
ensemblMAIT0DEG=rownames(ensemblGenes[ensemblGenes$external_gene_name %in% MAIT0DEG,])

ensemblMAIT0DEG= ensemblMAIT0DEG[ensemblMAIT0DEG %in% rownames(seurat.MAIT@data)]

ensemblMAIT1DEG=rownames(ensemblGenes[ensemblGenes$external_gene_name %in% MAIT1DEG,])

ensemblMAIT1DEG= ensemblMAIT1DEG[ensemblMAIT1DEG %in% rownames(seurat.MAIT@data)]

ensemblMAIT17bDEG=rownames(ensemblGenes[ensemblGenes$external_gene_name %in% MAIT17bDEG,])

ensemblMAIT17bDEG= ensemblMAIT17bDEG[ensemblMAIT17bDEG %in% rownames(seurat.MAIT@data)]


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors= c(gg_color_hue(8)[2],gg_color_hue(8)[8],gg_color_hue(8)[1],gg_color_hue(8)[5],
  gg_color_hue(8)[6],gg_color_hue(8)[7],gg_color_hue(8)[4],gg_color_hue(8)[3])



dfMAIT0 = data.frame(x=seurat.MAIT@dr$umap@cell.embeddings[,1],
                y=seurat.MAIT@dr$umap@cell.embeddings[,2],
                exprs=colMeans(as.matrix(seurat.MAIT@data[ensemblMAIT0DEG,])),
                clusters =factor(seurat.MAIT@ident, levels=paste0("M",seq(1,9))))

ggplot(dfMAIT0, aes(x=clusters, y = exprs, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  ylim(0,max(dfMAIT0$exprs)) +
  ylab("MAIT0 signature score") +
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
  ) + ggsave("MAIT0_ref_sig_violin_average_score.pdf", width=9, height=7)
dev.off()
dfMAIT1 = data.frame(x=seurat.MAIT@dr$umap@cell.embeddings[,1],
                     y=seurat.MAIT@dr$umap@cell.embeddings[,2],
                     exprs=colMeans(as.matrix(seurat.MAIT@data[ensemblMAIT1DEG,])),
                     clusters =factor(seurat.MAIT@ident, levels=paste0("M",seq(1,9))))

ggplot(dfMAIT1, aes(x=clusters, y = exprs, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  ) +
  scale_fill_manual(values=colors) +
  ylim(0,max(dfMAIT1$exprs)) +
  ylab("MAIT1 signature score") +
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
  ) + ggsave("MAIT1_ref_sig_violin_average_score.pdf", width=9, height=7)
dev.off()

dfMAIT17b = data.frame(x=seurat.MAIT@dr$umap@cell.embeddings[,1],
                     y=seurat.MAIT@dr$umap@cell.embeddings[,2],
                     exprs=colMeans(as.matrix(seurat.MAIT@data[ensemblMAIT17bDEG,])),
                     clusters =factor(seurat.MAIT@ident, levels=paste0("M",seq(1,9))))

ggplot(dfMAIT17b, aes(x=clusters, y = exprs, fill=clusters)) +
  geom_violin(
    scale ="width",
    adjust = 1,
    trim = TRUE
  )+
  scale_fill_manual(values=colors) +
  ylim(0,max(dfMAIT1$exprs)) +
  ylab("MAIT17b signature score") +
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
            legend.position = "none",
            # axis.line = element_blank(),
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=28, face="bold"),
            # axis.ticks=element_blank()
  ) + ggsave("MAIT17b_ref_sig_violin_average_score.pdf", width=9, height=7)

dev.off()










library(Stcr)



# ggplot(df)+ geom_point(aes(x=x,y=y,colour=exprs)) +
#   scale_colour_gradient(low="grey",high="red") +
#   theme_bw() +
#   theme(    plot.title = element_text(hjust = 0.5,size = 40),
#             axis.text.x=element_blank(),
#             axis.text.y=element_blank(),
#             axis.title.x=element_blank(),
#             axis.title.y=element_blank(),
#             panel.background=element_blank(),
#             panel.border=element_blank(),
#             plot.background=element_blank(),
#             panel.grid.major=element_blank(),
#             panel.grid.minor=element_blank(), 
#             axis.line = element_blank(),
#             axis.ticks=element_blank()
#   ) + 
#   ggsave("NatImmume_MAIT0_deg_average_umap.pdf", width=7,height=7) 


ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_MAIT = ImportSeurat(sset,seurat.MAIT, version="v2")

MAIT_cdr3_mat = SetCDR3Mat(sset_MAIT, cell_ind = colnames(seurat.MAIT@scale.data),
                           chain1 = c("TRA"), chain2 = c("TRB"))

MAIT_cdr3_pair_mat = SetCDR3PairMat(sset_MAIT, MAIT_cdr3_mat, chain1 = c("TRA"), chain2 = c("TRB"))


mat_cdr3_pair = MAIT_cdr3_pair_mat[[1]]
df_clonotype = MAIT_cdr3_pair_mat[[2]]
ident = as.factor(seurat.MAIT@ident)
df_mat_cdr3_pair = data.frame(mat_cdr3_pair)
df_mat_cdr3_pair$ident = ""
for(i in levels(ident)){
  cell_ind = names(ident[ident %in% c(i)])
  df_mat_cdr3_pair$ident[rownames(df_mat_cdr3_pair) %in% cell_ind] = i
}

for(i in 1:nrow(df_clonotype)){
  df_clonotype$ident[i] =  paste( unique(df_mat_cdr3_pair[mat_cdr3_pair[,i] > 0,'ident']), collapse = ",")
}
selected_ind = NULL
chain1_v_gene = "TRAV1\\*"
chain1_j_gene = "TRAJ33"
chain1_d_gene = NULL
chain2_d_gene = NULL
# chain2_v_gene = "(TRBV13|TRBV19)"
chain2_v_gene = NULL # 2nd revision
chain2_j_gene = NULL
df_clonotype2 = subset(df_clonotype, n > 0)
df_clonotype2$selected = FALSE
if(is.null(selected_ind)){
  if(is.null(chain1_v_gene) & is.null(chain1_d_gene) & is.null(chain1_j_gene) & 
     is.null(chain2_v_gene) & is.null(chain2_d_gene) & is.null(chain2_j_gene)){
    stop("neither selected_ind nor gene regex is not available")
  }
  
  combined_index = rep(TRUE,nrow(df_clonotype2))
  if(!is.null(chain1_v_gene)){
    grepl_chain1_v_gene = grepl(chain1_v_gene,df_clonotype2$chain1_v_gene)
    combined_index = combined_index & grepl_chain1_v_gene
  }
  if(!is.null(chain1_d_gene)){
    grepl_chain1_d_gene = grepl(chain1_d_gene,df_clonotype2$chain1_d_gene)
    combined_index = combined_index & grepl_chain1_d_gene
  }
  if(!is.null(chain1_j_gene)){
    grepl_chain1_j_gene = grepl(chain1_j_gene,df_clonotype2$chain1_j_gene)
    combined_index = combined_index & grepl_chain1_j_gene
  }
  if(!is.null(chain2_v_gene)){
    grepl_chain2_v_gene = grepl(chain2_v_gene,df_clonotype2$chain2_v_gene)
    combined_index = combined_index & grepl_chain2_v_gene
  }
  if(!is.null(chain2_d_gene)){
    grepl_chain2_d_gene = grepl(chain2_d_gene,df_clonotype2$chain2_d_gene)
    combined_index = combined_index & grepl_chain2_d_gene
  }
  if(!is.null(chain2_j_gene)){
    grepl_chain2_j_gene = grepl(chain2_j_gene,df_clonotype2$chain2_j_gene)
    combined_index = combined_index & grepl_chain2_j_gene
  }
  selected_ind = rownames(df_clonotype2)[combined_index]
}

df_clonotype2[selected_ind,]$selected = TRUE

df_clonotype_labeled = df_clonotype2
nonselected_clonotypes = as.vector(df_clonotype_labeled$clonotype)[df_clonotype_labeled$selected == FALSE]
selected_clonotypes = as.vector(df_clonotype_labeled$clonotype)[df_clonotype_labeled$selected == TRUE]

nonselected_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,nonselected_clonotypes]) > 0]
selected_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,selected_clonotypes]) > 0]

#Cd24a, Itm2a, Slamf6 and averages noncano cano boxplot

MAIT_box = data.frame(row.names=colnames(seurat.MAIT@data))
MAIT_box$Cd24a = as.matrix(seurat.MAIT@data["ENSMUSG00000047139",])
MAIT_box$Itm2a = as.matrix(seurat.MAIT@data["ENSMUSG00000031239",])
MAIT_box$Slamf6 = as.matrix(seurat.MAIT@data["ENSMUSG00000015314",])
MAIT_box$MAIT0 = colMeans(as.matrix(seurat.MAIT@data[ensemblMAIT0DEG,]))
MAIT_box$MAIT1 = colMeans(as.matrix(seurat.MAIT@data[ensemblMAIT1DEG,]))
MAIT_box$MAIT17b = colMeans(as.matrix(seurat.MAIT@data[ensemblMAIT17bDEG,]))
MAIT_box$cano[rownames(MAIT_box) %in% selected_ind] ="cano."
MAIT_box$cano[rownames(MAIT_box) %in% nonselected_ind] ="non-cano."

MAIT_box= MAIT_box[!is.na(MAIT_box$cano),]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_color_hue(8)[1]

library(grid)
library(gridExtra)


wilcox.test(MAIT0 ~ cano, data = M1_box)
pairwise.wilcox.test(M1_box$MAIT0,M1_box$cano)
M1_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M1"],]
M1_t = wilcox.test(MAIT0 ~ cano, data = M1_box)
# ggplot(M1_box) + geom_boxplot(aes(x=cano, y=Cd24a))
# ggplot(M1_box) + geom_boxplot(aes(x=cano, y=Itm2a))
# ggplot(M1_box) + geom_boxplot(aes(x=cano, y=Slamf6))
M1_plt = ggplot(M1_box) + geom_boxplot(aes(x=cano, y=MAIT0),fill=gg_color_hue(8)[2]) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(0,2) +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=15),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line = element_blank(),
            axis.ticks.x=element_blank()
  ) 

M2_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M2"],]
M2_t = wilcox.test(MAIT0 ~ cano, data = M2_box)
M2_plt = ggplot(M2_box) + geom_boxplot(aes(x=cano, y=MAIT0),fill=gg_color_hue(8)[8]) +
  ylim(0,2) +
  xlab(paste0("M2\n Pvalue = ",format(M2_t$p.value,digits = 2))) +
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
  ) 

M3_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M3"],]
M3_t = wilcox.test(MAIT0 ~ cano, data = M3_box)
M3_plt = ggplot(M1_box) + geom_boxplot(aes(x=cano, y=MAIT0),fill=gg_color_hue(8)[1]) +
  xlab(paste0("M3\n Pvalue = ",format(M3_t$p.value,digits = 2))) +
  ylim(0,2) +
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
  ) 


M4_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M4"],]
M4_t = wilcox.test(MAIT0 ~ cano, data = M4_box)
M4_plt = ggplot(M4_box) + geom_boxplot(aes(x=cano, y=MAIT0),fill=gg_color_hue(8)[5]) +
  xlab(paste0("M4\n Pvalue = ",format(M4_t$p.value,digits = 2))) +
  ylim(0,2) +
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
  ) 


M5_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M5"],]
M5_t = wilcox.test(MAIT0 ~ cano, data = M5_box)
M5_plt = ggplot(M5_box) + geom_boxplot(aes(x=cano, y=MAIT0),fill=gg_color_hue(8)[6]) +
  xlab(paste0("M5\n Pvalue = ",format(M5_t$p.value,digits = 2))) +
  ylim(0,2) +
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
  ) 

M6_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M6"],]
M6_t = wilcox.test(MAIT0 ~ cano, data = M6_box)
M6_plt = ggplot(M6_box) + geom_boxplot(aes(x=cano, y=MAIT0),fill=gg_color_hue(8)[7]) +
  xlab(paste0("M6\n Pvalue = ",format(M6_t$p.value,digits = 2))) +
  ylim(0,2) +
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
  ) 

M7_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M7"],]
M7_t = wilcox.test(MAIT0 ~ cano, data = M7_box)
M7_plt = ggplot(M7_box) + geom_boxplot(aes(x=cano, y=MAIT0),fill=gg_color_hue(8)[4]) +
  xlab(paste0("M7\n Pvalue = ",format(M7_t$p.value,digits = 2))) +
  ylim(0,2) +
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
  ) 

M8_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M8"],]
M8_t = wilcox.test(MAIT0 ~ cano, data = M8_box)
M8_plt = ggplot(M8_box) + geom_boxplot(aes(x=cano, y=MAIT0),fill=gg_color_hue(8)[3]) +
  xlab(paste0("M8\n Pvalue = ",format(M8_t$p.value,digits = 2))) +
  ylim(0,2) +
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
  ) 

pdf("MAIT0_marker_cano_noncano_comparison_wilcoxtest_alpha_beta.pdf",width=12,height =7)
grid.arrange(M1_plt,M2_plt,M3_plt,M4_plt,M5_plt,M6_plt,M7_plt,M8_plt, ncol=8)
dev.off()



allinone_plt = ggplot(MAIT_box) + geom_boxplot(aes(x=cano, y=MAIT0,fill=cano)) +
  # xlab(paste0("M8\n Pvalue = ",format(M8_t$p.value,digits = 2))) +
  ylim(0,2) +
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
  ) 

M1_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M1"],]
M1_t = wilcox.test(MAIT0 ~ cano, data = M1_box)
M1_plt = ggplot(M1_box) + geom_boxplot(aes(x=cano, y=MAIT0,fill=cano)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(0,2) +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=15),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line = element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "none"
  ) 

M2_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M2"],]
M2_t = wilcox.test(MAIT0 ~ cano, data = M2_box)
M2_plt = ggplot(M2_box) + geom_boxplot(aes(x=cano, y=MAIT0,fill=cano)) +
  ylim(0,2) +
  xlab(paste0("M2\n Pvalue = ",format(M2_t$p.value,digits = 2))) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M3_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M3"],]
M3_t = wilcox.test(MAIT0 ~ cano, data = M3_box)
M3_plt = ggplot(M1_box) + geom_boxplot(aes(x=cano, y=MAIT0,fill=cano)) +
  xlab(paste0("M3\n Pvalue = ",format(M3_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 


M4_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M4"],]
M4_t = wilcox.test(MAIT0 ~ cano, data = M4_box)
M4_plt = ggplot(M4_box) + geom_boxplot(aes(x=cano, y=MAIT0,fill=cano)) +
  xlab(paste0("M4\n Pvalue = ",format(M4_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 


M5_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M5"],]
M5_t = wilcox.test(MAIT0 ~ cano, data = M5_box)
M5_plt = ggplot(M5_box) + geom_boxplot(aes(x=cano, y=MAIT0,fill=cano)) +
  xlab(paste0("M5\n Pvalue = ",format(M5_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M6_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M6"],]
M6_t = wilcox.test(MAIT0 ~ cano, data = M6_box)
M6_plt = ggplot(M6_box) + geom_boxplot(aes(x=cano, y=MAIT0,fill=cano)) +
  xlab(paste0("M6\n Pvalue = ",format(M6_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M7_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M7"],]
M7_t = wilcox.test(MAIT0 ~ cano, data = M7_box)
M7_plt = ggplot(M7_box) + geom_boxplot(aes(x=cano, y=MAIT0,fill=cano)) +
  xlab(paste0("M7\n Pvalue = ",format(M7_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M8_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M8"],]
M8_t = wilcox.test(MAIT0 ~ cano, data = M8_box)
M8_plt = ggplot(M8_box) + geom_boxplot(aes(x=cano, y=MAIT0,fill=cano)) +
  xlab(paste0("M8\n Pvalue = ",format(M8_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

pdf("MAIT0_marker_cano_noncano_comparison_cano_alpha.pdf",width=12,height =7)
grid.arrange(M1_plt,M2_plt,M3_plt,M4_plt,M5_plt,M6_plt,M7_plt,M8_plt, ncol=8)
dev.off()


M1_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M1"],]
M1_t = wilcox.test(MAIT1 ~ cano, data = M1_box)
M1_plt = ggplot(M1_box) + geom_boxplot(aes(x=cano, y=MAIT1,fill=cano)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(0,2) +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=15),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line = element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "none"
  ) 

M2_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M2"],]
M2_t = wilcox.test(MAIT1 ~ cano, data = M2_box)
M2_plt = ggplot(M2_box) + geom_boxplot(aes(x=cano, y=MAIT1,fill=cano)) +
  ylim(0,2) +
  xlab(paste0("M2\n Pvalue = ",format(M2_t$p.value,digits = 2))) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M3_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M3"],]
M3_t = wilcox.test(MAIT1 ~ cano, data = M3_box)
M3_plt = ggplot(M1_box) + geom_boxplot(aes(x=cano, y=MAIT1,fill=cano)) +
  xlab(paste0("M3\n Pvalue = ",format(M3_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 


M4_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M4"],]
M4_t = wilcox.test(MAIT1 ~ cano, data = M4_box)
M4_plt = ggplot(M4_box) + geom_boxplot(aes(x=cano, y=MAIT1,fill=cano)) +
  xlab(paste0("M4\n Pvalue = ",format(M4_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 


M5_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M5"],]
M5_t = wilcox.test(MAIT1 ~ cano, data = M5_box)
M5_plt = ggplot(M5_box) + geom_boxplot(aes(x=cano, y=MAIT1,fill=cano)) +
  xlab(paste0("M5\n Pvalue = ",format(M5_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M6_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M6"],]
M6_t = wilcox.test(MAIT1 ~ cano, data = M6_box)
M6_plt = ggplot(M6_box) + geom_boxplot(aes(x=cano, y=MAIT1,fill=cano)) +
  xlab(paste0("M6\n Pvalue = ",format(M6_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M7_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M7"],]
M7_t = wilcox.test(MAIT1 ~ cano, data = M7_box)
M7_plt = ggplot(M7_box) + geom_boxplot(aes(x=cano, y=MAIT1,fill=cano)) +
  xlab(paste0("M7\n Pvalue = ",format(M7_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M8_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M8"],]
M8_t = wilcox.test(MAIT1 ~ cano, data = M8_box)
M8_plt = ggplot(M8_box) + geom_boxplot(aes(x=cano, y=MAIT1,fill=cano)) +
  xlab(paste0("M8\n Pvalue = ",format(M8_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

pdf("MAIT1_marker_cano_noncano_comparison_cano_alpha.pdf",width=12,height =7)
grid.arrange(M1_plt,M2_plt,M3_plt,M4_plt,M5_plt,M6_plt,M7_plt,M8_plt, ncol=8)
dev.off()



M1_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M1"],]
M1_t = wilcox.test(MAIT17b ~ cano, data = M1_box)
M1_plt = ggplot(M1_box) + geom_boxplot(aes(x=cano, y=MAIT17b,fill=cano)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(0,2) +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=15),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line = element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "none"
  ) 

M2_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M2"],]
M2_t = wilcox.test(MAIT17b ~ cano, data = M2_box)
M2_plt = ggplot(M2_box) + geom_boxplot(aes(x=cano, y=MAIT17b,fill=cano)) +
  ylim(0,2) +
  xlab(paste0("M2\n Pvalue = ",format(M2_t$p.value,digits = 2))) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M3_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M3"],]
M3_t = wilcox.test(MAIT17b ~ cano, data = M3_box)
M3_plt = ggplot(M1_box) + geom_boxplot(aes(x=cano, y=MAIT17b,fill=cano)) +
  xlab(paste0("M3\n Pvalue = ",format(M3_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 


M4_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M4"],]
M4_t = wilcox.test(MAIT17b ~ cano, data = M4_box)
M4_plt = ggplot(M4_box) + geom_boxplot(aes(x=cano, y=MAIT17b,fill=cano)) +
  xlab(paste0("M4\n Pvalue = ",format(M4_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 


M5_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M5"],]
M5_t = wilcox.test(MAIT17b ~ cano, data = M5_box)
M5_plt = ggplot(M5_box) + geom_boxplot(aes(x=cano, y=MAIT17b,fill=cano)) +
  xlab(paste0("M5\n Pvalue = ",format(M5_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M6_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M6"],]
M6_t = wilcox.test(MAIT17b ~ cano, data = M6_box)
M6_plt = ggplot(M6_box) + geom_boxplot(aes(x=cano, y=MAIT17b,fill=cano)) +
  xlab(paste0("M6\n Pvalue = ",format(M6_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M7_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M7"],]
M7_t = wilcox.test(MAIT17b ~ cano, data = M7_box)
M7_plt = ggplot(M7_box) + geom_boxplot(aes(x=cano, y=MAIT17b,fill=cano)) +
  xlab(paste0("M7\n Pvalue = ",format(M7_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

M8_box = MAIT_box[rownames(MAIT_box) %in% names(seurat.MAIT@ident)[seurat.MAIT@ident=="M8"],]
M8_t = wilcox.test(MAIT17b ~ cano, data = M8_box)
M8_plt = ggplot(M8_box) + geom_boxplot(aes(x=cano, y=MAIT17b,fill=cano)) +
  xlab(paste0("M8\n Pvalue = ",format(M8_t$p.value,digits = 2))) +
  ylim(0,2) +
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
            axis.ticks=element_blank(),
            legend.position = "none"
  ) 

pdf("MAIT17b_marker_cano_noncano_comparison_cano_alpha.pdf",width=12,height =7)
grid.arrange(M1_plt,M2_plt,M3_plt,M4_plt,M5_plt,M6_plt,M7_plt,M8_plt, ncol=8)
dev.off()
