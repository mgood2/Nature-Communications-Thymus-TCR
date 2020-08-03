# init
library(Stcr)
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.MAIT.RData") # seurat v2 file
load("../testdata/ident_MAIT.RData")
load("../testdata/ensemblGenes2018-10-15.RData")

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


df = data.frame(x = seurat.MAIT@dr$umap@cell.embeddings[,1],
                y = seurat.MAIT@dr$umap@cell.embeddings[,2],
                cano = colnames(seurat.MAIT@data) %in% selected_ind,
                noncano = colnames(seurat.MAIT@data) %in% nonselected_ind,
                TCR = colnames(seurat.MAIT@data) %in% selected_ind + colnames(seurat.MAIT@data) %in% nonselected_ind,
                cluster = seurat.MAIT@ident)


ggdf = as.data.frame(table(subset(df,noncano == TRUE)$cluster)/table((subset(df,TCR == 1)$cluster)))
library(ggplot2)
library(ggrepel)
plt1 = PlotFrequency(MAIT_cdr3_pair_mat, ident = ident_MAIT, reverse = TRUE, selected_ind = NULL,
                    chain1_v_gene = "TRAV1\\*", chain1_d_gene = NULL, chain1_j_gene = NULL, 
                    chain2_v_gene = NULL, chain2_d_gene = NULL, chain2_j_gene = NULL)

plt2 = PlotFrequency(MAIT_cdr3_pair_mat, ident = ident_MAIT, reverse = TRUE, selected_ind = NULL,
                    chain1_v_gene = NULL, chain1_d_gene = NULL, chain1_j_gene = "TRAJ33", 
                    chain2_v_gene = NULL, chain2_d_gene = NULL, chain2_j_gene = NULL)

plt3 = PlotFrequency(MAIT_cdr3_pair_mat, ident = ident_MAIT, reverse = TRUE, selected_ind = NULL,
                    chain1_v_gene = NULL, chain1_d_gene = NULL, chain1_j_gene = NULL, 
                    chain2_v_gene = "(TRBV13|TRBV19)", chain2_d_gene = NULL, chain2_j_gene = NULL)

plt4 = PlotFrequency(MAIT_cdr3_pair_mat, ident = ident_MAIT, reverse = TRUE, selected_ind = NULL,
                     chain1_v_gene = "TRAV1\\*", chain1_d_gene = NULL, chain1_j_gene = "TRAJ33", 
                     chain2_v_gene = NULL, chain2_d_gene = NULL, chain2_j_gene = NULL)

plt1 = plt1 + theme_bw() + ylab("") +
  geom_point(size=4) +
  theme(
        panel.background=element_blank(),
        panel.border=element_blank(),
        plot.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        
        axis.line = element_line(),
        axis.title=element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x =element_text(colour=c("black","darkred","purple", "lightgreen","lightgreen","blue","blue","blue"), angle = 90, vjust = 0.5),
        legend.position = "none") + ggsave("MAIT_noncano_freq_TRAV1.pdf",width=7,height=7)
dev.off()
plt2 = plt2 + theme_bw() + ylab("") +
  geom_point(size=4) +
  theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    plot.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), 
    axis.line = element_line(),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20),
    axis.text.x =element_text(colour=c("black","darkred","purple", "lightgreen","lightgreen","blue","blue","blue"), angle = 90, vjust = 0.5),
    legend.position = "none") + ggsave("MAIT_noncano_freq_TRAJ33.pdf",width=7,height=7)
dev.off()
plt3 = plt3 + theme_bw() + ylab("") +
  geom_point(size=4) +
  theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    plot.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), 
    axis.line = element_line(),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20),
    axis.text.x =element_text(colour=c("black","darkred","purple", "lightgreen","lightgreen","blue","blue","blue"), angle = 90, vjust = 0.5),
    legend.position = "none") + ggsave("MAIT_noncano_freq_TRBV13orTRBV19.pdf",width=7,height=7)
dev.off()

plt4 = plt4 + theme_bw() + ylab("") +
  geom_point(size=5) +
  geom_line(group = 0,size=1.5) +
  ylim(0,25) +
  theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    plot.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), 
    axis.line = element_line(),
    axis.title=element_text(size=20),
    axis.text=element_text(size=30),
    axis.text.x =element_text(colour=c("black","darkred","purple", "lightgreen","lightgreen","blue","blue","blue"), angle = 90, vjust = 0.5),
    legend.position = "none") + ggsave("MAIT_noncano_freq_alphachain.pdf",width=7,height=7)
dev.off()
print(plt2)

# write.csv(tmp, file="tmp.csv")
table(subset(df, cano == TRUE)$cluster)/sum(table(subset(df, cano == TRUE)$cluster))
table(subset(df, cano == FALSE)$cluster)/sum(table(subset(df, cano == FALSE)$cluster))
library(Seurat)
seurat.MAIT@data = as.matrix(seurat.MAIT@data)
seurat.MAIT@scale.data = as.matrix(seurat.MAIT@scale.data)
updated_seurat_MAIT = UpdateSeuratObject(seurat.MAIT)
updated_seurat_MAIT@meta.data$cano = colnames(seurat.MAIT@data) %in% selected_ind

deg_difference = FindMarkers(updated_seurat_MAIT, ident.1="TRUE", ident.2="FALSE", group.by="cano")
deg_difference$symbol = ensemblGenes[rownames(deg_difference),"external_gene_name"]

library(ggplot2)

M1_deg =  FindMarkers(updated_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M1", group.by="cano", logfc.threshold = 0)
M1_deg$symbol = ensemblGenes[rownames(M1_deg),"external_gene_name"]
M1_deg$sig = M1_deg$p_val_adj < 0.05 & abs(M1_deg$avg_logFC) > 1
table(M1_deg$sig)
M1_deg$symbol[!M1_deg$sig] = ""
M1_plt = ggplot(M1_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text_repel(size=10) +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M1") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(
        plot.title = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        # axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+ ggsave("M1_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M2_deg =  FindMarkers(updated_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M2", group.by="cano", logfc.threshold = 0)
M2_deg$symbol = ensemblGenes[rownames(M2_deg),"external_gene_name"]
M2_deg$sig = M2_deg$p_val_adj < 0.05 & abs(M2_deg$avg_logFC) > 1
table(M2_deg$sig)
M2_deg$symbol[!M2_deg$sig] = ""
M2_plt = ggplot(M2_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text(size=10) +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M2") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(        plot.title = element_blank(),
                axis.title = element_blank(),
                line = element_blank(),
        # axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M2_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M3_deg =  FindMarkers(updated_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M3", group.by="cano", logfc.threshold = 0)
M3_deg$symbol = ensemblGenes[rownames(M3_deg),"external_gene_name"]
M3_deg$sig = M3_deg$p_val_adj < 0.05 & abs(M3_deg$avg_logFC) > 1
table(M3_deg$sig)
M3_deg$symbol[!M3_deg$sig] = ""
M3_plt = ggplot(M3_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text_repel(size=10) +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant")  +
  ggtitle("M3") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(        plot.title = element_blank(),
                axis.title = element_blank(),
                line = element_blank(),
        # axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M3_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M4_deg =  FindMarkers(updated_seurat_MAIT, ident.1="TRUE", ident.2="FALSE", 
                      subset.ident = "M4", group.by="cano", logfc.threshold = 0)
M4_deg$symbol = ensemblGenes[rownames(M4_deg),"external_gene_name"]
M4_deg$sig = M4_deg$p_val_adj < 0.05 & abs(M4_deg$avg_logFC) > 1
table(M4_deg$sig)
M4_deg$symbol[!(M4_deg$sig)] = ""
M4_plt = ggplot(M4_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text(size=10) +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant")  +
  ggtitle("M4") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(        plot.title = element_blank(),
                axis.title = element_blank(),
                line = element_blank(),
        # axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M4_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M5_deg =  FindMarkers(updated_seurat_MAIT, ident.1="TRUE", ident.2="FALSE", 
                      subset.ident = "M5", group.by="cano", logfc.threshold = 0)
M5_deg$symbol = ensemblGenes[rownames(M5_deg),"external_gene_name"]
M5_deg$sig = M5_deg$p_val_adj < 0.05 & abs(M5_deg$avg_logFC) > 1
table(M5_deg$sig)
M5_deg$symbol[!M5_deg$sig] = ""
M5_plt = ggplot(M5_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text(size=10) +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant")  +
  ggtitle("M5") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(        plot.title = element_blank(),
                axis.title = element_blank(),
                line = element_blank(),
        # axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M5_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M6_deg =  FindMarkers(updated_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M6", group.by="cano", logfc.threshold = 0)
M6_deg$symbol = ensemblGenes[rownames(M6_deg),"external_gene_name"]
M6_deg$sig = M6_deg$p_val_adj < 0.05 & abs(M6_deg$avg_logFC) > 1
table(M6_deg$sig)
M6_deg$symbol[!M6_deg$sig] = ""
M6_plt = ggplot(M6_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text(size=10) +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M6") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(        plot.title = element_blank(),
                axis.title = element_blank(),
                line = element_blank(),
        # axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M6_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M7_deg =  FindMarkers(updated_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M7", group.by="cano", logfc.threshold = 0)
M7_deg$symbol = ensemblGenes[rownames(M7_deg),"external_gene_name"]
M7_deg$sig = M7_deg$p_val_adj < 0.05 & abs(M7_deg$avg_logFC) > 1
table(M7_deg$sig)
M7_deg$symbol[!M7_deg$sig] = ""
M7_plt = ggplot(M7_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text(size=10) +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M7") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(        plot.title = element_blank(),
                axis.title = element_blank(),
    line = element_blank(),
        # axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none")# +  ggsave("M7_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M8_deg =  FindMarkers(updated_seurat_MAIT, ident.1="TRUE", ident.2="FALSE", 
                      subset.ident = "M8", group.by="cano", logfc.threshold = 0)
M8_deg$symbol = ensemblGenes[rownames(M8_deg),"external_gene_name"]
M8_deg$sig = M8_deg$p_val_adj < 0.05 & abs(M8_deg$avg_logFC) > 1
table(M8_deg$sig)
M8_deg$symbol[!M8_deg$sig] = ""
M8_plt = ggplot(M8_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text(size=10) +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M8") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(        plot.title = element_blank(),
                axis.title = element_blank(),
                line = element_blank(),
        # axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M8_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)


# ggplot(M2_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig, label= symbol)) + geom_point() + geom_text()
# ggplot(M3_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig, label= symbol)) + geom_point() + geom_text()
# ggplot(M4_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig, label= symbol)) + geom_point() + geom_text()
# ggplot(M5_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig, label= symbol)) + geom_point() + geom_text()
# ggplot(M6_deg, aes(x=avg_logFC, y= -log10(p_val_adj), label= symbol)) + geom_point() + geom_text()
# ggplot(M7_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig, label= symbol)) + geom_point() + geom_text()
# ggplot(M8_deg, aes(x=avg_logFC, y= -log10(p_val_adj), label= symbol)) + geom_point() + geom_text()

library(grid)
library(gridExtra)



pdf("MAIT_deg_alpha_beta_cano_noncano_volcano_alpha_names.pdf",width=16,height =8)
grid.arrange(M1_plt,M2_plt,M3_plt,M4_plt,M5_plt,M6_plt,M7_plt,M8_plt, ncol=4)
dev.off()



# DP comparison
load("../testdata/seurat.DP.RData")

seurat.DP.v3 = UpdateSeuratObject(seurat.DP)
cano_seurat_MAIT = subset(updated_seurat_MAIT, cano == TRUE)
updated_seurat_MAIT@meta.data$ident = updated_seurat_MAIT@active.ident
noncano_cell_num = table(subset(updated_seurat_MAIT@meta.data, cano == FALSE)$ident)

set.seed(123)
M1_seurat.DP.v3 = subset(seurat.DP.v3, cells = sample(colnames(seurat.DP.v3),noncano_cell_num["M1"]))
cano_seurat_MAIT_DP_M1 = merge(cano_seurat_MAIT,M1_seurat.DP.v3)
cano_seurat_MAIT_DP_M1@meta.data[colnames(M1_seurat.DP.v3),]$cano = FALSE
cano_seurat_MAIT_DP_M1@active.ident[cano_seurat_MAIT_DP_M1@active.ident %in% c(11,14)] = "M1"
M1_deg =  FindMarkers(cano_seurat_MAIT_DP_M1, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = c("M1"), group.by="cano", logfc.threshold = 0)
M1_deg$symbol = ensemblGenes[rownames(M1_deg),"external_gene_name"]
M1_deg$sig = M1_deg$p_val_adj < 0.05 & abs(M1_deg$avg_logFC) > 1
table(M1_deg$sig)
M1_deg$symbol = ""
M1_plt = ggplot(M1_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M1") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=15),
        legend.position = "none") #+ ggsave("M1_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

set.seed(123)
M2_seurat.DP.v3 = subset(seurat.DP.v3, cells = sample(colnames(seurat.DP.v3),noncano_cell_num["M2"]))
cano_seurat_MAIT_DP_M2 = merge(cano_seurat_MAIT,M2_seurat.DP.v3)
cano_seurat_MAIT_DP_M2@meta.data[colnames(M2_seurat.DP.v3),]$cano = FALSE
cano_seurat_MAIT_DP_M2@active.ident[cano_seurat_MAIT_DP_M2@active.ident %in% c(11,14)] = "M2"
M2_deg =  FindMarkers(cano_seurat_MAIT_DP_M2, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = c("M2"), group.by="cano", logfc.threshold = 0)
M2_deg$symbol = ensemblGenes[rownames(M2_deg),"external_gene_name"]
M2_deg$sig = M2_deg$p_val_adj < 0.05 & abs(M2_deg$avg_logFC) > 1
table(M2_deg$sig)
M2_deg$symbol = ""
M2_plt = ggplot(M2_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M2") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=15),
        legend.position = "none")

set.seed(123)
M3_seurat.DP.v3 = subset(seurat.DP.v3, cells = sample(colnames(seurat.DP.v3),noncano_cell_num["M3"]))
cano_seurat_MAIT_DP_M3 = merge(cano_seurat_MAIT,M3_seurat.DP.v3)
cano_seurat_MAIT_DP_M3@meta.data[colnames(M3_seurat.DP.v3),]$cano = FALSE
cano_seurat_MAIT_DP_M3@active.ident[cano_seurat_MAIT_DP_M3@active.ident %in% c(11,14)] = "M3"
M3_deg =  FindMarkers(cano_seurat_MAIT_DP_M3, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = c("M3"), group.by="cano", logfc.threshold = 0)
M3_deg$symbol = ensemblGenes[rownames(M3_deg),"external_gene_name"]
M3_deg$sig = M3_deg$p_val_adj < 0.05 & abs(M3_deg$avg_logFC) > 1
table(M3_deg$sig)
M3_deg$symbol = ""
M3_plt = ggplot(M3_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M3") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=15),
        legend.position = "none")


set.seed(123)
M4_seurat.DP.v3 = subset(seurat.DP.v3, cells = sample(colnames(seurat.DP.v3),noncano_cell_num["M4"]))
cano_seurat_MAIT_DP_M4 = merge(cano_seurat_MAIT,M4_seurat.DP.v3)
cano_seurat_MAIT_DP_M4@meta.data[colnames(M4_seurat.DP.v3),]$cano = FALSE
cano_seurat_MAIT_DP_M4@active.ident[cano_seurat_MAIT_DP_M4@active.ident %in% c(11,14)] = "M4"
M4_deg =  FindMarkers(cano_seurat_MAIT_DP_M4, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = c("M4"), group.by="cano", logfc.threshold = 0)
M4_deg$symbol = ensemblGenes[rownames(M4_deg),"external_gene_name"]
M4_deg$sig = M4_deg$p_val_adj < 0.05 & abs(M4_deg$avg_logFC) > 1
table(M4_deg$sig)
M4_deg$symbol = ""
M4_plt = ggplot(M4_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M4") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=15),
        legend.position = "none")

set.seed(123)
M5_seurat.DP.v3 = subset(seurat.DP.v3, cells = sample(colnames(seurat.DP.v3),noncano_cell_num["M5"]))
cano_seurat_MAIT_DP_M5 = merge(cano_seurat_MAIT,M5_seurat.DP.v3)
cano_seurat_MAIT_DP_M5@meta.data[colnames(M5_seurat.DP.v3),]$cano = FALSE
cano_seurat_MAIT_DP_M5@active.ident[cano_seurat_MAIT_DP_M5@active.ident %in% c(11,14)] = "M5"
M5_deg =  FindMarkers(cano_seurat_MAIT_DP_M5, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = c("M5"), group.by="cano", logfc.threshold = 0)
M5_deg$symbol = ensemblGenes[rownames(M5_deg),"external_gene_name"]
M5_deg$sig = M5_deg$p_val_adj < 0.05 & abs(M5_deg$avg_logFC) > 1
table(M5_deg$sig)
M5_deg$symbol = ""
M5_plt = ggplot(M5_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M5") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=15),
        legend.position = "none")

set.seed(123)
M6_seurat.DP.v3 = subset(seurat.DP.v3, cells = sample(colnames(seurat.DP.v3),noncano_cell_num["M6"]))
cano_seurat_MAIT_DP_M6 = merge(cano_seurat_MAIT,M6_seurat.DP.v3)
cano_seurat_MAIT_DP_M6@meta.data[colnames(M6_seurat.DP.v3),]$cano = FALSE
cano_seurat_MAIT_DP_M6@active.ident[cano_seurat_MAIT_DP_M6@active.ident %in% c(11,14)] = "M6"
M6_deg =  FindMarkers(cano_seurat_MAIT_DP_M6, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = c("M6"), group.by="cano", logfc.threshold = 0)
M6_deg$symbol = ensemblGenes[rownames(M6_deg),"external_gene_name"]
M6_deg$sig = M6_deg$p_val_adj < 0.05 & abs(M6_deg$avg_logFC) > 1
table(M6_deg$sig)
M6_deg$symbol = ""
M6_plt = ggplot(M6_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M6") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=15),
        legend.position = "none")

set.seed(123)
M7_seurat.DP.v3 = subset(seurat.DP.v3, cells = sample(colnames(seurat.DP.v3),noncano_cell_num["M7"]))
cano_seurat_MAIT_DP_M7 = merge(cano_seurat_MAIT,M7_seurat.DP.v3)
cano_seurat_MAIT_DP_M7@meta.data[colnames(M7_seurat.DP.v3),]$cano = FALSE
cano_seurat_MAIT_DP_M7@active.ident[cano_seurat_MAIT_DP_M7@active.ident %in% c(11,14)] = "M7"
M7_deg =  FindMarkers(cano_seurat_MAIT_DP_M7, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = c("M7"), group.by="cano", logfc.threshold = 0)
M7_deg$symbol = ensemblGenes[rownames(M7_deg),"external_gene_name"]
M7_deg$sig = M7_deg$p_val_adj < 0.05 & abs(M7_deg$avg_logFC) > 1
table(M7_deg$sig)
M7_deg$symbol = ""
M7_plt = ggplot(M7_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M7") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=15),
        legend.position = "none")


set.seed(123)
M8_seurat.DP.v3 = subset(seurat.DP.v3, cells = sample(colnames(seurat.DP.v3),noncano_cell_num["M8"]))
cano_seurat_MAIT_DP_M8 = merge(cano_seurat_MAIT,M8_seurat.DP.v3)
cano_seurat_MAIT_DP_M8@meta.data[colnames(M8_seurat.DP.v3),]$cano = FALSE
cano_seurat_MAIT_DP_M8@active.ident[cano_seurat_MAIT_DP_M8@active.ident %in% c(11,14)] = "M8"
M8_deg =  FindMarkers(cano_seurat_MAIT_DP_M8, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = c("M8"), group.by="cano", logfc.threshold = 0)
M8_deg$symbol = ensemblGenes[rownames(M8_deg),"external_gene_name"]
M8_deg$sig = M8_deg$p_val_adj < 0.05 & abs(M8_deg$avg_logFC) > 1
table(M8_deg$sig)
M8_deg$symbol = ""
M8_plt = ggplot(M8_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M8") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=15),
        legend.position = "none")


library(grid)
library(gridExtra)



pdf("MAIT_deg_alpha_beta_cano_noncano_volcano_DP_comparison_alpha.pdf",width=16,height =8)
grid.arrange(M1_plt,M2_plt,M3_plt,M4_plt,M5_plt,M6_plt,M7_plt,M8_plt, ncol=4)
dev.off()



shuffle_seurat_MAIT = updated_seurat_MAIT


M1_canos = table(shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M1"],]$cano)
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M1"],]$cano = TRUE
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M1"],]$cano[sample(sum(M1_canos), M1_canos[1])] = FALSE

M1_deg =  FindMarkers(shuffle_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M1", group.by="cano", logfc.threshold = 0)
M1_deg$symbol = ensemblGenes[rownames(M1_deg),"external_gene_name"]
M1_deg$sig = M1_deg$p_val_adj < 0.05 & abs(M1_deg$avg_logFC) > 1
M1_deg$symbol[!M1_deg$sig] = ""
M1_plt = ggplot(M1_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M1") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+ ggsave("M1_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)


M2_canos = table(shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M2"],]$cano)
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M2"],]$cano = TRUE
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M2"],]$cano[sample(sum(M2_canos), M2_canos[1])] = FALSE


M2_deg =  FindMarkers(shuffle_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M2", group.by="cano", logfc.threshold = 0)
M2_deg$symbol = ensemblGenes[rownames(M2_deg),"external_gene_name"]
M2_deg$sig = M2_deg$p_val_adj < 0.05 & abs(M2_deg$avg_logFC) > 1
M2_deg$symbol[!M2_deg$sig] = ""
M2_plt = ggplot(M2_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M2") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M2_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M3_canos = table(shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M3"],]$cano)
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M3"],]$cano = TRUE
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M3"],]$cano[sample(sum(M3_canos), M3_canos[1])] = FALSE


M3_deg =  FindMarkers(shuffle_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M3", group.by="cano", logfc.threshold = 0)
M3_deg$symbol = ensemblGenes[rownames(M3_deg),"external_gene_name"]
M3_deg$sig = M3_deg$p_val_adj < 0.05 & abs(M3_deg$avg_logFC) > 1
M3_deg$symbol[!M3_deg$sig] = ""
M3_plt = ggplot(M3_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant")  +
  ggtitle("M3") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M3_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M4_canos = table(shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M4"],]$cano)
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M4"],]$cano = TRUE
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M4"],]$cano[sample(sum(M4_canos), M4_canos[1])] = FALSE


M4_deg =  FindMarkers(shuffle_seurat_MAIT, ident.1="TRUE", ident.2="FALSE", 
                      subset.ident = "M4", group.by="cano", logfc.threshold = 0)
M4_deg$symbol = ensemblGenes[rownames(M4_deg),"external_gene_name"]
M4_deg$sig = M4_deg$p_val_adj < 0.05 & abs(M4_deg$avg_logFC) > 1
M4_deg$symbol[!(M4_deg$sig)] = ""
M4_plt = ggplot(M4_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant")  +
  ggtitle("M4") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M4_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M5_canos = table(shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M5"],]$cano)
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M5"],]$cano = TRUE
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M5"],]$cano[sample(sum(M5_canos), M5_canos[1])] = FALSE


M5_deg =  FindMarkers(shuffle_seurat_MAIT, ident.1="TRUE", ident.2="FALSE", 
                      subset.ident = "M5", group.by="cano", logfc.threshold = 0)
M5_deg$symbol = ensemblGenes[rownames(M5_deg),"external_gene_name"]
M5_deg$sig = M5_deg$p_val_adj < 0.05 & abs(M5_deg$avg_logFC) > 1
M5_deg$symbol[!M5_deg$sig] = ""
M5_plt = ggplot(M5_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant")  +
  ggtitle("M5") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M5_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M6_canos = table(shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M6"],]$cano)
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M6"],]$cano = TRUE
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M6"],]$cano[sample(sum(M6_canos), M6_canos[1])] = FALSE


M6_deg =  FindMarkers(shuffle_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M6", group.by="cano", logfc.threshold = 0)
M6_deg$symbol = ensemblGenes[rownames(M6_deg),"external_gene_name"]
M6_deg$sig = M6_deg$p_val_adj < 0.05 & abs(M6_deg$avg_logFC) > 1
M6_deg$symbol[!M6_deg$sig] = ""
M6_plt = ggplot(M6_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M6") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M6_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M7_canos = table(shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M7"],]$cano)
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M7"],]$cano = TRUE
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M7"],]$cano[sample(sum(M7_canos), M7_canos[1])] = FALSE


M7_deg =  FindMarkers(shuffle_seurat_MAIT, ident.1="TRUE", ident.2="FALSE",
                      subset.ident = "M7", group.by="cano", logfc.threshold = 0)
M7_deg$symbol = ensemblGenes[rownames(M7_deg),"external_gene_name"]
M7_deg$sig = M7_deg$p_val_adj < 0.05 & abs(M7_deg$avg_logFC) > 1
M7_deg$symbol[!M7_deg$sig] = ""
M7_plt = ggplot(M7_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M7") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none")# +  ggsave("M7_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)

M8_canos = table(shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M8"],]$cano)
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M8"],]$cano = TRUE
shuffle_seurat_MAIT@meta.data[names(shuffle_seurat_MAIT@active.ident)[shuffle_seurat_MAIT@active.ident == "M8"],]$cano[sample(sum(M8_canos), M8_canos[1])] = FALSE


M8_deg =  FindMarkers(shuffle_seurat_MAIT, ident.1="TRUE", ident.2="FALSE", 
                      subset.ident = "M8", group.by="cano", logfc.threshold = 0)
M8_deg$symbol = ensemblGenes[rownames(M8_deg),"external_gene_name"]
M8_deg$sig = M8_deg$p_val_adj < 0.05 & abs(M8_deg$avg_logFC) > 1
M8_deg$symbol[!M8_deg$sig] = ""
M8_plt = ggplot(M8_deg, aes(x=avg_logFC, y= -log10(p_val_adj), colour=sig,label= symbol)) + 
  geom_point() + 
  geom_text() +
  labs(x= "Average Log Fold Change", y = "-Log10(adj. Pvalue)", colour="Significant") +
  ggtitle("M8") +
  ylim(0,45) +
  xlim(-3,3) +
  scale_colour_manual(values=c("grey", "red")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_bw() +
  theme(line = element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.position = "none") #+  ggsave("M8_canonical_noncanonical_DEG_volcano.pdf", width = 7, height=7)


pdf("MAIT_deg_alpha_beta_cano_noncano_volcano_shuffled_alpha.pdf",width=16,height =8)
grid.arrange(M1_plt,M2_plt,M3_plt,M4_plt,M5_plt,M6_plt,M7_plt,M8_plt, ncol=4)
dev.off()
