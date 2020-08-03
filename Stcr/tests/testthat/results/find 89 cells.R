# init
library(Stcr)
library(Seurat)
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
chain2_v_gene = "(TRBV13|TRBV19)"
# chain2_v_gene = NULL # 2nd revision
chain2_j_gene = NULL
df_clonotype2 = subset(df_clonotype, n > 0)
df_clonotype2$selected = FALSE

grepl_chain1_v_gene = sum(df_clonotype2[grepl(chain1_v_gene,df_clonotype2$chain1_v_gene),]$n)
grepl_chain1_j_gene = grepl(chain1_j_gene,df_clonotype2$chain1_j_gene)
grepl_chain2_v_gene = grepl(chain2_v_gene,df_clonotype2$chain2_v_gene)
not_grepl_chain1_v_gene = !grepl(chain1_v_gene,df_clonotype2$chain1_v_gene)
not_grepl_chain1_j_gene = !grepl(chain1_j_gene,df_clonotype2$chain1_j_gene)
not_grepl_chain2_v_gene = !grepl(chain2_v_gene,df_clonotype2$chain2_v_gene)


sum(df_clonotype2[!grepl(chain1_v_gene,df_clonotype2$chain1_v_gene),]$n)
sum(df_clonotype2[!grepl(chain1_j_gene,df_clonotype2$chain1_j_gene),]$n)
sum(df_clonotype2[!grepl(chain2_v_gene,df_clonotype2$chain2_v_gene),]$n)


sum(df_clonotype2[!grepl(chain1_v_gene,df_clonotype2$chain1_v_gene) & !grepl(chain1_j_gene,df_clonotype2$chain1_j_gene),]$n)
sum(df_clonotype2[!grepl(chain2_v_gene,df_clonotype2$chain2_v_gene) & !grepl(chain1_j_gene,df_clonotype2$chain1_j_gene),]$n)
sum(df_clonotype2[!grepl(chain2_v_gene,df_clonotype2$chain2_v_gene) & !grepl(chain1_v_gene,df_clonotype2$chain1_v_gene),]$n)

sum(df_clonotype2[!grepl(chain2_v_gene,df_clonotype2$chain2_v_gene) & !grepl(chain1_v_gene,df_clonotype2$chain1_v_gene) & !grepl(chain1_j_gene,df_clonotype2$chain1_j_gene),]$n)

df_clonotype2[!grepl(chain2_v_gene,df_clonotype2$chain2_v_gene) & !grepl(chain1_v_gene,df_clonotype2$chain1_v_gene) & !grepl(chain1_j_gene,df_clonotype2$chain1_j_gene),]$selected = TRUE


df_clonotype_labeled = df_clonotype2
nonselected_clonotypes = as.vector(df_clonotype_labeled$clonotype)[df_clonotype_labeled$selected == FALSE]
selected_clonotypes = as.vector(df_clonotype_labeled$clonotype)[df_clonotype_labeled$selected == TRUE]

nonselected_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,nonselected_clonotypes]) > 0]
selected_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,selected_clonotypes]) > 0]
seurat.MAIT@scale.data = as.matrix(seurat.MAIT@scale.data)
seurat.MAIT_v3 = UpdateSeuratObject(seurat.MAIT)
seurat.MAIT.removed89 = subset(seurat.MAIT_v3, cells =  colnames(seurat.MAIT@scale.data)[!(colnames(seurat.MAIT@scale.data) %in% selected_ind)])



sset_MAIT2 = ImportSeurat(sset,seurat.MAIT.removed89, version="v3")

MAIT_cdr3_mat2 = SetCDR3Mat(sset_MAIT2, cell_ind = colnames(seurat.MAIT.removed89),
                           chain1 = c("TRA"), chain2 = c("TRB"))

MAIT_cdr3_pair_mat2 = SetCDR3PairMat(sset_MAIT2, MAIT_cdr3_mat2, chain1 = c("TRA"), chain2 = c("TRB"))


mat_cdr3_pair2 = MAIT_cdr3_pair_mat2[[1]]
df2_clonotype = MAIT_cdr3_pair_mat2[[2]]
ident = as.factor(seurat.MAIT.removed89@active.ident)
df_mat_cdr3_pair2 = data.frame(mat_cdr3_pair2)
df_mat_cdr3_pair2$ident = ""
for(i in levels(ident)){
  cell_ind = names(ident[ident %in% c(i)])
  df_mat_cdr3_pair2$ident[rownames(df_mat_cdr3_pair2) %in% cell_ind] = i
}

for(i in 1:nrow(df2_clonotype)){
  df2_clonotype$ident[i] =  paste( unique(df_mat_cdr3_pair2[mat_cdr3_pair2[,i] > 0,'ident']), collapse = ",")
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




plt4 = PlotFrequency(MAIT_cdr3_pair_mat2, ident = seurat.MAIT.removed89@active.ident, reverse = TRUE, selected_ind = NULL,
                     chain1_v_gene = "TRAV1\\*", chain1_d_gene = NULL, chain1_j_gene = "TRAJ33", 
                     chain2_v_gene = NULL, chain2_d_gene = NULL, chain2_j_gene = NULL)

plt4 = plt4 + theme_bw() + ylab("") +
  geom_point(size=7) +
  geom_line(group = 0,size=2) +
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
    legend.position = "none") + ggsave("MAIT_noncano_freq_alphachain_removed89cells.pdf",width=7,height=7)
