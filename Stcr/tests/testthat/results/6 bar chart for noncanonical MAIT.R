# six bar chart
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

df_clonotype2_noncano = subset(df_clonotype_labeled, selected == FALSE)


grepl_chain1_v_gene = grepl(chain1_v_gene,df_clonotype2_noncano$chain1_v_gene)
grepl_chain1_j_gene = grepl(chain1_j_gene,df_clonotype2_noncano$chain1_j_gene)
grepl_chain2_v_gene = grepl(chain2_v_gene,df_clonotype2_noncano$chain2_v_gene)
TRAV_noncano_clonotype = as.vector(df_clonotype2_noncano[!grepl_chain1_v_gene,]$clonotype)
TRAV_noncano_idx = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,TRAV_noncano_clonotype]) > 0]
TRAJ_noncano_clonotype = as.vector(df_clonotype2_noncano[!grepl_chain1_j_gene,]$clonotype)
TRAJ_noncano_idx = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,TRAJ_noncano_clonotype]) > 0]
TRBV_noncano_clonotype =  as.vector(df_clonotype2_noncano[!grepl_chain2_v_gene,]$clonotype)
TRBV_noncano_idx = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,TRBV_noncano_clonotype]) > 0]

sum(df_clonotype2$n[df_clonotype2$clonotype %in% TRBV_noncano_clonotype])


# cell89_clonotype = intersect(intersect(TRAV_noncano_clonotype,TRAJ_noncano_clonotype),TRBV_noncano_clonotype)
# cell89_index = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,cell89_clonotype]) > 0]
# length(cell89_index)
cell89_index = intersect(intersect(TRAV_noncano_idx,TRAJ_noncano_idx),TRBV_noncano_idx)
int_TRAV_TRAJ = intersect(TRAV_noncano_idx,TRAJ_noncano_idx)
cell103_index = int_TRAV_TRAJ[!(int_TRAV_TRAJ %in% cell89_index)]
cell3_index = TRAV_noncano_idx[!(TRAV_noncano_idx %in% TRAJ_noncano_idx) & !(TRAV_noncano_idx %in% TRBV_noncano_idx)]
cell1_index = TRAV_noncano_idx[!(TRAV_noncano_idx %in% TRAJ_noncano_idx) & !(TRAV_noncano_idx %in% cell3_index)]
cell21_index = TRAJ_noncano_idx[!(TRAJ_noncano_idx %in% TRAV_noncano_idx) & !(TRAJ_noncano_idx %in% TRBV_noncano_idx)]
cell4_index = TRAJ_noncano_idx[!(TRAJ_noncano_idx %in% cell21_index) & !(TRAJ_noncano_idx %in% TRAV_noncano_idx)]

cell448_index = TRBV_noncano_idx[!(TRBV_noncano_idx %in% TRAV_noncano_idx) & !(TRBV_noncano_idx %in% TRAJ_noncano_idx)]
cell448_index = TRBV_noncano_idx[!(TRBV_noncano_idx %in% TRAV_noncano_idx) & !(TRBV_noncano_idx %in% TRAJ_noncano_idx)]
newclusters = as.vector(seurat.MAIT@ident)
names(newclusters) = names(seurat.MAIT@ident)
newclusters = factor(newclusters)
ggdf = data.frame(cluster = seurat.MAIT@ident,
                  cell3 = FALSE,
                  cell1 = FALSE,
                  cell103 = FALSE,
                  cell89 = FALSE,
                  cell21 = FALSE,
                  cell4 = FALSE,
                  cell448= FALSE)
ggdf[names(seurat.MAIT@ident) %in% cell3_index,]$cell3 = TRUE
ggdf[names(seurat.MAIT@ident) %in% cell1_index,]$cell1 = TRUE
ggdf[names(seurat.MAIT@ident) %in% cell103_index,]$cell103 = TRUE
ggdf[names(seurat.MAIT@ident) %in% cell89_index,]$cell89 = TRUE
ggdf[names(seurat.MAIT@ident) %in% cell21_index,]$cell21 = TRUE
ggdf[names(seurat.MAIT@ident) %in% cell4_index,]$cell4 = TRUE
ggdf[names(seurat.MAIT@ident) %in% cell448_index,]$cell448 = TRUE
table(ggdf$cell3)
library(ggplot2)
library(dplyr)
subggdf = subset(ggdf, cell448 == TRUE)
subggdf %>% group_by(cluster) %>% summarise(count = n())


cluster = c("M1","M2","M3","M4","M5","M6", "M7", "M8")
cell3 = c(1,0,2,0,0,0,0,0)
cell1 = c(0,0,0,1,0,0,0,0)
cell103=c(40,5,25,8,2,4,8,11)
cell89 =c(42,1,30,6,0,0,7,3)
cell21 =c(5,0,6,1,4,0,1,4)
cell4 = c(1,0,2,0,0,0,0,1)
cell448 = c(130,8,147,45,8,12,49,49)

gggdf = data.frame(cluster=cluster,
                   cell3=cell3,
                   cell1=cell1,
                   cell103=cell103,
                   cell89=cell89,
                   cell21=cell21,
                   cell4=cell4,
                   cell448=cell448)

plt1 = ggplot(gggdf) + geom_bar(aes(x=cluster, y= cell3),stat="identity") +
  ylim(0,5) +
  theme_bw() +
  theme(line = element_blank(),
        plot.margin = margin(0.2,0.2,0.2,1.2, "cm"),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=30),
        axis.text.x = element_blank(),
        legend.position = "none")
plt2 = ggplot(gggdf) + geom_bar(aes(x=cluster, y= cell1),stat="identity") +
  ylim(0,5) +
  theme_bw() +
  theme(line = element_blank(),
        plot.margin = margin(0.2,0.2,0.2,1.2, "cm"),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=30),
        axis.text.x = element_blank(),
        legend.position = "none")
plt3 = ggplot(gggdf) + geom_bar(aes(x=cluster, y= cell103),stat="identity") +
  ylim(0,45) +
  theme_bw() +
  theme(line = element_blank(),
        plot.margin = margin(0.2,0.2,0.2,.7, "cm"),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=30),
        axis.text.x = element_blank(),
        legend.position = "none")
plt4 = ggplot(gggdf) + geom_bar(aes(x=cluster, y= cell89),stat="identity") +
  ylim(0,45) +
  theme_bw() +
  theme(line = element_blank(),
        plot.margin = margin(0.2,0.2,0.2,.7, "cm"),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=30),
        axis.text.x = element_blank(),
        legend.position = "none")
plt5 = ggplot(gggdf) + geom_bar(aes(x=cluster, y= cell21),stat="identity") +
  ylim(0,10) +
  scale_y_continuous(breaks=c(0,2,4,6,8)) +
  theme_bw() +
  theme(line = element_blank(),
        plot.margin = margin(0.2,0.2,0.2,1.2, "cm"),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=30),
        axis.text.x = element_blank(),
        legend.position = "none")
plt6 = ggplot(gggdf) + geom_bar(aes(x=cluster, y= cell4),stat="identity") +
  ylim(0,5) +
  theme_bw() +
  theme(line = element_blank(),
        plot.margin = margin(0.2,0.2,0.2,1.2, "cm"),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=30),
        axis.text.x = element_blank(),
        legend.position = "none")
plt7 = ggplot(gggdf) + geom_bar(aes(x=cluster, y= cell448),stat="identity") +
  ylim(0,150) +
  theme_bw() +
  theme(line = element_blank(),
        plot.margin = margin(0.2,0.2,0.2,0, "cm"),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=30),
        axis.text.x = element_blank(),
        legend.position = "none")


library(grid)
library(gridExtra)

pdf("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat/bar_graph_noncanonical_MAIT.pdf", width=7, height =15)
grid.arrange(plt1,plt2,plt3,plt4,plt5,plt6,plt7, ncol=1)
dev.off()
