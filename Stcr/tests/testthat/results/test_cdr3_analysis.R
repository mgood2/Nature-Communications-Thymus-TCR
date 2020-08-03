# working dir
setwd("D:/espark/iNKT/data/Seuratset/")
NKT_subset <- readRDS("D:/espark/iNKT/data/Seuratset/af_remove_NKTSUbset.rds")
DimPlot(NKT_subset)
# library
library(Seurat)
library(grid)
library(gridExtra)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(pheatmap)
#index
NKT_ind = colnames(NKT_subset@assays$RNA@data)
#load contig file
G5_vdj = read.csv("D:/espark/iNKT/mNKT-NCD_S2/all_contig_annotations.csv", row.names = "contig_id")
G6_vdj = read.csv("D:/espark/iNKT/mNKT-HFD_S2_IMGT/all_contig_annotations.csv", row.names = "contig_id")
G7_vdj = read.csv("D:/espark/iNKT/mNKT-HFD-aGC_S2_IMGT/all_contig_annotations.csv", row.names = "contig_id")
#Quality control TCR
G5_vdj_prod_highconf = subset(G5_vdj, productive == "True" & high_confidence=="True") 
G6_vdj_prod_highconf = subset(G6_vdj, productive == "True" & high_confidence=="True")
G7_vdj_prod_highconf = subset(G7_vdj, productive == "True" & high_confidence=="True")
#TAGGing
G5_vdj[,"barcode"] = sprintf('%s_G5', G5_vdj[,"barcode"])
G6_vdj[,"barcode"] = sprintf('%s_G6', G6_vdj[,"barcode"])
G7_vdj[,"barcode"] = sprintf('%s_G7', G7_vdj[,"barcode"])
G5_vdj_prod_highconf$cellname = sprintf('%s_G5', G5_vdj_prod_highconf$barcode)
G6_vdj_prod_highconf$cellname = sprintf('%s_G6', G6_vdj_prod_highconf$barcode)
G7_vdj_prod_highconf$cellname = sprintf('%s_G7', G7_vdj_prod_highconf$barcode)



all_vdj_prod_highconf <- rbind(G5_vdj_prod_highconf,G6_vdj_prod_highconf,G7_vdj_prod_highconf)
all_vdj_prod_highconf_subset <- subset(all_vdj_prod_highconf, all_vdj_prod_highconf$cellname %in% NKT_ind)
all_vdj_prod_highconf_subset_TRA <- subset(all_vdj_prod_highconf_subset, all_vdj_prod_highconf_subset$chain =="TRA")
A_cdr <- unique(as.vector(all_vdj_prod_highconf_subset_TRA$cdr3))
df = rbind(G5_vdj_prod_highconf[,c("cdr3", "cellname")],G6_vdj_prod_highconf[,c("cdr3", "cellname")],G7_vdj_prod_highconf[,c("cdr3", "cellname")])
# only TCR in NKT cells
df = subset(df, cellname %in% NKT_ind)
row_cells = unique(df$cellname)
A_cdr3_list = A_cdr
mat_Acdr = matrix(data=0, ncol=length(A_cdr3_list), nrow=length(row_cells))
rownames(mat_Acdr) = row_cells
colnames(mat_Acdr) = A_cdr3_list
for(i in rownames(mat_Acdr)){
  for(j in colnames(mat_Acdr)){
    if(j %in% as.vector(subset(df, cellname %in% i)$cdr3)){
      mat_Acdr[i,j] = mat_Acdr[i,j] +1
    }
  }
}

for (i in unique(all_vdj_prod_highconf_subset$cellname[which(all_vdj_prod_highconf_subset$cellname %in% rownames(mat_Acdr)[rowSums(mat_Acdr)>1])])){
  mat_Acdr[i,] =0
  max_gene <- all_vdj_prod_highconf_subset$cdr3[which(max(all_vdj_prod_highconf_subset$reads[which(all_vdj_prod_highconf_subset$cellname == i & all_vdj_prod_highconf_subset$chain =="TRA")])==all_vdj_prod_highconf_subset$reads & all_vdj_prod_highconf_subset$cellname ==i)]
  mat_Acdr[i ,as.character(max_gene)] =1
}
table(rowSums(mat_Acdr))



all_vdj_prod_highconf_subset_TRB <- subset(all_vdj_prod_highconf_subset, all_vdj_prod_highconf_subset$chain =="TRB")
B_cdr <- unique(as.vector(all_vdj_prod_highconf_subset_TRB$cdr3))
df = rbind(G5_vdj_prod_highconf[,c("cdr3", "cellname")],G6_vdj_prod_highconf[,c("cdr3", "cellname")],G7_vdj_prod_highconf[,c("cdr3", "cellname")])
# only TCR in NKT cells
df = subset(df, cellname %in% NKT_ind)
row_cells = unique(df$cellname)
B_cdr3_list = B_cdr
mat_Bcdr = matrix(data=0, ncol=length(B_cdr3_list), nrow=length(row_cells))
rownames(mat_Bcdr) = row_cells
colnames(mat_Bcdr) = B_cdr3_list
for(i in rownames(mat_Bcdr)){
  for(j in colnames(mat_Bcdr)){
    if(j %in% as.vector(subset(df, cellname %in% i)$cdr3)){
      mat_Bcdr[i,j] = mat_Bcdr[i,j] +1
    }
  }
}

for (i in unique(all_vdj_prod_highconf_subset$cellname[which(all_vdj_prod_highconf_subset$cellname %in% rownames(mat_Bcdr)[rowSums(mat_Bcdr)>1])])){
  mat_Bcdr[i,] =0
  max_gene <- all_vdj_prod_highconf_subset$cdr3[which(max(all_vdj_prod_highconf_subset$reads[which(all_vdj_prod_highconf_subset$cellname == i & all_vdj_prod_highconf_subset$chain =="TRB")])==all_vdj_prod_highconf_subset$reads & all_vdj_prod_highconf_subset$cellname ==i)]
  mat_Bcdr[i ,as.character(max_gene)] =1
}



mat_ABcdr = ((t(mat_Acdr) %*% mat_Bcdr)/length(row_cells))

setwd("D:/espark/iNKT/clonotype/real/cdr3")
pheatmap(mat_ABcdr,cluster_cols = F,cluster_rows = F, color =colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         filename="Acdr vs Bcdr cellfreq.pdf")

table(mat_ABcdr[,1])
table(mat_ABcdr)
mat_Log_ABcdr = log10(mat_ABcdr+1)
table(mat_Log_ABcdr)

##########yogi
setwd("D:/espark/iNKT/clonotype/real/cdr3")
for(i in 0:(nlevels(NKT_subset@active.ident)-1)){
  cluster_cell_ind = names(NKT_subset@active.ident)[NKT_subset@active.ident == i]
  num_cluster_cell = length(cluster_cell_ind)
  c_mat_Acdr = mat_Acdr[rownames(mat_Acdr) %in% cluster_cell_ind,]
  c_mat_Bcdr = mat_Bcdr[rownames(mat_Bcdr) %in% cluster_cell_ind,]
  c_mat_ABcdr = (t(c_mat_Acdr) %*% c_mat_Bcdr)/num_cluster_cell
  log2FC_c_mat_ABcdr = log2((c_mat_ABcdr+1)/(mat_ABcdr+1))
  paletteLength <- 100
  if(max(log2FC_c_mat_ABcdr) > 0 ){
    myBreaks <- c(seq(min(log2FC_c_mat_ABcdr), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(log2FC_c_mat_ABcdr)/paletteLength, max(log2FC_c_mat_ABcdr), length.out=floor(paletteLength/2)))
    pheatmap(log2FC_c_mat_ABcdr,breaks = myBreaks,cluster_cols = F,cluster_rows = F, color =colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             filename=paste0("Acdr vs Bcdr log2FC_G",i,".pdf"))
  } else {
    myBreaks <- c(seq(min(log2FC_c_mat_ABcdr), 0, length.out=ceiling(paletteLength/2) + 1))
    pheatmap(log2FC_c_mat_ABcdr,breaks = myBreaks,cluster_cols = F,cluster_rows = F, color =colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)[1:50],
             
             filename=paste0("Acdr vs Bcdr log2FC_G",i,".pdf"))
  }
}

shannon_index <- function(vector){
  k = length(vector)
  return(-sum((vector/sum(vector))*log(vector/sum(vector)), na.rm=TRUE)/log(k))
}
plot_DI_cluster <- function(df, name){
  ggplot(df, aes(x=cluster,
                 y=cluster_shannon_indices,
                 group=1)) + geom_line() + geom_point(size=3)  +
    ggtitle(name) +
    xlab("") +
    ylab("Normalized diversity Index") +
    ylim(low=0, high=1) +
    # annotate(x=0, xend=0, y=0, yend=1, lwd=0.75, geom="segment") +
    theme(
      # axis.text.x=element_text(size=20),
      # axis.text.y=element_text(size=20),
      # axis.title.x=element_blank(),
      # axis.title.y=element_text(size=20),
      # plot.background=element_blank(),
      # panel.grid.major=element_blank(),
      # panel.grid.minor=element_blank(), 
      # axis.line.y = element_blank(),
      # axis.ticks.x =element_blank()
    ) + ggsave(paste0("Diversity_index_cluster_",name,"_NKT.pdf"))
  dev.off()
}
plot_DI_TCR <- function(df, name){
  ggplot(df, aes(x=cluster,
                 y=tcr_shannon_indices,
                 group=1)) + geom_line() + geom_point(size=3)  +
    ggtitle(name) +
    xlab("") +
    ylab("Normalized diversity Index") +
    ylim(low=0, high=1) +
    # annotate(x=0, xend=0, y=0, yend=1, lwd=0.75, geom="segment") +
    theme(
      axis.text.x=element_text(angle=60,hjust = 1),
      # axis.text.y=element_text(angle=90),
      # axis.title.x=element_blank(),
      # axis.title.y=element_text(size=20),
      # plot.background=element_blank(),
      # panel.grid.major=element_blank(),
      # panel.grid.minor=element_blank(), 
      # axis.line.y = element_blank(),
      # axis.ticks.x =element_blank()
    ) + ggsave(paste0("Diversity_index_TCR_",name,"_NKT.pdf"))
  dev.off()
}
plot_bar_TCR <- function(df,name){
  ggplot(df ,aes(y = value, x = Var2, fill=Var1)) + 
    geom_bar(stat="identity",position="fill") +
    # scale_x_discrete(limits=seq(0,(nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_clean@ident)-1))) +
    ggtitle(name) +
    xlab("") +
    ylab("TCR type frequency (%)") +
    scale_y_continuous(labels = scales::percent) +
    annotate(x=0, xend=0, y=0, yend=1, lwd=0.75, geom="segment") +
    theme(
      axis.text.x=element_text(size=20),
      axis.text.y=element_text(size=20),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=20),
      plot.background=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(), 
      axis.line = element_blank(),
      axis.ticks.x =element_blank()
    ) + ggsave(paste0(name,"_bar_graph_by_cluster_NKT.pdf"), width =9)
  dev.off()
}


# DV diversity index
# number of cells > 10  that has TCR

Acdr_used = names(colSums(mat_Acdr))
mat_cluster_TCR = matrix(data=0, ncol=nlevels(NKT_subset@active.ident),
                         nrow = length(Acdr_used))
rownames(mat_cluster_TCR) = Acdr_used
colnames(mat_cluster_TCR) = paste0("A",levels(NKT_subset@active.ident))
for(i in 0:(nlevels(NKT_subset@active.ident)-1)){
  cluster_cell_ind = names(NKT_subset@active.ident)[NKT_subset@active.ident == i]
  c_mat_Acdr = mat_Acdr[rownames(mat_Acdr) %in% cluster_cell_ind,]
  c_mat_Acdr = c_mat_Acdr[,Acdr_used]
  mat_cluster_TCR[,(i+1)] = colSums(c_mat_Acdr)
}

cluster_shannon_indices = apply(mat_cluster_TCR, 2,shannon_index)
cluster_shannon_indices1 = cluster_shannon_indices
df = data.frame(cluster_shannon_indices,
                cluster = names(cluster_shannon_indices))

plot_DI_cluster(df,"Acdr")

dev.off()

tcr_shannon_indices = apply(mat_cluster_TCR, 1,shannon_index)
df = data.frame(tcr_shannon_indices,
                cluster = names(tcr_shannon_indices))
plot_DI_TCR(df,"Acdr")

mat_cluster_TCR_long = melt(mat_cluster_TCR)
plot_bar_TCR(mat_cluster_TCR_long,"Acdr")

#Beta chain CDR

Bcdr_used = names(colSums(mat_Bcdr))
mat_cluster_TCR = matrix(data=0, ncol=nlevels(NKT_subset@active.ident),
                         nrow = length(Bcdr_used))
rownames(mat_cluster_TCR) = Bcdr_used
colnames(mat_cluster_TCR) = paste0("B",levels(NKT_subset@active.ident))
for(i in 0:(nlevels(NKT_subset@active.ident)-1)){
  cluster_cell_ind = names(NKT_subset@active.ident)[NKT_subset@active.ident == i]
  c_mat_Bcdr = mat_Bcdr[rownames(mat_Bcdr) %in% cluster_cell_ind,]
  c_mat_Bcdr = c_mat_Bcdr[,Bcdr_used]
  mat_cluster_TCR[,(i+1)] = colSums(c_mat_Bcdr)
}

cluster_shannon_indices = apply(mat_cluster_TCR, 2,shannon_index)
cluster_shannon_indices1 = cluster_shannon_indices
df = data.frame(cluster_shannon_indices,
                cluster = names(cluster_shannon_indices))

plot_DI_cluster(df,"Bcdr")

dev.off()

tcr_shannon_indices = apply(mat_cluster_TCR, 1,shannon_index)
df = data.frame(tcr_shannon_indices,
                cluster = names(tcr_shannon_indices))
plot_DI_TCR(df,"Bcdr")

mat_cluster_TCR_long = melt(mat_cluster_TCR)
plot_bar_TCR(mat_cluster_TCR_long,"Bcdr")




#paired cluster 
paired_mat_ABcdr = t(mat_Acdr) %*% mat_Bcdr
#using cells > 10 TCR
table(paired_mat_ABcdr)
mat_ABcdr_ = matrix(data = 0, nrow = nrow(mat_Acdr), ncol = table(paired_mat_ABcdr > 0)[2])
rownames(mat_ABcdr_) = rownames(mat_Acdr)
cnt = 1
pair_names = c()
for(i in rownames(paired_mat_ABcdr)){
  # print(names((paired_mat_GVBV > 50)[i,])[(paired_mat_GVBV > 50)[i,] == TRUE])
  if(length(names((paired_mat_ABcdr > 0)[i,])[(paired_mat_ABcdr > 0)[i,] == TRUE])> 0){
    print(i)
    for(j in names((paired_mat_ABcdr > 0)[i,])[(paired_mat_ABcdr > 0)[i,] == TRUE]){
      print(j)
      pair_ABcdr = mat_Acdr[,i] * mat_Bcdr[,j]
      mat_ABcdr_[,cnt] = pair_ABcdr
      pair_names = c(pair_names,paste0(i,"-",j))
      cnt = cnt + 1
    }
    print("---------")
  }
}
colnames(mat_ABcdr_) =  pair_names

# GVBV diversity index
ABcdr_used = names(colSums(mat_ABcdr_))[colSums(mat_ABcdr_) > 0]
mat_cluster_TCR = matrix(data=0, ncol=nlevels(NKT_subset@active.ident),
                         nrow = length(ABcdr_used))
rownames(mat_cluster_TCR) = ABcdr_used
colnames(mat_cluster_TCR) = paste0("A",levels(NKT_subset@active.ident))
for(i in 0:(nlevels(NKT_subset@active.ident)-1)){
  cluster_cell_ind = names(NKT_subset@active.ident)[NKT_subset@active.ident == i]
  c_mat_ABcdr = mat_ABcdr_[rownames(mat_ABcdr_) %in% cluster_cell_ind,]
  c_mat_ABcdr = c_mat_ABcdr[,ABcdr_used]
  mat_cluster_TCR[,(i+1)] = colSums(c_mat_ABcdr)
}
cluster_shannon_indices = apply(mat_cluster_TCR, 2,shannon_index)
cluster_shannon_indices3 = cluster_shannon_indices
df = data.frame(cluster_shannon_indices,
                cluster = names(cluster_shannon_indices)
)
plot_DI_cluster(df,"ABcdr")



# other pair

ABcdr_used = names(colSums(mat_ABcdr_))[colSums(mat_ABcdr_) > 1]
mat_cluster_TCR = matrix(data=0, ncol=nlevels(NKT_subset@active.ident),
                         nrow = length(ABcdr_used))
rownames(mat_cluster_TCR) = ABcdr_used
colnames(mat_cluster_TCR) = paste0("A",levels(NKT_subset@active.ident))
for(i in 0:(nlevels(NKT_subset@active.ident)-1)){
  cluster_cell_ind = names(NKT_subset@active.ident)[NKT_subset@active.ident == i]
  c_mat_ABcdr = mat_ABcdr_[rownames(mat_ABcdr_) %in% cluster_cell_ind,]
  c_mat_ABcdr = c_mat_ABcdr[,ABcdr_used]
  mat_cluster_TCR[,(i+1)] = colSums(c_mat_ABcdr)
}

ABcdr_other = names(colSums(mat_ABcdr_))[colSums(mat_ABcdr_) <= 1]
mat_cluster_TCR_other = matrix(data=0, ncol=nlevels(NKT_subset@active.ident),
                               nrow = length(ABcdr_other))
rownames(mat_cluster_TCR_other) = ABcdr_other
colnames(mat_cluster_TCR_other) = paste0("A",levels(NKT_subset@active.ident))
for(i in 0:(nlevels(NKT_subset@active.ident)-1)){
  cluster_cell_ind = names(NKT_subset@active.ident)[NKT_subset@active.ident == i]
  c_mat_ABcdr = mat_ABcdr_[rownames(mat_ABcdr_) %in% cluster_cell_ind,]
  c_mat_ABcdr = c_mat_ABcdr[,ABcdr_other]
  mat_cluster_TCR_other[,(i+1)] = colSums(c_mat_ABcdr)
}

tcr_shannon_indices = apply(mat_cluster_TCR, 1,shannon_index)
df = data.frame(tcr_shannon_indices,
                cluster = factor(names(tcr_shannon_indices),levels=c(ABcdr_used,"other")))
plot_DI_TCR(df,"ABcdr_over_1")



mat_cluster_TCR = rbind(mat_cluster_TCR, colSums(mat_cluster_TCR_other))
rownames(mat_cluster_TCR)[length(rownames(mat_cluster_TCR))] = "other"



# mat_cluster_TCR = mat_cluster_TCR[,1:7]




gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors = gg_color_hue(30)
colors[30] = "grey"

mat_cluster_TCR
mat_cluster_TCR_long = melt(mat_cluster_TCR)
# plot_bar_TCR(mat_cluster_TCR_long,"GVDV")
ggplot(mat_cluster_TCR_long ,aes(y = value, x = Var2, fill=Var1)) + 
  geom_bar(stat="identity",position="fill") +
  scale_fill_manual(values=colors) +
  ggtitle("AVBV") +
  xlab("") +
  ylab("TCR type frequency (%)") +
  scale_y_continuous(labels = scales::percent) +
  annotate(x=0, xend=0, y=0, yend=1, lwd=0.75, geom="segment") +
  theme(
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=20),
    plot.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), 
    axis.line = element_blank(),
    axis.ticks.x =element_blank()
  ) + ggsave(paste0("AVBV","_bar_graph_by_cluster_NKT.pdf"), width =9)
dev.off()

df = data.frame(x1=cluster_shannon_indices1,
                x2=cluster_shannon_indices2,
                x3=cluster_shannon_indices3,
                cluster = names(cluster_shannon_indices1))
ggplot(df) + geom_line(aes(x=cluster,y=x1,group=1,colour="BV")) + geom_point(aes(x=cluster,y=x1,group=1,colour="BV"),size=3)  +
  geom_line(aes(x=cluster,y=x2,group=2,colour="AV")) + geom_point(aes(x=cluster,y=x2,group=2,colour="AV"),size=3)  +
  geom_line(aes(x=cluster,y=x3,group=3,colour="AVBV")) + geom_point(aes(x=cluster,y=x3,group=3,colour="AVBV"),size=3)  +
  scale_colour_manual(values =c("BV"="#66c2a5",
                                "AV"="#fc8d62",
                                "AVBV"="#8da0cb")) +
  xlab("") +
  labs(colour="TCR") +
  ggtitle("AV,BV,AVBV") +
  ylab("Normalized diversity Index") +
  ylim(low=0, high=1) +
  # annotate(x=0, xend=0, y=0, yend=1, lwd=0.75, geom="segment") +
  theme(
    # axis.text.x=element_text(size=20),
    # axis.text.y=element_text(size=20),
    # axis.title.x=element_blank(),
    # axis.title.y=element_text(size=20),
    # plot.background=element_blank(),
    # panel.grid.major=element_blank(),
    # panel.grid.minor=element_blank(), 
    # axis.line.y = element_blank(),
    # axis.ticks.x =element_blank()
  ) + ggsave(paste0("Diversity_index_cluster_AV_BV_AVBV_NKT.pdf"))
dev.off()



# clonotype counting

# use gdT.R to count in the scell cluster

load("mat_alpha_NKT.RData")
load("mat_beta_NKT.RData")
load("mat_NKT_cdr3.RData")

B_cdr3_ind = names(colSums(mat_NKT_cdr3))[colSums(mat_NKT_cdr3) > 0]
A_cdr3_ind = names(rowSums(mat_NKT_cdr3))[rowSums(mat_NKT_cdr3) > 0]

mat_gdT_cdr3_2 =mat_gdT_cdr3[G_cdr3_ind,D_cdr3_ind]

# G - D pairs
cdr3_pairs = c()
for(i in 1:nrow(mat_gdT_cdr3_2)){
  for(j in 1:ncol(mat_gdT_cdr3_2)){
    if(mat_gdT_cdr3_2[i,j] > 0){
      cdr3_pairs = c(cdr3_pairs, paste0(rownames(mat_gdT_cdr3_2)[i],"-",colnames(mat_gdT_cdr3_2)[j]))
    }
  }
}

# length(gdT_ind) # cell 2681
# length(rownames(mat_gamma)) # cell 2637
# mat_cdr3_pair_gdT =  matrix(data = 0 , nrow = length(rownames(mat_gamma)), ncol = length(cdr3_pairs))
# rownames(mat_cdr3_pair_gdT) = rownames(mat_gamma)
# colnames(mat_cdr3_pair_gdT) = cdr3_pairs
# for(i in 1:length(cdr3_pairs)){
#   cdr3_pair = strsplit(cdr3_pairs[i], "-")[[1]]
#   G_TCR = cdr3_pair[1]
#   D_TCR = cdr3_pair[2]
#   for(j in rownames(mat_cdr3_pair_gdT)){
#     TCRS = subset(gd_vdj_prod_highconf, cellname %in% j)$cdr3
#     if(G_TCR %in% TCRS & D_TCR %in% TCRS){
#       mat_cdr3_pair_gdT[j,i] = mat_cdr3_pair_gdT[j,i] + 1
#     }
#   }
# }
# save(mat_cdr3_pair_gdT, file="mat_cdr3_pair_gdT.RData")
load(file="mat_cdr3_pair_gdT.RData")
use_one_clonotype_ind = rownames(mat_cdr3_pair_gdT)[rowSums(mat_cdr3_pair_gdT) == 1]

mat_delta2 = mat_delta[use_one_clonotype_ind,]
mat_gamma2 = mat_gamma[use_one_clonotype_ind,]

mat_cluster_clonotype_counts = matrix(data=0, ncol=nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident),
                                      nrow =3)
colnames(mat_cluster_clonotype_counts) = levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)
rownames(mat_cluster_clonotype_counts) = c(1,2,">=3")

for(i in 0:(nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)-1)){
  cluster_cell_ind = names(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)[seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident == i]
  mat_gdT_cdr3_cluster = t(mat_delta2[rownames(mat_delta2) %in% cluster_cell_ind,] ) %*% mat_gamma2[rownames(mat_gamma2) %in% cluster_cell_ind,]
  count_gdT_cdr3_cluster = table(mat_gdT_cdr3_cluster)
  if (length(count_gdT_cdr3_cluster) > 3){
    mat_cluster_clonotype_counts[,(i+1)] = c(count_gdT_cdr3_cluster[names(count_gdT_cdr3_cluster) == 1],
                                             count_gdT_cdr3_cluster[names(count_gdT_cdr3_cluster) == 2],
                                             sum(count_gdT_cdr3_cluster[!(names(count_gdT_cdr3_cluster) %in% c(0,1,2))]))
  } else if(length(count_gdT_cdr3_cluster) > 2){
    mat_cluster_clonotype_counts[,(i+1)] = c(count_gdT_cdr3_cluster[names(count_gdT_cdr3_cluster) == 1],
                                             count_gdT_cdr3_cluster[names(count_gdT_cdr3_cluster) == 2],
                                             0)
  } else if(length(count_gdT_cdr3_cluster) > 1){
    mat_cluster_clonotype_counts[,(i+1)] = c(count_gdT_cdr3_cluster[names(count_gdT_cdr3_cluster) == 1],
                                             0,
                                             0)
  }
}

colnames(mat_cluster_clonotype_counts) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
df_long= melt(mat_cluster_clonotype_counts)
colnames(df_long) = c("cell","cluster","counts")


ggplot(df_long)+
  geom_bar(stat ="identity", position="fill", aes(x = cluster, y= counts, fill= cell)) +
  scale_fill_manual(values=c("#e0f3db", "#a8ddb5", "#43a2ca")) +
  xlab("Cluster") +
  ylab("Clonotype type frequency") +
  scale_y_continuous(labels = scales::percent) +
  annotate(x=0, xend=0, y=0, yend=1, lwd=0.75, geom="segment") +
  theme(
    panel.grid = element_blank(), 
    axis.line = element_blank(),
    axis.ticks.x =element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y=element_text(size=20),
    legend.title = element_blank(),
    legend.text = element_text(size=20)
  ) + ggsave("Figure2_gdT_clontype_barplot.pdf")
dev.off()


mat_cdr3_pair_gdT2 = mat_cdr3_pair_gdT[use_one_clonotype_ind,]
df_clonotype = data.frame(pair =colnames(mat_cdr3_pair_gdT2), clonotype = paste0("clonotype",seq(1,ncol(mat_cdr3_pair_gdT2))))
View(df_clonotype)

mat_cluster_clonotypes = matrix(data=0, ncol=nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident),
                                nrow =ncol(mat_cdr3_pair_gdT2))
colnames(mat_cluster_clonotypes) = levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)
rownames(mat_cluster_clonotypes) = colnames(mat_cdr3_pair_gdT2)
for(i in 0:(nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)-1)){
  cluster_cell_ind = names(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)[seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident == i]
  mat_cdr3_pair_gdT2_cluster = mat_cdr3_pair_gdT2[rownames(mat_cdr3_pair_gdT2) %in% cluster_cell_ind,]
  mat_cluster_clonotypes[,(i+1)] = colSums(mat_cdr3_pair_gdT2_cluster)
}
rownames(mat_cluster_clonotypes) = df_clonotype$clonotype

colnames(mat_cluster_clonotypes) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))

i=1
j=2




# geo mean
# A, B cluster index
# sqrt( sum( A ) sum( B) )
geo_mean = matrix(data = 0, ncol = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident),
                  nrow = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
rownames(geo_mean) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
colnames(geo_mean) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))


for(i in 1:nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)){
  for(j in 1:nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)){
    if (i != j){
      sumA = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),i]/sum(mat_cluster_clonotypes[,i]))
      sumB = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),j]/sum(mat_cluster_clonotypes[,j]))
      geo_mean[i,j] = sqrt(sumA*sumB)
    }
  }
}


pheatmap(geo_mean, color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                             "RdBu")))(100),
         cluster_rows = F, cluster_cols = F, filename = "Figure2_gdT_clonotype_cluster_heatmap_geo_mean.pdf" )



# P value 
# 10000 sampling
set.seed(42)
mat_pvalue = matrix(data = 0, ncol = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident),
                    nrow = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
rownames(mat_pvalue) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
colnames(mat_pvalue) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
mat_accum = matrix(data = 0, ncol = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident),
                   nrow = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
rownames(mat_accum) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
colnames(mat_accum) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))

totaliter =10000
for(iter in 1:totaliter){
  sampled_cluster = seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@cell.names
  names(sampled_cluster) = seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@cell.names
  sampled_cluster =sampled_cluster[sample(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@cell.names)]
  
  cluster_counts = table(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)
  sampled_cluster_value = c()
  for(i in 1:(nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))){
    sampled_cluster_value = c(sampled_cluster_value, rep(as.numeric(names(cluster_counts)[i]),cluster_counts[[i]]))
  }  
  names(sampled_cluster_value) = names(sampled_cluster)
  sampled_cluster = sampled_cluster_value
  
  sampled_mat_cluster_clonotypes = matrix(data=0, ncol=nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident),
                                          nrow =ncol(mat_cdr3_pair_gdT2))
  colnames(sampled_mat_cluster_clonotypes) = levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)
  rownames(sampled_mat_cluster_clonotypes) = colnames(mat_cdr3_pair_gdT2)
  
  for(i in 0:(nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)-1)){
    cluster_cell_ind = names(sampled_cluster)[sampled_cluster == i]
    mat_cdr3_pair_gdT2_cluster = mat_cdr3_pair_gdT2[rownames(mat_cdr3_pair_gdT2) %in% cluster_cell_ind,]
    sampled_mat_cluster_clonotypes[,(i+1)] = colSums(mat_cdr3_pair_gdT2_cluster)
  }
  
  rownames(sampled_mat_cluster_clonotypes) = df_clonotype$clonotype
  
  colnames(sampled_mat_cluster_clonotypes) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
  
  sampled_geo_mean = matrix(data = 0, ncol = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident),
                            nrow = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
  rownames(sampled_geo_mean) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
  colnames(sampled_geo_mean) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
  
  for(i in 1:nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)){
    for(j in 1:nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)){
      if (i != j){
        sumA = sum(sampled_mat_cluster_clonotypes[(sampled_mat_cluster_clonotypes[,i]*sampled_mat_cluster_clonotypes[,j] > 0),i]/sum(sampled_mat_cluster_clonotypes[,i]))
        sumB = sum(sampled_mat_cluster_clonotypes[(sampled_mat_cluster_clonotypes[,i]*sampled_mat_cluster_clonotypes[,j] > 0),j]/sum(sampled_mat_cluster_clonotypes[,j]))
        sampled_geo_mean[i,j] = sqrt(sumA*sumB)
      }
    }
  }
  for(i in 1:nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)){
    for(j in 1:nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)){
      if (i != j){
        mat_accum[i,j] = mat_accum[i,j] + sampled_geo_mean[i,j]
        if(sampled_geo_mean[i,j] >= geo_mean[i,j])
          mat_pvalue[i,j] = mat_pvalue[i,j] + 1
      }
    }
  }
}
P = (mat_pvalue +1) /(totaliter+1)
mat_pvalue_twotail = matrix(data = 0, ncol = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident),
                            nrow = nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
rownames(mat_pvalue_twotail) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
colnames(mat_pvalue_twotail) = paste0("G",levels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident))
for(i in 1:nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)){
  for(j in 1:nlevels(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)){
    mat_pvalue_twotail[i,j] = 2*min(P[i,j],(1-P[i,j]))
  }
}
mat_pvalue_twotail
P
P < 0.05
P < 0.01
mat_avg_geo_mean = mat_accum/10000

mat_log2FC = log2((geo_mean+1)/(mat_avg_geo_mean+1))

paletteLength <- 100
myColor <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(mat_log2FC), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat_log2FC)/paletteLength, max(mat_log2FC), length.out=floor(paletteLength/2)))
pheatmap(mat_log2FC, cluster_rows = F, cluster_cols = F, color = myColor, breaks = myBreaks, filename = "Figure2_gdT_clonotype_cluster_heatmap_logFC.pdf" )
pheatmap(P, cluster_rows = F, cluster_cols = F)
# table(seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT@ident)
mean_0 <- mean(nchar(as.character(all_vdj_prod_highconf_subset_TRB$cdr3[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==0]])))
mean_1 <- mean(nchar(as.character(all_vdj_prod_highconf_subset_TRB$cdr3[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==1]])))
mean_2 <- mean(nchar(as.character(all_vdj_prod_highconf_subset_TRB$cdr3[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==2]])))
mean_3 <- mean(nchar(as.character(all_vdj_prod_highconf_subset_TRB$cdr3[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==3]])))
mean_4 <- mean(nchar(as.character(all_vdj_prod_highconf_subset_TRB$cdr3[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==4]])))
mean_5 <- mean(nchar(as.character(all_vdj_prod_highconf_subset_TRB$cdr3[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==5]])))

plot(x=c(0,1,2,3,4,5),y=c(mean_0,mean_1,mean_2,mean_3,mean_4,mean_5),type = "l")

table(all_vdj_prod_highconf_subset_TRB$v_gene[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==1]] %in% all_vdj_prod_highconf_subset_TRB$v_gene[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==2]])[1]
mat_per <-matrix(nrow =6, ncol = 6 )
for (i in levels(NKT_subset@active.ident)){
  print(i)
  
  for (j in levels(NKT_subset@active.ident)){
    print(j)
    mat_per[as.numeric(i)+1,as.numeric(j)+1] <- as.numeric(table(all_vdj_prod_highconf_subset_TRA$v_gene[all_vdj_prod_highconf_subset_TRA$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==j]] %in% all_vdj_prod_highconf_subset_TRA$v_gene[all_vdj_prod_highconf_subset_TRA$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==i]])[2])/(as.numeric(table(all_vdj_prod_highconf_subset_TRA$v_gene[all_vdj_prod_highconf_subset_TRA$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==j]] %in% all_vdj_prod_highconf_subset_TRA$v_gene[all_vdj_prod_highconf_subset_TRA$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==i]]))[1]+as.numeric(table(all_vdj_prod_highconf_subset_TRA$v_gene[all_vdj_prod_highconf_subset_TRA$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==j]] %in% all_vdj_prod_highconf_subset_TRA$v_gene[all_vdj_prod_highconf_subset_TRA$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==i]]))[2])
    print(mat_per[as.numeric(i),as.numeric(j)])
    }
}


mat_per[is.na(mat_per)] <-1
c0 <- mat_per[,1]
c1 <- mat_per[,2]
c2 <- mat_per[,3]
c3 <- mat_per[,4]
c4 <- mat_per[,5]
c5 <- mat_per[,6]

df_per <- data.frame(C0=c0,C1=c1,C2=c2,C3=c3,C4=c4,C5=c5)
barplot(as.matrix(df_per), main="TCR_Alpha_chain_overlap",beside = T, col = rainbow(nrow(df_per)),ylim=c(0,2))
legend(30,2,c("cluster0","cluster1", "cluster2","cluster3","cluster4","cluster5"),fill = rainbow(nrow(df_per)))




mat_per <-matrix(nrow =6, ncol = 6 )
for (i in levels(NKT_subset@active.ident)){
  print(i)
  
  for (j in levels(NKT_subset@active.ident)){
    print(j)
    mat_per[as.numeric(i)+1,as.numeric(j)+1] <- as.numeric(table(all_vdj_prod_highconf_subset_TRB$v_gene[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==j]] %in% all_vdj_prod_highconf_subset_TRB$v_gene[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==i]])[2])/(as.numeric(table(all_vdj_prod_highconf_subset_TRB$v_gene[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==j]] %in% all_vdj_prod_highconf_subset_TRB$v_gene[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==i]]))[1]+as.numeric(table(all_vdj_prod_highconf_subset_TRB$v_gene[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==j]] %in% all_vdj_prod_highconf_subset_TRB$v_gene[all_vdj_prod_highconf_subset_TRB$cellname %in% colnames(NKT_subset)[NKT_subset@active.ident==i]]))[2])
    print(mat_per[as.numeric(i),as.numeric(j)])
  }
}


mat_per[is.na(mat_per)] <-1
c0 <- mat_per[,1]
c1 <- mat_per[,2]
c2 <- mat_per[,3]
c3 <- mat_per[,4]
c4 <- mat_per[,5]
c5 <- mat_per[,6]

df_per <- data.frame(C0=c0,C1=c1,C2=c2,C3=c3,C4=c4,C5=c5)
barplot(as.matrix(df_per), main="TCR_Beta_chain_overlap",beside = T, col = rainbow(nrow(df_per)),ylim=c(0,2))
legend(30,2,c("cluster0","cluster1", "cluster2","cluster3","cluster4","cluster5"),fill = rainbow(nrow(df_per)))

mat_per <- NKT_subset@active.ident

barplot(100*mat_per[,6])
