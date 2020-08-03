# init
library(ggplot2)
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


df_clonotype_noncano = subset(df_clonotype_labeled, selected == FALSE)

df = data.frame(x = seurat.MAIT@dr$umap@cell.embeddings[,1],
                y = seurat.MAIT@dr$umap@cell.embeddings[,2],
                cano = colnames(seurat.MAIT@data) %in% selected_ind,
                noncano = colnames(seurat.MAIT@data) %in% nonselected_ind,
                TCR = colnames(seurat.MAIT@data) %in% selected_ind + colnames(seurat.MAIT@data) %in% nonselected_ind,
                cluster = seurat.MAIT@ident)

TRAV_MAIT = data.frame(subset(df_clonotype_labeled, selected == FALSE )$chain1_v_gene, 
           subset(df_clonotype_labeled, selected == FALSE )$n)
TRAJ_MAIT = data.frame(subset(df_clonotype_labeled, selected == FALSE )$chain1_j_gene, 
                       subset(df_clonotype_labeled, selected == FALSE )$n)

library(dplyr)

TRAV_list = unique(as.vector(TRAV_MAIT[,1]))
TRAV_list = gsub("(TRAV[0-9D-]+)(.*)","\\1",TRAV_list)
TRAV_counter = data.frame(row.names = unique(TRAV_list))
TRAV_counter$counts = 0
for(j in 1:nrow(TRAV_MAIT)){
  for(i in 1:nrow(TRAV_counter)){
    # if(as.vector(TRAV_MAIT[j,1]) == rownames(TRAV_counter)[i]){
    if(grepl(paste0(rownames(TRAV_counter)[i],"[^0-9D-].*"), as.vector(TRAV_MAIT[j,1]))){
      TRAV_counter[i,1] = TRAV_counter[i,1] + TRAV_MAIT[j,2]
    }
  }
}

TRAV_counter$TRAV = rownames(TRAV_counter)
TRAV_counter = TRAV_counter[order(TRAV_counter$counts, decreasing = T),]
TRAV_counter$TRAV = factor(TRAV_counter$TRAV, levels = rownames(TRAV_counter))

# remove canonical tcr
TRAV_counter = TRAV_counter[!(rownames(TRAV_counter) %in% c("TRAV1")),]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycolors = gg_color_hue(nlevels(TRAV_counter$TRAV))
length(mycolors)
color_seq = c()
for(i in 1:(length(mycolors)/2)){
  color_seq = c(color_seq,seq(1,((length(mycolors)/2)),1)[i], (seq((length(mycolors)/2) +1, length(mycolors),1))[i])
}

ggplot(TRAV_counter, aes(x="", y=counts, fill= TRAV)) + 
  scale_fill_manual(values = mycolors[color_seq]) +
  geom_bar(stat="identity") + coord_polar("y", start =0) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title=element_text(size=14, face="bold"),
  ) + ggsave("pi chart MAIT TRAV.pdf", width =7, height= 7)

ggplot(TRAV_counter, aes(x=TRAV, y=counts), color="grey") + 
  # scale_fill_manual(values = mycolors[color_seq]) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    legend.position = "none",
    # axis.ticks = element_blank(),
    axis.text.y = element_text(size=30),
    axis.text.x = element_text(angle=90, hjust = 0.0, vjust = 0.5, size=30),
    plot.title=element_text(size=14, face="bold"),
  )+ ggsave("bar chart MAIT noncano TRAV 2.pdf", width =22, height= 7)



TRAJ_list = unique(as.vector(TRAJ_MAIT[,1]))
TRAJ_list = gsub("(TRAJ[0-9D-]+)(.*)","\\1",TRAJ_list)
TRAJ_counter = data.frame(row.names = unique(TRAJ_list))
TRAJ_counter$counts = 0

for(j in 1:nrow(TRAJ_MAIT)){
  for(i in 1:nrow(TRAJ_counter)){
    # if(as.vector(TRAV_MAIT[j,1]) == rownames(TRAV_counter)[i]){
    if(grepl(paste0(rownames(TRAJ_counter)[i],"[^0-9D-].*"), as.vector(TRAJ_MAIT[j,1]))){
      TRAJ_counter[i,1] = TRAJ_counter[i,1] + TRAJ_MAIT[j,2]
    }
  }
}

TRAJ_counter$TRAJ = rownames(TRAJ_counter)
TRAJ_counter = TRAJ_counter[order(TRAJ_counter$counts, decreasing = T),]
TRAJ_counter$TRAJ = factor(TRAJ_counter$TRAJ, levels = rownames(TRAJ_counter))


# remove canonical tcr
TRAJ_counter = TRAJ_counter[!(rownames(TRAJ_counter) %in% c("TRAJ33")),]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycolors = gg_color_hue(nlevels(TRAJ_counter$TRAJ))
length(mycolors)
color_seq = c()
for(i in 1:(length(mycolors)/2)){
  color_seq = c(color_seq,seq(1,((length(mycolors)/2)),1)[i], (seq((length(mycolors)/2) +1, length(mycolors),1))[i])
}

ggplot(TRAJ_counter, aes(x="", y=counts, fill= TRAJ)) + 
  scale_fill_manual(values = mycolors[color_seq]) +
  geom_bar(stat="identity") + coord_polar("y", start =0) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title=element_text(size=14, face="bold"),
  )+ ggsave("pi chart MAIT TRAJ.pdf", width =7, height= 7)

ggplot(TRAJ_counter, aes(x=TRAJ, y=counts), color="grey") + 
  # scale_fill_manual(values = mycolors[color_seq]) +
  geom_bar(stat="identity") + #coord_polar("y", start =0) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    legend.position = "none",
    # axis.ticks = element_blank(),
    axis.text.y = element_text(size=30),
    axis.text.x = element_text(angle=90, hjust = 0.1, vjust = 0.5, size=30),
    plot.title=element_text(size=14, face="bold"),
  )+ ggsave("bar chart MAIT noncano TRAJ 2.pdf", width =17, height= 7)

dev.off()
