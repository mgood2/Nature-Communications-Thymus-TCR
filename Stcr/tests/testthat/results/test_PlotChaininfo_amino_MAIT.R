# init
setwd("D:/elee4472/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.MAIT.RData") # seurat v2 file
load("../testdata/ident_MAIT.RData")


ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_MAIT = ImportSeurat(sset,seurat.MAIT)

MAIT_cdr3_mat = SetCDR3Mat(sset_MAIT, cell.ind = colnames(seurat.MAIT@scale.data), chain1 = c("TRA"), chain2 = c("TRB"))


# alpha
df_chain1 = subset(MAIT_cdr3_mat[[2]], chain %in% c("TRA") & cellName %in% colnames(seurat.MAIT@scale.data))
df_chain1_cdr3_amino_len = data.frame(row.names = df_chain1$cellName, cdr3 = df_chain1$cdr3,
                                v_gene = df_chain1$v_gene)
df_chain1_cdr3_amino_len$v_gene = gsub(".*(AV[0-9]+).*","\\1",as.vector(df_chain1_cdr3_amino_len$v_gene))
df_chain1_cdr3_amino_len$canonical = "noncanonical"
df_chain1_cdr3_amino_len$canonical[df_chain1_cdr3_amino_len$v_gene %in% c("AV1")] = "canonical"
df_chain1_cdr3_amino_len$cdr3Length = ""
for( i in 1:nrow(df_chain1_cdr3_amino_len)){
  df_chain1_cdr3_amino_len[i,"cdr3Length"] = nchar(sset_MAIT@contig[sset_MAIT@contig$cdr3 == as.vector(df_chain1_cdr3_amino_len[i,"cdr3"]),"cdr3"])[1]
}
breaks= unique(sort(c("30","33","36","39","42","45","48","51","54","57","60","63","66","69","72",
                      "30","33","36","39","42","45","48","51","54","57",
                      "15","21","24","27","30","33","36","39","42","45","48","51","54","57","60","69")))

# ggplot(df_chain1_cdr3_amino_len) +
#   geom_bar(aes(x=cdr3Length, fill = canonical)) +
#   xlim(breaks) +
#   ggsave("MAIT_alpha_canonical_cdr3Length.pdf", width=10, height=7)
# dev.off()

write.csv(df_chain1_cdr3_amino_len, file="MAIT_alpha_canonical_cdr3Length_amino_count.csv")


for(i in levels(factor(df_chain1_cdr3_amino_len$canonical))){
  targetGene = i
  printer1<- file(paste0("weblogo_MAIT_alpha_",targetGene,".txt"),"w")
  targetdf = subset(df_chain1_cdr3_amino_len, canonical == targetGene)
  countdf = data.frame(counts = table(factor(targetdf$cdr3Length)))
  write.csv(countdf, file=paste0("MAIT_alpha_",targetGene,"_amino_cdr3Length.csv"))
  ggplot(countdf) + 
    geom_point(aes(x=factor(counts.Var1, levels= sort(as.numeric(levels(factor(countdf$counts.Var1))))), y=counts.Freq)) + 
    geom_line(aes(x=factor(counts.Var1, levels= sort(as.numeric(levels(factor(countdf$counts.Var1))))), y=counts.Freq),group = 1) +
    # xlim(breaks) +
    ggsave(paste0("MAIT_alpha_",targetGene,"_amino_cdr3Length.pdf"), width=10, height=7)
  targetLen = names(which.max(table(targetdf$cdr3Length)))
  maxtargetdf = subset(df_chain1_cdr3_amino_len, canonical == targetGene & cdr3Length == targetLen)
  
  for(i in 1:nrow(maxtargetdf)){
    writeLines(paste0(">",maxtargetdf$canonical[i]),con=printer1,sep="\n")
    writeLines(as.vector(maxtargetdf$cdr3)[i],con=printer1, sep="\n")
  }
  close(printer1)
}


#weblogo options
# weblogo < weblogo_gdT_gamma_GV2.txt > weblogo_gdT_gamma_GV2.png --format "png" --errorbars FALSE --units "probability" --show-xaxis FALSE --show-yaxis FALSE




# beta
df_chain2 = subset(MAIT_cdr3_mat[[2]], chain %in% c("TRB") & cellName %in% colnames(seurat.MAIT@scale.data))
df_chain2_cdr3_amino_len = data.frame(row.names = df_chain2$cellName, cdr3 = df_chain2$cdr3,
                                v_gene = df_chain2$v_gene)
df_chain2_cdr3_amino_len$v_gene = gsub(".*(BV[0-9]+).*","\\1",as.vector(df_chain2_cdr3_amino_len$v_gene))

df_chain2_cdr3_amino_len$canonical = "noncanonical"
df_chain2_cdr3_amino_len$canonical[df_chain2_cdr3_amino_len$v_gene %in% c("BV13","BV19")] = "canonical"


df_chain2_cdr3_amino_len$cdr3Length = ""
for( i in 1:nrow(df_chain2_cdr3_amino_len)){
  df_chain2_cdr3_amino_len[i,"cdr3Length"] = nchar(sset_MAIT@contig[sset_MAIT@contig$cdr3 == as.vector(df_chain2_cdr3_amino_len[i,"cdr3"]),"cdr3"])[1]
}

write.csv(df_chain2_cdr3_amino_len, file="MAIT_beta_canonical_cdr3Length_amino_count.csv")

for(i in levels(factor(df_chain2_cdr3_amino_len$canonical))){
  targetGene = i
  printer1<- file(paste0("weblogo_MAIT_beta_",targetGene,".txt"),"w")
  targetdf = subset(df_chain2_cdr3_amino_len, canonical == targetGene)
  countdf = data.frame(counts = table(factor(targetdf$cdr3Length)))
  write.csv(countdf, file=paste0("MAIT_beta_",targetGene,"_amino_cdr3Length.csv"))
  ggplot(countdf) + 
    geom_point(aes(x=factor(counts.Var1, levels= sort(as.numeric(levels(factor(countdf$counts.Var1))))), y=counts.Freq)) + 
    geom_line(aes(x=factor(counts.Var1, levels= sort(as.numeric(levels(factor(countdf$counts.Var1))))), y=counts.Freq),group = 1) +
    # xlim(breaks) +
    ggsave(paste0("MAIT_beta_",targetGene,"_amino_cdr3Length.pdf"), width=10, height=7)
  targetLen = names(which.max(table(targetdf$cdr3Length)))
  maxtargetdf = subset(df_chain2_cdr3_amino_len, canonical == targetGene & cdr3Length == targetLen)
  
  for(i in 1:nrow(maxtargetdf)){
    writeLines(paste0(">",maxtargetdf$canonical[i]),con=printer1,sep="\n")
    writeLines(as.vector(maxtargetdf$cdr3)[i],con=printer1, sep="\n")
  }
  close(printer1)
}



#weblogo options
# weblogo < weblogo_gdT_gamma_GV2.txt > weblogo_gdT_gamma_GV2.png --format "png" --errorbars FALSE --units "probability" --show-xaxis FALSE --show-yaxis FALSE
