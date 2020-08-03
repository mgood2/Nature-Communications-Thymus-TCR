# init
setwd("D:/Dropbox/R projects/Stcr/tests/testthat")
load("../testdata/seurat.gdT.RData") # seurat v2 file
load("../testdata/ident_gdT.RData")
ssetF3 = Load10xVDJ('../testdata/mTCR-F12-gdT/out_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-G1-gdT/out_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5)

sset_gdT = ImportSeurat(sset,seurat.gdT)

gdT_cdr3_mat = SetCDR3Mat(sset_gdT, cell.ind = colnames(seurat.gdT@scale.data), chain1 = c("TRG"), chain2 = c("TRD","Multi"), v_gene2 = "DV",
                          j_gene2="DJ")



gdT_cdr3_pair_mat = SetCDR3PairMat(sset_gdT, gdT_cdr3_mat, chain1 = c("TRG"), chain2 = c("TRD","Multi"))

# gamma
df_chain1 = subset(gdT_cdr3_mat[[2]], chain %in% c("TRG") & cellName %in% colnames(seurat.gdT@scale.data))
df_chain1_cdr3_amino_len = data.frame(row.names = df_chain1$cellName, cdr3 = df_chain1$cdr3,
                                v_gene = df_chain1$v_gene)
df_chain1_cdr3_amino_len$v_gene = gsub(".*(GV[0-9]+).*","\\1",as.vector(df_chain1_cdr3_amino_len$v_gene))
df_chain1_cdr3_amino_len$cdr3Length = ""
for( i in 1:nrow(df_chain1_cdr3_amino_len)){
  df_chain1_cdr3_amino_len[i,"cdr3Length"] = nchar(sset_gdT@contig[sset_gdT@contig$cdr3 == as.vector(df_chain1_cdr3_amino_len[i,"cdr3"]),"cdr3"])[1]
}

write.csv(df_chain1_cdr3_amino_len, file="gdT_gamma_vgene_cdr3Length_amino_count.csv")

for(i in levels(factor(df_chain1_cdr3_amino_len$v_gene))){
  targetGene = i
  printer1<- file(paste0("weblogo_gdT_gamma_",targetGene,".txt"),"w")
  targetdf = subset(df_chain1_cdr3_amino_len, v_gene == targetGene)
  countdf = data.frame(counts = table(factor(targetdf$cdr3Length)))
  write.csv(countdf, file= paste0("gdT_gamma_vgene_",targetGene,"_amino_cdr3Length.csv"))
  ggplot(countdf) + 
    geom_point(aes(x=factor(counts.Var1, levels= sort(as.numeric(levels(factor(countdf$counts.Var1))))), y=counts.Freq)) + 
    geom_line(aes(x=factor(counts.Var1, levels= sort(as.numeric(levels(factor(countdf$counts.Var1))))), y=counts.Freq),group = 1) +
    # xlim(breaks) +
    ggsave(paste0("gdT_gamma_vgene_",targetGene,"_amino_cdr3Length.pdf"), width=10, height=7)
  targetLen = names(which.max(table(targetdf$cdr3Length)))
  maxtargetdf = subset(df_chain1_cdr3_amino_len, v_gene == targetGene & cdr3Length == targetLen)
  
  for(i in 1:nrow(maxtargetdf)){
    writeLines(paste0(">",maxtargetdf$v_gene[i]),con=printer1,sep="\n")
    writeLines(as.vector(maxtargetdf$cdr3)[i],con=printer1, sep="\n")
  }
  close(printer1)
}


#weblogo options
# weblogo < weblogo_gdT_gamma_GV2.txt > weblogo_gdT_gamma_GV2.png --format "png" --errorbars FALSE --units "probability" --show-xaxis FALSE --show-yaxis FALSE



# delta
df_chain2 = subset(gdT_cdr3_mat[[2]], chain %in% c("TRD","Multi") & cellName %in% colnames(seurat.gdT@scale.data))
df_chain2_cdr3_amino_len = data.frame(row.names = df_chain2$cellName, cdr3 = df_chain2$cdr3,
                                      v_gene = df_chain2$v_gene)
df_chain2_cdr3_amino_len$v_gene = gsub(".*(DV[0-9]+).*","\\1",as.vector(df_chain2_cdr3_amino_len$v_gene))
df_chain2_cdr3_amino_len$cdr3Length = ""
for( i in 1:nrow(df_chain2_cdr3_amino_len)){
  df_chain2_cdr3_amino_len[i,"cdr3Length"] = nchar(sset_gdT@contig[sset_gdT@contig$cdr3 == as.vector(df_chain2_cdr3_amino_len[i,"cdr3"]),"cdr3"])[1]
}

write.csv(df_chain2_cdr3_amino_len, file="gdT_delta_vgene_cdr3Length_amino_count.csv")

for(i in levels(factor(df_chain2_cdr3_amino_len$v_gene))){
  targetGene = i
  printer1<- file(paste0("weblogo_gdT_delta_",targetGene,".txt"),"w")
  targetdf = subset(df_chain2_cdr3_amino_len, v_gene == targetGene)
  countdf = data.frame(counts = table(factor(targetdf$cdr3Length)))
  write.csv(countdf, file= paste0("gdT_delta_vgene_",targetGene,"_amino_cdr3Length.csv"))
  ggplot(countdf) + 
    geom_point(aes(x=factor(counts.Var1, levels= sort(as.numeric(levels(factor(countdf$counts.Var1))))), y=counts.Freq)) + 
    geom_line(aes(x=factor(counts.Var1, levels= sort(as.numeric(levels(factor(countdf$counts.Var1))))), y=counts.Freq),group = 1) +
    # xlim(breaks) +
    ggsave(paste0("gdT_delta_vgene_",targetGene,"_amino_cdr3Length.pdf"), width=10, height=7)
  targetLen = names(which.max(table(targetdf$cdr3Length)))
  maxtargetdf = subset(df_chain2_cdr3_amino_len, v_gene == targetGene & cdr3Length == targetLen)
  
  for(i in 1:nrow(maxtargetdf)){
    writeLines(paste0(">",maxtargetdf$v_gene[i]),con=printer1,sep="\n")
    writeLines(as.vector(maxtargetdf$cdr3)[i],con=printer1, sep="\n")
  }
  close(printer1)
}


#weblogo options
# weblogo < weblogo_gdT_delta_DV1.txt > weblogo_gdT_delta_DV1.png --format "png" --errorbars FALSE --units "probability" --show-xaxis FALSE --show-yaxis FALSE

