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
df_chain1 = subset(gdT_cdr3_mat[[2]], chain %in% c("TRG"))
df_chain1_cdr3_len = data.frame(row.names = df_chain1$cellName, cdr3 = df_chain1$cdr3,
                                v_gene = df_chain1$v_gene)
df_chain1_cdr3_len$v_gene = gsub(".*(GV[0-9]+).*","\\1",as.vector(df_chain1_cdr3_len$v_gene))
df_chain1_cdr3_len$cdr3Length = ""
for( i in 1:nrow(df_chain1_cdr3_len)){
  df_chain1_cdr3_len[i,"cdr3Length"] = nchar(sset_gdT@contig[sset_gdT@contig$cdr3 == as.vector(df_chain1_cdr3_len[i,"cdr3"]),"cdr3_nt"])[1]
}

breaks= unique(sort(c("30","33","36","39","42","45","48","51","54","57","60","63","66","69","72",
                      "30","33","36","39","42","45","48","51","54","57",
                      "15","21","24","27","30","33","36","39","42","45","48","51","54","57","60","69")))

ggplot(df_chain1_cdr3_len) + 
  geom_bar(aes(x=cdr3Length, fill = v_gene)) + 
  xlim(breaks) +
  ggsave("gdT_gamma_vgene_cdr3Length.pdf", width=10, height=7)
dev.off()




library(dplyr)
d3 = df_chain1_cdr3_len %>% group_by(cdr3Length,v_gene) %>% summarise(count=n()) %>% as.data.frame(.)
d3$percent = apply(d3,1,function(x){as.numeric(x[3])/table(df_chain1_cdr3_len$v_gene)[x[2]]})
write.csv(d3,file="gdT_cdrLen.csv")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
gcolor = gg_color_hue(nlevels(factor(d3$v_gene)))
cnt = 0
for(i in levels(factor(d3$v_gene))){
  cnt = cnt + 1
  ggplot(subset(d3, v_gene == i)) + 
    geom_bar(aes(x = factor(cdr3Length), y = percent*100, fill = v_gene),stat="identity", width = 0.7) + 
    scale_fill_manual(values=c(gcolor[cnt])) +
    xlim(breaks) + 
    ggsave(paste0("gdT_gamma_vgene_cdr3Length",i,".pdf"), width=10, height=7)
  
}
ggplot(subset(d3, v_gene == "GV1")) + 
  geom_bar(aes(x = factor(cdr3Length), y = percent*100, fill = v_gene),stat="identity", width = 0.7) + 
  xlim(breaks) + 
  ggsave("gdT_gamma_vgene_cdr3Length.pdf", width=10, height=7)
dev.off()

d2 <- df_chain1_cdr3_len %>% 
  group_by(cdr3Length,v_gene) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

ggplot(d2) + 
  geom_bar(aes(x = factor(cdr3Length), y = perc*100, fill = v_gene),stat="identity", width = 0.7) + 
  xlim(breaks) +
  theme_bw() + 
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    plot.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), 
    axis.line = element_blank(),
  )
ggsave("gdT_gamma_vgene_cdr3Length_percentage.pdf", width=10, height=7)

write.csv(df_chain1_cdr3_len, file="gdT_gamma_vgene_cdr3Length_count.csv")

# delta
df_chain2 = subset(gdT_cdr3_mat[[2]], chain %in% c("TRD","Multi"))
df_chain2_cdr3_len = data.frame(row.names = df_chain2$cellName, cdr3 = df_chain2$cdr3,
                                v_gene = df_chain2$v_gene)
df_chain2_cdr3_len$v_gene = gsub(".*(DV[0-9]+).*","\\1",as.vector(df_chain2_cdr3_len$v_gene))
df_chain2_cdr3_len$v_gene = factor(df_chain2_cdr3_len$v_gene,levels=c( "DV1", "DV2","DV4","DV5","DV6","DV7","DV9","DV10", "DV12" ))
df_chain2_cdr3_len$cdr3Length = ""
for( i in 1:nrow(df_chain2_cdr3_len)){
  df_chain2_cdr3_len[i,"cdr3Length"] = nchar(sset_gdT@contig[sset_gdT@contig$cdr3 == as.vector(df_chain2_cdr3_len[i,"cdr3"]),"cdr3_nt"])[1]
}

ggplot(df_chain2_cdr3_len) + 
  geom_bar(aes(x=cdr3Length, fill = v_gene)) + 
  xlim(breaks) +
  ggsave("gdT_delta_vgene_cdr3Length.pdf", width=10, height=7)
dev.off()


library(dplyr)
d3 = df_chain2_cdr3_len %>% group_by(cdr3Length,v_gene) %>% summarise(count=n()) %>% as.data.frame(.)
d3$percent = apply(d3,1,function(x){as.numeric(x[3])/table(df_chain2_cdr3_len$v_gene)[x[2]]})
write.csv(d3,file="gdT_delta_cdrLen.csv")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
gcolor = gg_color_hue(nlevels(factor(d3$v_gene)))
cnt = 0
for(i in levels(factor(d3$v_gene))){
  cnt = cnt + 1
  ggplot(subset(d3, v_gene == i)) + 
    geom_bar(aes(x = factor(cdr3Length), y = percent*100, fill = v_gene),stat="identity", width = 0.7) + 
    scale_fill_manual(values=c(gcolor[cnt])) +
    xlim(breaks) + 
    ggsave(paste0("gdT_delta_vgene_cdr3Length",i,".pdf"), width=10, height=7)
  
}


d2 <- df_chain2_cdr3_len %>% 
  group_by(cdr3Length,v_gene) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

ggplot(d2) + 
  geom_bar(aes(x = factor(cdr3Length), y = perc*100, fill = v_gene),stat="identity", width = 0.7) + 
  xlim(breaks) +
  theme_bw() + 
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    plot.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), 
    axis.line = element_blank(),
  )
ggsave("gdT_delta_vgene_cdr3Length_percentage.pdf", width=10, height=7)

write.csv(df_chain2_cdr3_len, file="gdT_delta_vgene_cdr3Length_count.csv")
