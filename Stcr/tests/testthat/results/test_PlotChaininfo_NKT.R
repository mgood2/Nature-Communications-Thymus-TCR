# init
setwd("D:/Dropbox/R projects/Stcr/tests/testthat")
load("../testdata/seurat.NKT.RData") # seurat v2 file
load("../testdata/ident_NKT.RData")


ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_NKT = ImportSeurat(sset,seurat.NKT)

NKT_cdr3_mat = SetCDR3Mat(sset_NKT, cell.ind = colnames(seurat.NKT@scale.data), chain1 = c("TRA"), chain2 = c("TRB"))


# alpha
df_chain1 = subset(NKT_cdr3_mat[[2]], chain %in% c("TRA"))
df_chain1_cdr3_len = data.frame(row.names = df_chain1$cellName, cdr3 = df_chain1$cdr3,
                                v_gene = df_chain1$v_gene)
df_chain1_cdr3_len$v_gene = gsub(".*(AV[0-9]+).*","\\1",as.vector(df_chain1_cdr3_len$v_gene))
df_chain1_cdr3_len$canonical = "noncanonical"
df_chain1_cdr3_len$canonical[df_chain1_cdr3_len$v_gene %in% c("AV11")] = "canonical"
df_chain1_cdr3_len$cdr3Length = ""
for( i in 1:nrow(df_chain1_cdr3_len)){
  df_chain1_cdr3_len[i,"cdr3Length"] = nchar(sset_NKT@contig[sset_NKT@contig$cdr3 == as.vector(df_chain1_cdr3_len[i,"cdr3"]),"cdr3_nt"])[1]
}
breaks= unique(sort(c("30","33","36","39","42","45","48","51","54","57","60","63","66","69","72",
                      "30","33","36","39","42","45","48","51","54","57",
                      "15","21","24","27","30","33","36","39","42","45","48","51","54","57","60","69")))

ggplot(df_chain1_cdr3_len) +
  geom_bar(aes(x=cdr3Length, fill = canonical)) + 
  xlim(breaks) +
  ggsave("NKT_alpha_canonical_cdr3Length.pdf", width=10, height=7)
dev.off()

library(dplyr)
d3 = df_chain1_cdr3_len %>% group_by(cdr3Length,canonical) %>% summarise(count=n()) %>% as.data.frame(.)
d3$percent = apply(d3,1,function(x){as.numeric(x[3])/table(df_chain1_cdr3_len$canonical)[x[2]]})
write.csv(d3,file="NKT_alpha_cdrLen.csv")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
gcolor = gg_color_hue(nlevels(factor(d3$canonical)))
cnt = 0
for(i in levels(factor(d3$canonical))){
  cnt = cnt + 1
  ggplot(subset(d3, canonical == i)) + 
    geom_bar(aes(x = factor(cdr3Length), y = percent*100, fill = canonical),stat="identity", width = 0.7) + 
    scale_fill_manual(values=c(gcolor[cnt])) +
    xlim(breaks) + 
    ggsave(paste0("NKT_alpha_cdr3Length",i,".pdf"), width=10, height=7)
}
d2 <- df_chain1_cdr3_len %>% 
  group_by(cdr3Length,canonical) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

ggplot(d2) + 
  geom_bar(aes(x = factor(cdr3Length), y = perc*100, fill = canonical),stat="identity", width = 0.7) + 
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
  ) +
  ggsave("NKT_alpha_canonical_cdr3Length_percentage.pdf", width=10, height=7)

write.csv(df_chain1_cdr3_len, file="NKT_alpha_canonical_cdr3Length_count.csv")

# beta
df_chain2 = subset(NKT_cdr3_mat[[2]], chain %in% c("TRB"))
df_chain2_cdr3_len = data.frame(row.names = df_chain2$cellName, cdr3 = df_chain2$cdr3,
                                v_gene = df_chain2$v_gene)
df_chain2_cdr3_len$v_gene = gsub(".*(BV[0-9]+).*","\\1",as.vector(df_chain2_cdr3_len$v_gene))

df_chain2_cdr3_len$canonical = "noncanonical"
df_chain2_cdr3_len$canonical[df_chain2_cdr3_len$v_gene %in% c("BV1","BV13","BV29")] = "canonical"


df_chain2_cdr3_len$cdr3Length = ""
for( i in 1:nrow(df_chain2_cdr3_len)){
  df_chain2_cdr3_len[i,"cdr3Length"] = nchar(sset_NKT@contig[sset_NKT@contig$cdr3 == as.vector(df_chain2_cdr3_len[i,"cdr3"]),"cdr3_nt"])[1]
}




ggplot(df_chain2_cdr3_len) + 
  geom_bar(aes(x=cdr3Length, fill = canonical)) + 
  xlim(breaks) +
  ggsave("NKT_beta_canonical_cdr3Length.pdf", width=10, height=7)

library(dplyr)
d3 = df_chain2_cdr3_len %>% group_by(cdr3Length,canonical) %>% summarise(count=n()) %>% as.data.frame(.)
d3$percent = apply(d3,1,function(x){as.numeric(x[3])/table(df_chain2_cdr3_len$canonical)[x[2]]})
write.csv(d3,file="NKT_beta_cdrLen.csv")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
gcolor = gg_color_hue(nlevels(factor(d3$canonical)))
cnt = 0
for(i in levels(factor(d3$canonical))){
  cnt = cnt + 1
  ggplot(subset(d3, canonical == i)) + 
    geom_bar(aes(x = factor(cdr3Length), y = percent*100, fill = canonical),stat="identity", width = 0.7) + 
    scale_fill_manual(values=c(gcolor[cnt])) +
    xlim(breaks) + 
    ggsave(paste0("NKT_beta_cdr3Length",i,".pdf"), width=10, height=7)
}

d2 <- df_chain2_cdr3_len %>% 
  group_by(cdr3Length,canonical) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

ggplot(d2) + 
  geom_bar(aes(x = factor(cdr3Length), y = perc*100, fill = canonical),stat="identity", width = 0.7) + 
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
  ) +
  ggsave("NKT_beta_canonical_cdr3Length_percentage.pdf", width=10, height=7)


write.csv(df_chain2_cdr3_len, file="NKT_beta_canonical_cdr3Length_count.csv")