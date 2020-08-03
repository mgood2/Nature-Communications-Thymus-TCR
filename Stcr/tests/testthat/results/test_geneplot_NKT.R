# init
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.NKT.RData") # seurat v2 file
load("../testdata/ident_NKT.RData")
load('../testdata/ensemblGenes2018-10-15.RData')

gene = c("Hes1","Bex1","Taf9","Hnrnpab","Snrpc","Pcbp2","Zfp644","Smarcb1","Aes")

S.genes = c("Dnajc2",
            "Mcm2",
            "Mcm3",
            "Mcm4",
            "Mki67",
            "Mre11a",
            "Msh2",
            "Pcna",
            "Rad17",
            "Rad51",
            "Sumo1")

N1_gene = c("Pclaf",
            "Birc5",
            "Mki67",
            "Tk1",
            "Cdk1",
            "Ccna2",
            "Spc24",
            "Pbk",
            "Clspn",
            "Ncapg",
            "Gm15428",
            "Cenpm",
            "Cdca3",
            "Tpx2",
            "Knl1",
            "Shcbp1",
            "Rad51ap1",
            "Aurkb",
            "Ncaph",
            "Cep55",
            "Esco2",
            "Racgap1",
            "Bub1b",
            "Sgo1",
            "Rrm2",
            "Ska1",
            "Spc25",
            "Hmmr",
            "Ndc80",
            "Ckap2l",
            "Ccnb2",
            "Fbxo5",
            "Plk1",
            "Rad51",
            "Cenpe",
            "Neil3",
            "Cenpf",
            "Nuf2",
            "Mxd3",
            "Cit",
            "Gm11223",
            "Ube2c",
            "Uhrf1",
            "Kif11",
            "Ccne1",
            "Nusap1",
            "Slc43a3",
            "Cdca8",
            "Prc1",
            "Asf1b",
            "Ttk",
            "E2f8",
            "Cdc6",
            "Tacc3",
            "Kif15",
            "Cdc45",
            "Mybl2",
            "Top2a",
            "Kif22",
            "Ncapg2",
            "Kif4",
            "Kif20a",
            "Hist1h2ak",
            "Cdca5",
            "Ccnb1",
            "Sgo2a",
            "Mad2l1",
            "Cenpw",
            "Cdca2",
            "Ccnf",
            "Lockd",
            "Trip13",
            "Tyms",
            "Chaf1a",
            "Hist1h1b",
            "Diaph3",
            "Lig1",
            "Aurka",
            "Cdc20",
            "Kif23",
            "Knstrn",
            "Dhfr",
            "Ncapd2",
            "Pole",
            "Stmn1",
            "Mis18bp1",
            "Ect2",
            "Cenps",
            "Arhgap11a",
            "Nrm",
            "E2f2",
            "Zfp367",
            "Fignl1",
            "Gm4739",
            "Gins2",
            "Hist1h2ac",
            "Cenph",
            "Rrm1",
            "Ccdc34",
            "Plk4",
            "Dscc1",
            "Lmnb1",
            "Syce2",
            "Ezh2",
            "Pkmyt1",
            "Tcf19",
            "Hist1h2ap",
            "Cks1b",
            "Prim1",
            "Chtf18",
            "Hist1h2ae",
            "Nsd2",
            "H2afz",
            "Gm8203",
            "Ptma",
            "Gmnn",
            "Cdkn2c",
            "Hmgb2",
            "Fen1",
            "Kif20b",
            "Atad5",
            "Mcm5",
            "Ube2t",
            "Rfc4",
            "Smc2",
            "Hells",
            "E2f1",
            "Mcm3",
            "Brca2",
            "Gm6594",
            "Hmgn2",
            "Rfc5",
            "Gm7125",
            "Trim59",
            "Pttg1",
            "Psmc3ip",
            "Pola1",
            "Gm3756",
            "Hmgb3",
            "Gm15452",
            "Tipin",
            "Fam111a",
            "Cks2",
            "Gm10282",
            "Chaf1b",
            "Gm4316",
            "Tuba1b",
            "Tmpo",
            "Ran",
            "Hist1h1a",
            "Dut",
            "Dtl",
            "Haus4",
            "Mcm7",
            "Suv39h1",
            "Rmi2",
            "Ppia",
            "Timeless",
            "Prim2",
            "Prdx4",
            "Cenpl",
            "Orc6",
            "Dlgap5",
            "Mcm6",
            "Smc4",
            "Tubb5",
            "Pmf1",
            "Hmgb1",
            "Pole2",
            "Selenoh",
            "Hnrnpa3",
            "Rfc3",
            "Hist1h1d",
            "Hirip3",
            "Ska2",
            "Hnrnpab",
            "Incenp",
            "Anp32b",
            "Cbx5",
            "Hint1",
            "Slbp",
            "Tubg1",
            "Nasp",
            "Ranbp1",
            "Pola2",
            "Pcna",
            "Mms22l",
            "Hist1h3e",
            "Rpa2",
            "Dck",
            "Marcksl1",
            "Csrp1",
            "Gm4617",
            "Ung",
            "Dnajc9",
            "Srsf3",
            "Pold1",
            "Anp32e",
            "Gm6166",
            "Atad2",
            "Lsm2",
            "Prdx1",
            "Cdca7",
            "Cdkn2a",
            "Dek",
            "Nup85",
            "Dctpp1",
            "Hat1",
            "H2afv",
            "Tfdp1",
            "Carhsp1",
            "Fabp5",
            "Rcc1",
            "Siva1",
            "Mcm4",
            "Snrpd1",
            "Hsp90aa1",
            "Rfwd3",
            "Txn1",
            "Mcm2",
            "Nap1l1",
            "Cdt1",
            "Dbf4",
            "Ube2e3",
            "H2afx",
            "Gapdh",
            "Srsf2",
            "Ybx1",
            "Nucks1",
            "Kpna2",
            "Cmc2",
            "Nelfe",
            "Atp5g3",
            "Alyref",
            "Hmgn1",
            "Banf1",
            "Psat1",
            "Gm20628",
            "Pdlim1",
            "Dnmt1",
            "Ube2s",
            "Casp3",
            "Cdk2",
            "Nrgn",
            "Rfc2",
            "Snrpb",
            "Gm10184",
            "Lsm3",
            "Tceal9",
            "Sumo2",
            "Lsm5",
            "Mrpl18",
            "Cdc25b",
            "Hist1h4k",
            "Snrpe",
            "Ldha",
            "Hnrnpa2b1",
            "Smtn",
            "Usp1",
            "Cycs",
            "Raly",
            "Eri1",
            "U2af1",
            "Ssrp1",
            "Gm4204",
            "Ncapd3",
            "1500009L16Rik",
            "Tpi1",
            "Rbm3",
            "Pold2",
            "Srsf7",
            "Pkm",
            "Impa2",
            "Hist2h2ac",
            "Nutf2",
            "Erh",
            "AI662270",
            "Dtymk",
            "Rbbp4",
            "Anapc5",
            "Topbp1",
            "Tubb4b",
            "Snrpf",
            "Xpo1",
            "Hist1h4d",
            "Ppil1",
            "Slc25a5",
            "Hpf1",
            "Rpa3",
            "Cenpa",
            "Naa40",
            "Nudc",
            "Lbr",
            "Emp3",
            "Rbbp7",
            "Igfbp4",
            "Smc1a",
            "Cdk4",
            "Hnrnpd",
            "Sae1",
            "Hmga1",
            "Atpif1",
            "Plp2",
            "Mthfd2",
            "Rpa1",
            "Pa2g4",
            "Ddx39",
            "Nudt5",
            "Dpy30",
            "Slc29a1",
            "Eif5a",
            "Uchl5",
            "Il13",
            "Rps27l",
            "Cox7a2",
            "Hist1h1e",
            "Ncaph2",
            "Nup62",
            "Mif",
            "Snrpert",
            "Dbi",
            "Supt16",
            "Hprt",
            "Ifi30",
            "Rangap1",
            "Tuba1c",
            "Ilf2",
            "Rad21",
            "Hist2h2bb",
            "Ewsr1",
            "Phgdh",
            "Lyar",
            "Set",
            "Bex3",
            "Lsm8",
            "Eif4a1",
            "Hjurp",
            "Cbx1",
            "Prelid3b",
            "Dynll2",
            "Tagln2",
            "Ywhah",
            "Anapc11",
            "Ptges3",
            "Higd1a",
            "Fus",
            "Klf2",
            "Tkt",
            "Actn4",
            "Nsmce4a",
            "Psip1",
            "Rfc1",
            "Anxa2",
            "Ncl",
            "Cmtm7",
            "Nop58",
            "Paics",
            "Lgals1",
            "Ly6e",
            "Calm2",
            "Cdkn2d",
            "Vim")

ensembl_g1 = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% gene]
ensembl_N1_gene = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N1_gene]
ensembl_s_gene = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% S.genes]


for(i in ensembl_s_gene){
df = data.frame(x=seurat.NKT@dr$umap@cell.embeddings[,1], 
                y=seurat.NKT@dr$umap@cell.embeddings[,2], 
                expression=seurat.NKT@scale.data[i,])
df = df[order(df$expression),]
ggplot(df,aes(x=x, y=y, colour=expression)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$expression)),
                         labels=c(0,round(as.numeric(max(df$expression)), digits = 2))) +
  ylab("UMAP2") +
  xlab("UMAP1") +
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
  ) + ggsave(paste0("GenePlot_NKT_",ensemblGenes[i,"external_gene_name"],".pdf"), width = 7, height=7)
}

ensembl_g1 = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% gene]
ensembl_N1_gene = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N1_gene]

score_g1 = colMeans(as.matrix(seurat.NKT@scale.data[rownames(seurat.NKT@scale.data) %in% ensembl_g1,]))
df = data.frame(x=seurat.NKT@dr$umap@cell.embeddings[,1], 
                y=seurat.NKT@dr$umap@cell.embeddings[,2], 
                expression=score_g1)
df = df[order(df$expression),]
ggplot(df,aes(x=x, y=y, colour=expression)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$expression)),
                         labels=c(0,round(as.numeric(max(df$expression)), digits = 2))) +
  ylab("UMAP2") +
  xlab("UMAP1") +
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
  ) + ggsave(paste0("GenePlot_NKT_bulkgene_marker.pdf"), width = 7, height=7)

score_N1_gene  = colMeans(as.matrix(seurat.NKT@scale.data[rownames(seurat.NKT@scale.data) %in% ensembl_N1_gene,]))
df = data.frame(x=seurat.NKT@dr$umap@cell.embeddings[,1], 
                y=seurat.NKT@dr$umap@cell.embeddings[,2], 
                expression=score_N1_gene)
df = df[order(df$expression),]
ggplot(df,aes(x=x, y=y, colour=expression)) +
  geom_point() + 
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(df$expression)),
                         labels=c(0,round(as.numeric(max(df$expression)), digits = 2))) +
  ylab("UMAP2") +
  xlab("UMAP1") +
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
  ) + ggsave(paste0("GenePlot_NKT_sc_new_marker.pdf"), width = 7, height=7)




load(file="../testdata/NKTp_gdTP.RData")
colnames(NKTp_gdTP)[7] = "Pval"
NKTp_gdTP_FC2_p0.001 = subset(NKTp_gdTP,log2FC > log2(2) & Pval < 0.001)[,"GeneID"]
score_NKTp  = colMeans(as.matrix(seurat.NKT@scale.data[rownames(seurat.NKT@scale.data) %in% NKTp_gdTP_FC2_p0.001,]))

cor(score_N1_NKTp,score_N1_gene, method="spearman")
table(ensembl_N1_gene %in% NKTp_gdTP_FC2_p0.001)
table(NKTp_gdTP_FC2_p0.001 %in% ensembl_N1_gene)
plot(scale(score_N1_NKTp),scale(score_N1_gene))
lines(x=c(-3,0,5),y=c(-3,0,5), col="red")
dev.off()




NKT_box = data.frame(cellnames=c(colnames(seurat.NKT@data),colnames(seurat.NKT@data)) )

scale_score_NKTp = scale(score_NKTp)
scale_score_N1 = scale(score_N1_gene)


NKT_box$exprs = c(scale_score_NKTp,scale_score_N1)

NKT_box$marker = c(rep("NKTp", length(score_NKTp)), rep("N1", length(score_N1_gene)))



library(grid)
library(gridExtra)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

N1_box = NKT_box[NKT_box$cellnames %in% names(seurat.NKT@ident)[seurat.NKT@ident=="N1"],]
N1_plt = ggplot(N1_box) + geom_boxplot(aes(x=marker, y=exprs),fill=gg_color_hue(2)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(range(NKT_box$exprs)) +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=30),
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

N2_box = NKT_box[NKT_box$cellnames %in% names(seurat.NKT@ident)[seurat.NKT@ident=="N2"],]
N2_plt = ggplot(N2_box) + geom_boxplot(aes(x=marker, y=exprs),fill=gg_color_hue(2)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(range(NKT_box$exprs)) +
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

N3_box = NKT_box[NKT_box$cellnames %in% names(seurat.NKT@ident)[seurat.NKT@ident=="N3"],]
N3_plt = ggplot(N3_box) + geom_boxplot(aes(x=marker, y=exprs),fill=gg_color_hue(2)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(range(NKT_box$exprs)) +
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
N4_box = NKT_box[NKT_box$cellnames %in% names(seurat.NKT@ident)[seurat.NKT@ident=="N4"],]
N4_plt = ggplot(N4_box) + geom_boxplot(aes(x=marker, y=exprs),fill=gg_color_hue(2)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(range(NKT_box$exprs)) +
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
N5_box = NKT_box[NKT_box$cellnames %in% names(seurat.NKT@ident)[seurat.NKT@ident=="N5"],]
N5_plt = ggplot(N5_box) + geom_boxplot(aes(x=marker, y=exprs),fill=gg_color_hue(2)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(range(NKT_box$exprs)) +
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
N6_box = NKT_box[NKT_box$cellnames %in% names(seurat.NKT@ident)[seurat.NKT@ident=="N6"],]
N6_plt = ggplot(N6_box) + geom_boxplot(aes(x=marker, y=exprs),fill=gg_color_hue(2)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(range(NKT_box$exprs)) +
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
N7_box = NKT_box[NKT_box$cellnames %in% names(seurat.NKT@ident)[seurat.NKT@ident=="N7"],]
N7_plt = ggplot(N7_box) + geom_boxplot(aes(x=marker, y=exprs),fill=gg_color_hue(2)) +
  # xlab(paste0("M1\n Pvalue = ",format(M1_t$p.value,digits = 2))) +
  ylim(range(NKT_box$exprs)) +
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




pdf("NKT_NKTp_N1_marker_2color.pdf",width=12,height =7)
grid.arrange(N1_plt,N2_plt,N3_plt,N4_plt,N5_plt,N6_plt,N7_plt, ncol=7)
dev.off()











