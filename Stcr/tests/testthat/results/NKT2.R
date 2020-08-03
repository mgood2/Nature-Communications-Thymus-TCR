# init
library(Seurat)
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")

load("../testdata/ensemblGenes2018-10-15.RData")



##### Code for Tcell Project
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/06/04


library(topGO)
library(xlsx)

N1_geneset = read.xlsx2("../testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                       startRow = 5, colIndex = 2)
test = N1_geneset
test = unique(as.vector(test[,1]))
N1_geneset = test[test != ""]


N2_geneset = read.xlsx2("../testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 3)
test = N2_geneset
test = unique(as.vector(test[,1]))
N2_geneset = test[test != ""]

N3_geneset = read.xlsx2("../testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 4)
test = N3_geneset
test = unique(as.vector(test[,1]))
N3_geneset = test[test != ""]

N4_geneset = read.xlsx2("../testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 5)
test = N4_geneset
test = unique(as.vector(test[,1]))
N4_geneset = test[test != ""]

N5_geneset = read.xlsx2("../testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 6)
test = N5_geneset
test = unique(as.vector(test[,1]))
N5_geneset = test[test != ""]

N6_geneset = read.xlsx2("../testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 7)
test = N6_geneset
test = unique(as.vector(test[,1]))
N6_geneset = test[test != ""]

N7_geneset = read.xlsx2("../testdata/Table S2. List of signature genes.xlsx", sheetIndex = 1,
                        startRow = 5, colIndex = 8)
test = N7_geneset
test = unique(as.vector(test[,1]))
N7_geneset = test[test != ""]

N1_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N1_geneset]
N2_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N2_geneset]
N3_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N3_geneset]
N4_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N4_geneset]
N5_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N5_geneset]
N6_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N6_geneset]
N7_geneset_ensem = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% N7_geneset]




# gdT1_genecluster= read.csv(file="pheatmap_gdT1_genecluster_info.csv", row.names = 1)
# for(i in 1:nrow(gdT1_genecluster)){
#   print(rownames(ensemblGenes)[ensemblGenes$external_gene_name == rownames(gdT1_genecluster)[i]])
#   rownames(gdT1_genecluster)[i] = rownames(ensemblGenes)[ensemblGenes$external_gene_name == rownames(gdT1_genecluster)[i]]
# }
backgroundGenes = rownames(seurat.MAIT@scale.data)
# gdT1_genecluster[gdT1_genecluster$Cluster %in% c(3,7),] = "gdT1_GO_stage1"
# gdT1_genecluster[gdT1_genecluster$Cluster %in% c(6),] = "gdT1_GO_stage2"
# gdT1_genecluster[gdT1_genecluster$Cluster %in% c(5,8,4,1),] = "gdT1_GO_stage3"
# # gdT1_genecluster[gdT1_genecluster$Cluster %in% c(3,4),] = "gdT1_GO_stage4"


NKT2_genes= list(N3_geneset_ensem,
N4_geneset_ensem,
N5_geneset_ensem,
N6_geneset_ensem)
cluster = c("N3", "N4","N5", "N6")
cnt = 1
for(j in NKT2_genes){
  # backgroundGenes = read.csv2(file = "backgroundGenes.txt", header=FALSE)['V1'][,1]
  targetGenes = j
  allGene = factor(as.integer(backgroundGenes %in% targetGenes))
  names(allGene) = backgroundGenes
  library(plyr)
  onts= c("BP")
  tab = as.list(onts)
  names(tab) = onts
  for (i in 1:length(onts)) {
    tgd = new("topGOdata", ontology = onts[i], allGenes=allGene, nodeSize=10,
              annot=annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
    resultTopGO.elim = runTest(tgd, algorithm="elim", statistic = "Fisher")
    tab[[i]] = GenTable(tgd, Fisher.elim = resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes=200)
  }
  
  topGOResults = rbind.fill(tab)
  write.csv(topGOResults, file=paste("topGOResults_NKT2_GO_stage",cluster[cnt],".csv",sep=""))
  cnt = cnt + 1
}
