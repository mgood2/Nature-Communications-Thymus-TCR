ComputeShannonIndex <- function(vector){
  # k = length(vector)
  k = sum(vector > 0)
  return(-sum((vector/sum(vector))*log(vector/sum(vector)), na.rm=TRUE)/log(k))
}


# nkT
library(Stcr)
setwd("D:/elee4472/OneDrive - dgist.ac.kr/1.Projects/R projects/Stcr/tests/testthat")
load("../testdata/seurat.NKT.RData") # seurat v2 file
load("../testdata/ident_NKT.RData")
load("../testdata/ensemblGenes2018-10-15.RData")

ssetF3 = Load10xVDJ('../testdata/mTCR-S1/outs_IMGT/', filtered = FALSE)
ssetF5 = Load10xVDJ('../testdata/mTCR-S2/outs_IMGT/', filtered = FALSE)

ssetF3 = SetSample(ssetF3,"F3")
ssetF5 = SetSample(ssetF5,"F5")


sset = MergeSamples(ssetF3,ssetF5) 

sset_NKT = ImportSeurat(sset,seurat.NKT, version="v2")

NKT_cdr3_mat = SetCDR3Mat(sset_NKT, cell_ind = colnames(seurat.NKT@scale.data),
                          chain1 = c("TRA"), chain2 = c("TRB"))

NKT_cdr3_pair_mat = SetCDR3PairMat(sset_NKT, NKT_cdr3_mat, chain1 = c("TRA"), chain2 = c("TRB"))


mat_cdr3_pair = NKT_cdr3_pair_mat[[1]]
df_clonotype = NKT_cdr3_pair_mat[[2]]
ident = as.factor(seurat.NKT@ident)
df_mat_cdr3_pair = data.frame(mat_cdr3_pair)

ggdf = data.frame(table(colSums(df_mat_cdr3_pair)))
NKTNN = table(colSums(df_mat_cdr3_pair))

ks.test(NKTNN,gdtNN)
ks.test(gdtNN,MAITNN)
ks.test(NKTNN,MAITNN)
p.adjust(c(0.7499,0.1294,0.0004329),method = "bonferroni", n= 3)


#ks test
NKTNNall = c()
NKTNNm1 = NKTNN[2:length(NKTNN)]

for(i in names(NKTNNm1)){
 NKTNNall = c(NKTNNall, rep(i, as.numeric(NKTNNm1[i])))
}
NKTNNall = as.numeric(NKTNNall)

gdtNNall = c()
gdtNNm1 = gdtNN[2:length(gdtNN)]

for(i in names(gdtNNm1)){
  gdtNNall = c(gdtNNall, rep(i, as.numeric(gdtNNm1[i])))
}
gdtNNall = as.numeric(gdtNNall)

MAITNNall = c()
MAITNNm1 = MAITNN[2:length(MAITNN)]

for(i in names(MAITNNm1)){
  MAITNNall = c(MAITNNall, rep(i, as.numeric(MAITNNm1[i])))
}
MAITNNall = as.numeric(MAITNNall)


ks.test(NKTNNall,gdtNNall)
ks.test(gdtNNall,MAITNNall)
ks.test(NKTNNall,MAITNNall)
p.adjust(c(0.7499,0.1294,0.0004329),method = "bonferroni", n= 3)




plot(y=ecdf(gdtNNall)(c(unique(gdtNNall))),x=c(unique(gdtNNall)), col="blue",log="x")
points(y=ecdf(NKTNNall)(c(unique(NKTNNall))),x=c(unique(NKTNNall)), col="red",log="x")
points(y=ecdf(MAITNNall)(c(unique(MAITNNall))),x=c(unique(MAITNNall)), col="green",log="x")

ggdf1 = data.frame(MAITNN)
ggdf2 = data.frame(gdtNN)
# NKTNN.backup = NKTNN
# MAITNN.backup = MAITNN
# gdtNN.backup = gdtNN

# NKTNN = NKTNN[order(NKTNN, decreasing = T)]
# MAITNN = MAITNN[order(MAITNN, decreasing = T)]
# gdtNN = gdtNN[order(gdtNN, decreasing = T)]


nkt_test = table(dfNKT$n)[order(as.numeric(names(table(dfNKT$n))))]
mait_test = table(dfMAIT$n)[order(as.numeric(names(table(dfMAIT$n))))]
gdt_test = table(dfgdT$n)[order(as.numeric(names(table(dfgdT$n))))]

ggdf_all = rbind(ggdf, ggdf1,ggdf2)
ggdf_all$Var1 = as.numeric(ggdf_all$Var1)
ggdf_all$type = c(rep("NKT",nrow(ggdf)),rep("MAIT",nrow(ggdf1)),rep("gdT",nrow(ggdf2)))
ggdf_all$type = factor(ggdf_all$type, levels = c("NKT","MAIT","gdT"))
library(ggplot2)

ggdf_all

ggdf_all$factors = c(rep(1,nrow(ggdf)),rep(3285/2287,nrow(ggdf1)),rep(3285/2667,nrow(ggdf2)))
ggdf_all$Norm = ggdf_all$Freq * ggdf_all$factors
ncol(seurat.NKT@data) #3285 ==> 1
ncol(seurat.MAIT@data) #3285/2287
ncol(seurat.gdT@data) #3285/2667

library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}



ggplot(ggdf_all) + geom_point(aes(x=Var1,y=Norm,colour=type),shape=1, size=6) + geom_line(aes(x=Var1,y=Norm,colour=type),size=1.5) +
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20))+
  scale_y_continuous(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    trans=reverselog_trans(10)
    
  )  +
  
  #scale_y_continuous(labels = trans_format("log10", math_format(10^.x)),trans=reverselog_trans(10)) +
  xlim(0,30) +
  xlab("Number of cells") +
  ylab("Number of clonotype\n(Normalized)") +
  # ggplot2::annotation_logticks(side="l") +
  scale_colour_manual(values=c("red","green","blue")) +
  theme_bw() +
  theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    plot.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), 
    axis.line = element_line(),
    axis.title=element_text(size=45),
    axis.text=element_text(size=30),
    legend.position = "none") + ggsave("figure4c.pdf",width=7,height=7)

p +   
ggplot(ggdf_all) + geom_point(aes(x=Var1,y=-log10(Freq),colour=type)) + geom_line(aes(x=Var1,y=-log10(Freq),colour=type)) +
  scale_x_continuous(breaks=c(50, 100, 150)) + xlim(50,100)

p = ggplot(ggdf, aes(x=Norm,y=-log10(Freq))) + geom_point(  color="red") + geom_line(color="red", group="1")
p = p + geom_point( ggdf1, aes(x=Var1,y=-log10(Freq)), color="green")
ggplot(ggdf2, aes(x=Var1,y=-log10(Freq))) + geom_point(  color="blue") + geom_line(color="blue", group="1")

