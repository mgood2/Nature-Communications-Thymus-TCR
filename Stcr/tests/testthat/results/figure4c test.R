library(xlsx)

dfNKT = read.xlsx2("D:/Dropbox/paper/Manuscript/Ncomm submission/Table S3. List of clonotypes.xlsx", sheetName = "NKT")
dfMAIT = read.xlsx2("D:/Dropbox/paper/Manuscript/Ncomm submission/Table S3. List of clonotypes.xlsx", sheetName = "MAIT")
dfgdT = read.xlsx2("D:/Dropbox/paper/Manuscript/Ncomm submission/Table S3. List of clonotypes.xlsx", sheetName = "gamma delta T")


sort(names(table(df$n)), decreasing = T)
table(dfNKT$n)[order(as.numeric(names(table(dfNKT$n))))]

nkt_test = table(dfNKT$n)[order(as.numeric(names(table(dfNKT$n))))]
mait_test = table(dfMAIT$n)[order(as.numeric(names(table(dfMAIT$n))))]
gdt_test = table(dfgdT$n)[order(as.numeric(names(table(dfgdT$n))))]


ggdf_all = rbind(as.data.frame(nkt_test), as.data.frame(mait_test),as.data.frame(gdt_test))
ggdf_all$Var1 = as.numeric(ggdf_all$Var1)
ggdf_all$type = c(rep("NKT",length(table(dfNKT$n))),rep("MAIT",length(table(dfMAIT$n))),rep("gdT",length(table(dfgdT$n))))
ggdf_all$type = factor(ggdf_all$type, levels = c("NKT","MAIT","gdT"))
library(ggplot2)

ggdf_all

ggdf_all$factors = c(rep(1,nrow(ggdf)),rep(3285/2287,nrow(ggdf1)),rep(3285/2667,nrow(ggdf2)))
ggdf_all$Norm = ggdf_all$Var1 * ggdf_all$factors
ncol(seurat.NKT@data) #3285 ==> 1
ncol(seurat.MAIT@data) #3285/2287
ncol(seurat.gdT@data) #3285/2667

ggplot(ggdf_all) + geom_point(aes(x=Var1,y=-log10(Freq),colour=type)) + geom_line(aes(x=Var1,y=-log10(Freq),colour=type)) +
  scale_x_continuous(breaks=c(0, 5, 10, 15)) #+ xlim(0,20)
