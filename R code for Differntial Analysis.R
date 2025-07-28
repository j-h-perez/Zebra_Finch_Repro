#####Zebra_Finch#####
###Created July 9th 2024###
###Olivia Breisacher; Jonathan Perez###

BiocManager::install(c("edgeR","limma","clusterProfiler"))

#####Library Load######
library(tidyverse)
library(edgeR)
library(clusterProfiler)
library(org.Gg.eg.db)
library(enrichplot)

getwd()
setwd("/Users/ob/Desktop/SURF")
d1<-read.csv("HTseq_ZFinch_Gonads.csv")
d2<-read.csv("graph_data2.csv") #gonad size data from jess for making plots#
str(d1)

####data processing####
#splitting data set into male and females#

df<-d1[,1:16]
str(df)

dm<-d1[,17:31]
str(dm)

#adding a column for the males#
dm$GENE<-d1$GENE

#reordering the columns rather than adding a column#
dm <- dm[c(16, 1:15)]
str(dm)

#starting pre=processing for Limma Voom#
df1<-DGEList(df)

dm1<-DGEList(dm)

#calculating normalization
df1<-calcNormFactors(df1)
df1

dm1<-calcNormFactors(dm1)
dm1

#Filtering low expressed genes#
cutoff <- 1
drop <- which(apply(cpm(df1), 1, max) < cutoff)
df2 <- df1[-drop,] 
dim(df2) # number of genes left
  
drop <- which(apply(cpm(dm1), 1, max) < cutoff)
dm2 <- dm1[-drop,] 
dim(dm2) # number of genes left

#testing new cutoff for males 16 July 2024#
cutoff <- 1
drop <- which(apply(cpm(dm1), 1, max) < cutoff)
dm2 <- dm1[-drop,] 

dim(dm2)


#treatment metadata
group<-c(rep("T",8),rep("C",7))
group  

#mds plots
plotMDS(df2, col = as.numeric(group))
plotMDS(dm2, col = as.numeric(group))

#creating mean-variance plots
mm <- model.matrix(~0 + group)
yf <- voom(df2, mm, plot = T)
ym <- voom(dm2, mm, plot = T)

#fitting linear model in Limma
fit.f <- lmFit(yf, mm)
head(coef(fit.f))

fit.m <- lmFit(ym, mm)
head(coef(fit.m))

contr.f <- makeContrasts(groupC - groupT, levels = colnames(coef(fit.f)))
contr.f
  
contr.m <- makeContrasts(groupC - groupT, levels = colnames(coef(fit.m)))
contr.m

#females
tmp.f <- contrasts.fit(fit.f, contr.f)
tmp.f <- eBayes(tmp.f)
top.table.f <- topTable(tmp.f, sort.by = "P", n = Inf)
head(top.table.f, 20)
View(top.table.f)


#males
tmp.m <- contrasts.fit(fit.m, contr.m)
tmp.m <- eBayes(tmp.m)
top.table.m <- topTable(tmp.m, sort.by = "P", n = Inf)
head(top.table.m, 20)

#####Gene ontogeny analysis#####
BiocManager::install ("org.Gg.eg.db", force = F)#use force true to force install
library(org.Gg.eg.db)


#females#
df.de<-dplyr::filter(top.table.f,abs(logFC) >=1.5 & adj.P.Val <= 0.05)
df.de.id <- bitr (df.de$GENE, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Gg.eg.db)
ego.f<-enrichGO(df.de.id$ENTREZID, OrgDb ="org.Gg.eg.db", ont="ALL", 
                pAdjustMethod = "fdr")
library(enrichplot)
ego.f<- pairwise_termsim(ego.f)

emapplot(ego.f, cex_label_category=.8, cex_line=.5)+coord_cartesian()
#goplot(ego.f) not working troubleshot#




#males#
#did not work#
dm.de<-dplyr::filter(top.table.m,abs(logFC) >=1.5 & P.Value <= 0.05)
BiocManager::install ("org.Gg.eg.db", force = T)
library(org.Gg.eg.db)
dm.de.id <- bitr (dm.de$GENE, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Gg.eg.db)
ego.m<-enrichGO(dm.de.id$ENTREZID, OrgDb ="org.Gg.eg.db", ont="ALL")
library(enrichplot)
ego.m<- pairwise_termsim(ego.m)

emapplot(ego.m, cex_label_category=.8, cex_line=.5)+coord_cartesian()



#16 July 2024#
######creating graphs#####
d2f<-dplyr::filter(d2, Sex=="F" )
d2m<-dplyr::filter(d2, Sex=="M")

#female
g2f<-ggplot(data=d2f, aes (x=Treatment, y=Largest.Follicle.size..mm.,fill=Treatment))+
geom_boxplot()+
theme_bw()+
  scale_fill_manual(values=c("#add3ff","#ffcaa1"))+
theme(axis.title = element_text(size = 16, face = "bold"))+
theme(axis.text = element_text(size=14, face= "bold"))+
labs(y="Largest Follicle Size (mm)")

g2f
ggsave("Fig2.pdf")

#males
g2m<-ggplot(data=d2m, aes (x=Treatment, y=Testes.Mass..g.,fill=Treatment))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c("#add3ff","#ffcaa1"))+
  theme(axis.title = element_text(size = 16, face = "bold"))+
  theme(axis.text = element_text(size=14, face= "bold"))+
  labs(y="Testes Mass (g)")

g2m
ggsave("Fig3.pdf")

#####creating 
star.males <- filter(dm, GENE=="STAR")
#turn into column/pivoting/wide to tall format
star.males<-star.males%>%pivot_longer(
cols = -GENE,
names_to = "Treatment",
values_to = "Value")


star.males$Treatment<-c(rep("T",8),rep("C",7))

#name of graph <- ggplot(data=d2m, aes (x=Treatment, y=CPM)
STAR <-ggplot(data=star.males, aes (x=Treatment, y=Value, fill = Treatment))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c("#add3ff","#ffcaa1"))+
  theme(axis.title = element_text(size = 16, face = "bold"))+
  theme(axis.text = element_text(size=14, face= "bold"))+
  labs(y="CPM of StAR")

 STAR
