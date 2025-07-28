#############################################################################
#################Zebra Finch Neural Repro Mechanisms#########################
##############################################################################
#Created by Jonathan H Perez 
#Based on work by Jessica Moodie (31-8-2023)
#This file started 12-5-2025


#Purpose:Graphical rework of figures from Jessica Moodie's
#Master's thesis to make composite figures for publiction purposes
#Revisiting and refinement of statisics, confirmation for JHP
#of statistics prior to submission

#####Library Load#####
library(tidyverse)
library(plotrix)
library(lmerTest)
library(car)
library(performance)
library(emmeans)
library(gridExtra)
library(cowplot)

#####Data Loading and Working Directory######
dir()
setwd("C:/Users/jhperez/Documents/Research/Zebra Finch Neural Mechanisms of Aseasonal Breeding")

####Load data for updated plots####
#same data from analyses just composited for easier plotting

d1 <- read.csv("ZebraFinch_WaterDep_MasterData-qPCR_Morph.csv")
str(d1)
d1.m <- filter(d1, Sex=="M")
d1.f <- filter(d1, Sex=="F")

#planned color choices #CC5500 for burnt orange #D3D3D3 for light gray

#####Graphing#####
#testes mass
g.testesm <- ggplot(d1.m,aes(x=Treatment, y=Testes.Mass..g.))+
  geom_boxplot(fill="#D3D3D3")+
  theme_bw()+
  labs(y="Testes Mass (g)")+
  theme(axis.title = element_text(face=2,size=16), 
        axis.text = element_text(face=2, size=14))
g.testesm

#ovary mass
g.ovarym <- ggplot(d1.f,aes(x=Treatment, y=Ovary.Mass..mg.))+
  geom_boxplot(fill="#CC5500")+
  theme_bw()+
  labs(y="Ovary Mass (mg)")+
  theme(axis.title = element_text(face=2,size=16), 
        axis.text = element_text(face=2, size=14))
g.ovarym

#Largest follicle size
g.foll <- ggplot(d1.f,aes(x=Treatment, y=Largest.Follicle.size..mm.))+
  geom_boxplot(fill="#CC5500")+
  geom_jitter(width=0.2, height=0)+
  #geom_point()+
  theme_bw()+
  labs(y="Largest Follicle Diameter (mm)")+
  theme(axis.title = element_text(face=2,size=16), 
        axis.text = element_text(face=2, size=14))
g.foll

Anova(lm(Largest.Follicle.size..mm.~Treatment, d1.f), type="III")

#GnRH
g.gnrh <- ggplot(d1, aes(x=Treatment, y=GnRH, fill=Sex, col=Sex))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0)+
  #Need to figure out the color and spacing on the jitters
  theme_bw()+
  scale_fill_manual(values=c("#CC5500","#D3D3D3"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="GnRH mRNA Expression")+
  theme(axis.title = element_text(face=2,size=16), 
        axis.text = element_text(face=2, size=14))+
  theme(legend.position = "None")
g.gnrh

#GnIH
g.gnih <- ggplot(d1, aes(x=Treatment, y=GnIH, fill=Sex, col=Sex))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0)+
  theme_bw()+
  scale_fill_manual(values=c("#CC5500","#D3D3D3"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="GnIH mRNA Expression")+
  theme(axis.title = element_text(face=2,size=16), 
        axis.text = element_text(face=2, size=14))+
  ylim(0,1500)+
  theme(legend.position = "None")
g.gnih
#excluding outlier in plotting and likely in analysis as well

#FSH
g.fsh <- ggplot(d1, aes(x=Treatment, y=FSH, fill=Sex, col=Sex))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0)+
  theme_bw()+
  scale_fill_manual(values=c("#CC5500","#D3D3D3"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="FSH mRNA Expression")+
  theme(axis.title = element_text(face=2,size=16), 
        axis.text = element_text(face=2, size=14))+
  theme(legend.position = c(0.95, 0.95),  # (x, y) coordinates in [0,1]
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black")
  )
g.fsh

#DIO2
g.dio2 <- ggplot(d1, aes(x=Treatment, y=DIO2, fill=Sex, col=Sex))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0)+
  theme_bw()+
  scale_fill_manual(values=c("#CC5500","#D3D3D3"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="DIO2 mRNA Expression")+
  theme(axis.title = element_text(face=2,size=16), 
        axis.text = element_text(face=2, size=14))+
  theme(legend.position = "None")
g.dio2

#DIO3

#DIO2
g.dio3 <- ggplot(d1, aes(x=Treatment, y=DIO3, fill=Sex, col=Sex))+
  geom_boxplot()+
  geom_jitter(width=0.2, height=0)+
  theme_bw()+
  scale_fill_manual(values=c("#CC5500","#D3D3D3"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="DIO3 mRNA Expression")+
  theme(axis.title = element_text(face=2,size=16), 
        axis.text = element_text(face=2, size=14))+
  theme(
    legend.position = c(0.15, 0.95),  # (x, y) coordinates in [0,1]
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black")
  )
g.dio3

#Combining Graphs
##Figure 1
#figure 1 is a schematic of the experiment and the timeline. Might add egg count
#tables to that as well.

###Figure 2
Fig.2 <- cowplot::plot_grid(g.dio2, g.dio3, g.gnrh,g.gnih, ncol=2, align="hv")
ggsave("C:/Users/jhperez/Documents/Manuscripts/Zebra_Finch_Water_Deprivation/Figures/Figure_2.pdf"
       , dpi=300, width =14, height=10)
ggsave("C:/Users/jhperez/Documents/Manuscripts/Zebra_Finch_Water_Deprivation/Figures/Figure_2.png"
       , dpi=300, width =14, height=10)


##Figure 3

Fig.3 <- cowplot::plot_grid(g.fsh, g.foll, g.testesm, g.ovarym, ncol=2, align="hv")
Fig.3
ggsave("C:/Users/jhperez/Documents/Manuscripts/Zebra_Finch_Water_Deprivation/Figures/Figure_3.pdf"
       , dpi=300, width =14, height=10)
ggsave("C:/Users/jhperez/Documents/Manuscripts/Zebra_Finch_Water_Deprivation/Figures/Figure_3.png"
       , dpi=300, width =14, height=10)

########Follow up or exploratory stats########
#For my sanity doing a re-run of all stas performed by Moodie J here.

#GnrH
m.gnrh<-lm(GnRH~Treatment*Sex, data=d1)
Anova(m.gnrh)
check_model(m.gnrh)
check_outliers(m.gnrh) #case 20, likely a second
hist(d1$GnRH)

d2[20, 6] <- NA
#removing outlier
m.gnrh.1<-lm((GnRH*100)~Treatment*Sex, data=d2)
Anova(m.gnrh.1)
check_model(m.gnrh.1)
check_outliers(m.gnrh.1) #case 10 now added, 
#switching to gamma

m.gnrh.3 <- glm(GnRH~Treatment*Sex, data=d2, family = Gamma(link = "inverse"))
summary(m.gnrh.3) #won't converge with log link forced with inverse
check_model(m.gnrh.3)
#none of the gnrh models are particularly pretty looking but they all say no
#difference between estimators

#GnIH
m.gnih<-lm(GnIH~Treatment*Sex, data=d1)
Anova(m.gnih)
check_model(m.gnih)
check_outliers(m.gnih) #case 4
d2[4,7] <- NA

m.gnih<-lm(log(GnIH)~Treatment*Sex, data=d2)
Anova(m.gnih)
check_model(m.gnih)
check_outliers(m.gnih)

#FSH
m.fsh <- lm(FSH~Treatment*Sex, data=d1)
Anova(m.fsh, type="III")
check_model(m.fsh)
#normality of residuals is a mess
boxcox(m.fsh)#suggestes squareroot transformation

m.fsh2 <- lm(sqrt(FSH)~Treatment*Sex,data=d1)
Anova(m.fsh2, type="III")
check_model(m.fsh2)
#well that made model fit worse for normality interestingly
#same inference

m.fsh3 <- lm((FSH^0.75)~Treatment*Sex,data=d1)
Anova(m.fsh3, type="III")
check_model(m.fsh3)
check_outliers(m.fsh3)
#still floats around the normality line no matter what we try

#inference of treatment as the main effect remains the same across all models
#is also consistent with the graphical presentation of the data

#Trying gamma model
m.fsh.4 <- glm(FSH~Treatment*Sex, family=Gamma(link = "log"), data=d1)
summary(m.fsh.4)
check_model(m.fsh.4)#best looking fit of the models

################################################################################
#DIO2 reanalysis
m.dio2 <- lm(DIO2~Treatment*Sex, data=d1)
Anova(m.dio2, type="III")
check_model(m.dio2)
check_outliers(m.dio2)
#outlier detected case 17 equates to row 20 in data due to NAs
#outlier removed and data reanalyzed
#with outlier #1 removed we have a significant effect, but the testing flags

d2 <- d1
d2[20,4] <- NA #removing outlier from conrol group
m.dio2.2 <- lm(DIO2~Treatment*Sex, data=d2)
Anova(m.dio2.2, type="III")
check_model(m.dio2.2)
#single point that is off line in residuals normality plot
#suggests log transform
boxcox(m.dio2.2) #95% includes 0 so trying log

m.dio2.3 <- lm(log(DIO2)~Treatment*Sex, data=d2)
Anova(m.dio2.3, type="III") #significance disappears here
check_model(m.dio2.3)
check_outliers(m.dio2.3)
#log transformation addressed the outlier problem 
#but graphicl goodness of fit tests are now a hot hot mess
#think an alternative model may be appropriate to consider

hist(d1$DIO2) #Gamma seems appropriate

#Trying gamma glm, trying with non-outlier removed dataset first

m.dio2.4 <- glm(DIO2~Treatment*Sex, family=Gamma(link="log"), data=d1)
summary(m.dio2.4)
check_outliers(m.dio2.4)
#still picks that C 17 female (row 20) up as an outlier

#retrying reduced dataset
m.dio2.5 <- glm(DIO2~Treatment*Sex, family=Gamma(link="log"), data=d2)
summary(m.dio2.5)
check_outliers(m.dio2.5) #passes outlier tests
check_model(m.dio2.5) #not pefect but much better than others.
Anova(m.dio2.5, type = "III")
################################################################################
#DIO3

m.dio3<-lm(DIO3~Treatment*Sex, data=d1)
Anova(m.dio3, type="III")
#no effect detected
check_model(m.dio3)
check_outliers(m.dio3)
#case 5 is an outlider
d1[,5]
d2[5,5] <- NA #setting outlier to NA

m.dio3.2<-lm(DIO3~Treatment*Sex, data=d2)
Anova(m.dio3.2, type="III")
#still no effect
check_outliers(m.dio3.2) #same problem with sequential outliers
check_model(m.dio3.2)
hist(d1$DIO3) #based on this the initial 35k outlier needs to go regardless
#rest should perform better with a gamma
hist(d2$DIO3)
d3 <- d2[-5,]
m.dio3.3 <- glm(DIO3 ~ Treatment * Sex,
                family = Gamma(link = "log"),
                data = d3)
#finite values seem to be the issue will have to stick to non-gamma





###############################################################################