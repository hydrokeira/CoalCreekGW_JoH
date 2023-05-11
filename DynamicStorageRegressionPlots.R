#dynamic storage plots
require(ggplot2)
require(ggpmisc)
require(ggpubr)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")
#read in dynamic storage csv
DS<-read.csv("DynamicStoragePlots.csv")
colnames(DS)<-c("Date", "SOI", "MedGWFlux", "MeanGWFlux", "DS", "PrecipSum_7days", "Q", "Qavg_7days")

GW<-read.csv("ModelOutputGWSummary.csv")

GW<-GW[,c(3:6, 10,13:15)]

DS<-merge(DS, GW, by="Date")

#make plots

p1<-ggplot(DS, aes(DS, SOI))+geom_point(size=3)+geom_smooth(method = "lm", se=FALSE, col="black")+theme_classic()+
  stat_correlation(mapping = use_label(c("R2", "P")), parse=TRUE, size=7, label.x="right")+
  theme(text = element_text(size=25))+labs(x="Dynamic Storage (mm)", y="Seasonal Origin Index")
  

p2<-ggplot(DS, aes(DS, MeanGWFlux))+geom_point(size=3)+geom_smooth(method = "lm", se=FALSE, col="black")+theme_classic()+
  stat_correlation(mapping = use_label(c("R2", "P")), parse=TRUE, size=7, label.x="left")+
  theme(text = element_text(size=25))+labs(x="Dynamic Storage (mm)", y="Mean GW Flux (m/s)")

p3<-ggplot(DS, aes(DS, Q))+geom_point(size=3)+geom_smooth(method = "lm", se=FALSE, col="black")+theme_classic()+
  stat_correlation(mapping = use_label(c("R2", "P")), parse=TRUE, size=7, label.x="left")+
  theme(text = element_text(size=25))+labs(x="Dynamic Storage (mm)", y="Discharge (mm)")

p4<-ggplot(DS, aes(DS, fracture_alluvial))+geom_point(size=3)+geom_smooth(method = "lm", se=FALSE, col="black")+theme_classic()+
  stat_correlation(mapping = use_label(c("R2", "P")), parse=TRUE, size=7, label.x="left")+
  theme(text = element_text(size=25))+labs(x="Dynamic Storage (mm)", y="Fracture:Alluvial GW Ratio")


pdf("DS_Plots.pdf", family = "Times", width = 15, height = 10)

ggarrange(p1, NULL, p3, p2, NULL, p4, widths = c(1, 0.05, 1, 1, 0.05, 1))+theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "in"))

dev.off()

#test sig
lm1<-lm(DS~fracture_alluvial, DS)
summary(lm1)

max(DS$DS)-min(DS$DS)


