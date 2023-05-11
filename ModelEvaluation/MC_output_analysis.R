install.packages("dplyr")
install.packages("ggplot2")
require(dplyr)
require(ggplot2)
setwd("C:/Users/johnkeir/RadonModel")
AICcsv_list<-list.files(path = ".", pattern = "AIC")
Radoncsv_list<-list.files(path = ".", pattern = "Radon_")
Gascsv_list<-list.files(path = ".", pattern = "Gas")

AIC_all <- AICcsv_list %>% 
  lapply(readLines) %>% 
  as.data.frame %>%
  bind_rows %>%
  t

Radon_all <- Radoncsv_list %>% 
  lapply(readLines) %>% 
  as.data.frame %>%
  bind_rows %>%
  t

Gas_all <- Gascsv_list %>% 
  lapply(readLines) %>% 
  as.data.frame %>%
  bind_rows %>%
  t

MC_all<-as.data.frame(cbind(AIC_all, Radon_all, Gas_all))
rownames(MC_all)<-paste("Run", seq(1,3000,1))
colnames(MC_all)<-c("AIC", "Radon", "Gas")
MC_all[,c(1:3)]<-lapply(MC_all[,c(1:3)], as.numeric)
MC_all$radon_gas<-MC_all$Radon*MC_all$Gas

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

MC_files<-list.files(path = ".", pattern = "MC_output_all")

MC_vect<-list()

MC_10_vect<-list()

for (i in 1:length(MC_files)) {
  
  MC_all<-read.csv(MC_files[i])
  
  MC_10<-subset(MC_all, MC_all$AIC < sort(MC_all$AIC)[11])
  
  MC_10$date<-MC_files[i]
  
  # m1<-min(MC_10$Radon)
  # m2<-max(MC_10$Radon)
  # m3<-median(MC_10$Radon)
  # 
  # m4<-min(MC_10$Gas)
  # m5<-max(MC_10$Gas)
  # m6<-median(MC_10$Gas)
  # 
  # m7<-mean(MC_10$AIC)
  # 
  # MC_vect[[i]]<-c(m1, m2, m3, m4, m5, m6, m7)
  
  MC_10_vect[[i]]<-MC_10
  
}

# MC_10_df<-do.call(bind_rows, MC_10_vect)
# 
# write.csv(MC_10_df, "MC_10_df.csv")

MC_df<-as.data.frame(do.call(rbind, MC_vect))

MC_df$datafile<-MC_files

colnames(MC_df)<-c("min_rn", "max_rn", "med_rn", "min_gas", "max_gas", "med_gas", "mean_AIC", "datafile")

write.csv(MC_df, "MC_min_max_values.csv")

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek")

##gas exchange
GE<-read.csv("GasExchangeValues.csv")

ggplot(GE, aes(Depth.at.Coal.15, Gas.Median))+theme_classic()+
  geom_point(size=3)+
  xlab("Depth at Coal-15 (m)")+ylab("Gas Exchange Coefficient")+
  geom_smooth(method = "lm", se=FALSE, col="black", lwd=1)+
  theme(text = element_text(size = 30))

lmQ<-lm(Gas.Median~Q.at.Coal.15, GE)
summary(lmQ)

lmD<-lm(Gas.Median~Depth.at.Coal.15, GE)
summary(lmD)

pdf("MC_run_params_plots.pdf")

ggplot(MC_all)+geom_point(aes(Radon, AIC, col=Gas))+theme_classic()

ggplot(MC_all)+geom_point(aes(Gas, AIC, col=Radon))+theme_classic()

ggplot(MC_10)+geom_point(aes(Radon, AIC, col=Gas))+theme_classic()

ggplot(MC_10)+geom_point(aes(Gas, AIC, col=Radon))+theme_classic()

dev.off()

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

MC_all<-read.csv("MC_output_all_0623.csv")

AIC_70<-subset(MC_all, MC_all$AIC < sort(MC_all$AIC)[70])

mean(AIC_70$Radon)
mean(AIC_70$Gas)

median(AIC_70$Radon)
median(AIC_70$Gas)

min(AIC_70$Radon)
min(AIC_70$Gas)

max(AIC_70$Radon)
max(AIC_70$Gas)

ggplot(AIC_70)+geom_histogram(aes(x=Radon))

AIC_top_10<-order(AIC_70, AIC_70$AIC)


