setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

old<-list.files(pattern = "GWDischarge_exported")
old<-old[grep("_medMC.csv", old)]

new<-list.files(pattern = "GWDischarge_exported")
new<-new[grep("_top5AICmed.csv", new)]

date_list<-as.numeric(gsub("\\D", "", old))

p_list<-list()

for (i in 1:length(old)) {
  
  old_1<-read.csv(old[i])
  
  new_1<-read.csv(new[i])
  
  p_list[[i]]<-ggplot()+geom_line(old_1, mapping=aes(y=GWDischarge, x=seq(1,10000)), col="red")+
    geom_line(new_1, mapping=aes(y=GWDischarge, x=seq(1,10000)), col="blue")+
    theme_bw()+ggtitle(date_list[i])
  
}

pdf("Old_New_MC_Compare.pdf", width = 15, height = 10)

ggarrange(plotlist = p_list, nrow=2, ncol=3)

dev.off()
