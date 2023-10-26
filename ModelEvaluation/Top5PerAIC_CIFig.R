setwd("/Users/keirajohnson/CoalCreekGW_HESS/GWFluxModel/FinalRuns")

MC_files<-list.files(path = ".", pattern = "top5AIC")

dates<-c("0623", "0629", "0714", "0720", "0802", "0830")

nice_dates<-c("06/23/21", "06/29/21", "07/14/21", "07/20/21", "08/02/21", "08/30/21")

i=3

upstream_streamlength<-2848/10000

cc6_streamlength<-2687/10000

#Q_all$SM_cumsum<-Q_all$run*Q_all$SM_step

#Q_all$SM<-Q_all$start-Q_all$SM_cumsum

letters_tag<-c("a", "b", "c", "d", "e", "f")

plot_list<-list()

for (i in 1:length(dates)) {
  
  one_date<-MC_files[MC_files %like% dates[i]]
  
  q<-one_date[one_date %like% "GWDischarge"]
  
  lower<-read.csv(q[1])
  
  med<-read.csv(q[2])
  
  upper<-read.csv(q[3])
  
  all_runs<-bind_cols(lower, med, upper)
  
  colnames(all_runs)<-c("lower", "med", "upper")
  
  all_runs$run<-seq(1,10000,1)
  
  all_runs$start<-ifelse(dates[i]=="0830", 161, 0)
  
  all_runs<-all_runs %>%
    mutate(SM_step = ifelse(dates[i]=="0830", cc6_streamlength, 
                            upstream_streamlength))
  
  all_runs$SM_cumsum<-all_runs$run*all_runs$SM_step
  
  all_runs$SM<-all_runs$start+all_runs$SM_cumsum
  
  if(i==1){
    
    plot_list[[i]]<-ggplot(all_runs, aes(x=SM))+
      geom_line(aes(y=lower), lwd=0)+
      geom_line(aes(y=upper), lwd=0)+
      geom_ribbon(aes(ymin=lower, ymax=med), fill="grey")+
      geom_ribbon(aes(ymin=med, ymax=upper), fill="grey")+
      geom_line(aes(y=med))+
      theme_bw()+
      ggtitle(nice_dates[i])+
      theme(text = element_text(size=20))+
      labs(x="", y="Groundwater Flux (m/s)", tag = letters_tag[i])+
      xlim(0,3000)+
      ylim(0,5.1E-5)
    
  }else if(i %in% c(2,3)){
    
    plot_list[[i]]<-ggplot(all_runs, aes(x=SM))+
      geom_line(aes(y=lower), lwd=0)+
      geom_line(aes(y=upper), lwd=0)+
      geom_ribbon(aes(ymin=lower, ymax=med), fill="grey")+
      geom_ribbon(aes(ymin=med, ymax=upper), fill="grey")+
      geom_line(aes(y=med))+
      theme_bw()+
      ggtitle(nice_dates[i])+
      theme(text = element_text(size=20))+
      labs(x="", y="",tag = letters_tag[i])+
      xlim(0,3000)+
      ylim(0,5.1E-5)
    
  }else if(i==4){
    
    plot_list[[i]]<-ggplot(all_runs, aes(x=SM))+
      geom_line(aes(y=lower), lwd=0)+
      geom_line(aes(y=upper), lwd=0)+
      geom_ribbon(aes(ymin=lower, ymax=med), fill="grey")+
      geom_ribbon(aes(ymin=med, ymax=upper), fill="grey")+
      geom_line(aes(y=med))+
      theme_bw()+
      ggtitle(nice_dates[i])+
      theme(text = element_text(size=20))+
      labs(x="Stream Distance (m)", y="Groundwater Flux (m/s)", tag = letters_tag[i])+
      xlim(0,3000)+
      ylim(0,5.1E-5)
    
  }else{
    
    plot_list[[i]]<-ggplot(all_runs, aes(x=SM))+
      geom_line(aes(y=lower), lwd=0)+
      geom_line(aes(y=upper), lwd=0)+
      geom_ribbon(aes(ymin=lower, ymax=med), fill="grey")+
      geom_ribbon(aes(ymin=med, ymax=upper), fill="grey")+
      geom_line(aes(y=med))+
      theme_bw()+
      ggtitle(nice_dates[i])+
      theme(text = element_text(size=20))+
      labs(x="Stream Distance (m)", y="", tag = letters_tag[i])+
      xlim(0,3000)+
      ylim(0,5.1E-5)
    
  }
  
}

pdf("CI_minmaxCItop5AIC.pdf", width = 12, height = 8, family = "Times")

ggarrange(plotlist = plot_list)

dev.off()

