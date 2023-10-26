require(stringr)

setwd("/Users/keirajohnson/CoalCreekGW_HESS/GWFluxModel/FinalRuns")

MC_files<-list.files(path = ".", pattern = "top5AIC")

MC_files<-MC_files[!MC_files %like% "95"]

Q<-MC_files[MC_files %like% "Discharge"]

Q<-Q[!Q %like% "GWDischarge"]

Radon<-MC_files[MC_files %like% "Radon"]

dates<-str_split(Radon, 'exported_', simplify = TRUE)[,2]
dates<-sub("\\_top5AICmed.csv*", "", dates)

Q_list<-list()

Rn_list<-list()

for (i in 1:length(dates)) {
  
  Q_file<-read.csv(Q[i])
  
  Q_file$date<-dates[i]
  
  Q_file$run<-seq(1,10000, 1)
  
  Rn_file<-read.csv(Radon[i])
  
  Rn_file$date<-dates[i]
  
  Rn_file$run<-seq(1,10000, 1)
  
  Q_list[[i]]<-Q_file
  
  Rn_list[[i]]<-Rn_file
  
}

Q_all<-bind_rows(Q_list)

Rn_all<-bind_rows(Rn_list)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

model_input<-read.csv("RadonModelInput_FINAL.csv")

model_input<-model_input[complete.cases(model_input$modeledQ),]

model_input$unique<-paste(model_input$NewSite, model_input$Date)

model_input<-model_input[!duplicated(model_input$unique),]

upstream_streamlength<-2848/10000

cc6_streamlength<-2687/10000

Q_meas<-model_input[c("Date","modeledQ","Distance.Downstream.from.Upstream", "time.corrected.Rn.x")]

Q_meas$Date<-format(as.Date(Q_meas$Date, "%m/%d/%y"), "%m/%d/%y")

Q_meas$date<-gsub("/", "", Q_meas$Date)

Q_meas$SM<-Q_meas$Distance.Downstream.from.Upstream

#apply to Q

Q_all$start<-ifelse(Q_all$date=="083021", 161, 0)

Q_all<-Q_all %>%
  mutate(SM_step = ifelse(Q_all$date=="083021", cc6_streamlength, 
                          upstream_streamlength))

Q_all$SM_cumsum<-Q_all$run*Q_all$SM_step

Q_all$SM<-Q_all$start+Q_all$SM_cumsum

Q_p_list<-list()

Q_letter_list<-c("a", "c", "e", "g", "i", "k")

date_list<-c("06/23/21", "06/29/21", "07/14/21", "07/20/21","08/02/21", "08/30/21")

for (i in 1:length(dates)) {
  
  Q_one<-subset(Q_all, Q_all$date==dates[i])
  
  Q_meas_one<-subset(Q_meas, Q_meas$date==dates[i])
  
  if(i==6) {
    
    Q_p_list[[i]]<-ggplot()+geom_point(Q_meas_one, mapping=aes(SM, modeledQ))+
      geom_line(Q_one, mapping=aes(SM, Discharge))+
      #scale_x_continuous(trans = "reverse")+
      theme_bw()+labs(x="Stream Distance (m)", y="Discharge (cms)", tag = Q_letter_list[i])+
      xlim(0, 3000)+
      theme(text = element_text(size=15))+
      geom_text(x=100, y=0.9*max(modeledQ), label=paste(date_list[i]))
    
  } else{
    
    Q_p_list[[i]]<-ggplot()+geom_point(Q_meas_one, mapping=aes(SM, modeledQ))+
      geom_line(Q_one, mapping=aes(SM, Discharge))+
      #scale_x_continuous(trans = "reverse")+
      theme_bw()+labs(x="", y="Discharge (cms)")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size=15))+xlim(0, 3000)+
      labs(tag = Q_letter_list[i])+
      geom_text(x=1000, y=0.9*max(Q_meas_one$modeledQ), label=paste(date_list[i]))
    
  }
  
  
  
  
}

Q_half<-ggarrange(plotlist=Q_p_list, ncol=1, nrow=6, align="v")

Q_half

#apply to Rn

Rn_all$start<-ifelse(Rn_all$date=="083021", 161, 0)

Rn_all<-Rn_all %>%
  mutate(SM_step = ifelse(Rn_all$date=="083021", cc6_streamlength, 
                          upstream_streamlength))

Rn_all$SM_cumsum<-Rn_all$run*Rn_all$SM_step

Rn_all$SM<-Rn_all$start+Rn_all$SM_cumsum

Rn_p_list<-list()

Rn_letter_list<-c("b", "d", "f", "h", "j", "l")

for (i in 1:length(dates)) {
  
  Rn_one<-subset(Rn_all, Rn_all$date==dates[i])
  
  Q_meas_one<-subset(Q_meas, Q_meas$date==dates[i])
  
  if(i==6) {
    
    Rn_p_list[[i]]<-ggplot()+geom_point(Q_meas_one, mapping=aes(SM, time.corrected.Rn.x))+
      geom_line(Rn_one, mapping=aes(SM, C))+
      #scale_x_continuous(trans = "reverse")+
      theme_bw()+labs(x="Stream Distance (m)", y=expression(paste(" "^{222},"Rn (piC/L)")),tag = Rn_letter_list[i])+
      xlim(0, 3000)+
      theme(text = element_text(size=15))
    
  } else{
    
    Rn_p_list[[i]]<-ggplot()+geom_point(Q_meas_one, mapping=aes(SM, time.corrected.Rn.x))+
      geom_line(Rn_one, mapping=aes(SM, C))+
      #scale_x_continuous(trans = "reverse")+
      theme_bw()+labs(x="", y=expression(paste(" "^{222},"Rn (piC/L)")))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size=15))+xlim(0, 3000)+
      labs(tag = Rn_letter_list[i])
    
  }
  
  
  
  
}

Rn_half<-ggarrange(plotlist=Rn_p_list, ncol=1, nrow=6, align="v")

date1 <- ggplot() + 
  geom_text(aes(x=0, y=0, label = "06/23/21"), 
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()
date2 <- ggplot() + 
  geom_text(aes(x=0, y=0, label = "06/29/21"), 
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()
date3 <- ggplot() + 
  geom_text(aes(x=0, y=0, label = "07/14/21"), 
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()
date4 <- ggplot() + 
  geom_text(aes(x=0, y=0, label = "07/20/21"), 
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()
date5 <- ggplot() + 
  geom_text(aes(x=0, y=0, label = "08/02/21"), 
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()
date6 <- ggplot() + 
  geom_text(aes(x=0, y=0, label = "08/30/21"), 
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()

plot_dates<-ggarrange(date1, date2, date3, date4, date5, date6, nrow = 6, ncol = 1, align = "v")


pdf("Rn_Q_ModelFit.pdf", width = 14, height = 12, family = "Times")

ggarrange(plot_dates, Q_half, Rn_half, ncol = 3, nrow=1, widths = c(0.1, 0.45, 0.45))

dev.off()




