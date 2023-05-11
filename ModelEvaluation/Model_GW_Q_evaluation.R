#evaluate GW discharge across different model runs
require(dplyr)
require(reshape2)
require(ggplot2)
require(data.table)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

GW_csv_list<-list.files(path = ".", pattern = "GWDischarge")
GW_csv_list<-GW_csv_list[grep("medMC.csv", GW_csv_list)]

date_list<-c("06-23","06-29", "07-14", "07-20", "08-02", "08-30")

Q_csv_list<-list.files(path = ".", pattern = "Discharge")
Q_csv_list<-Q_csv_list[grep("medMC.csv", Q_csv_list)]
Q_csv_list<-setdiff(Q_csv_list, GW_csv_list)

#csv_list<-csv_list[-grep("MC2", csv_list)]

csv_names<-gsub(".csv.*","",Q_csv_list)
colnames_list<-gsub(".*exported_", "", csv_names)

gw_df <- GW_csv_list %>% 
  lapply(read.csv) %>% 
  as.data.frame %>%
  bind_cols

gw_df$run<-seq(1,10000,1)
colnames(gw_df)<-c(date_list, "run")

# gw_df_median<-lapply(gw_df[,c(1:6)], median)
# gw_df_mean<-lapply(gw_df[,c(1:6)], mean)

gw_df_melt<-melt(gw_df, id.vars = c("run"))

model_input<-read.csv("RadonModelInput_FINAL.csv")

Qinput<-as.data.frame(model_input %>%
  group_by(Date) %>%
  summarise(min_Q=min(modeledQ, na.rm=T), max_Q=max(modeledQ, na.rm=T)))

is.na(Qinput) <- sapply(Qinput, is.infinite)

Qinput<-Qinput[complete.cases(Qinput),]

Qinput$diff<-Qinput$max_Q-Qinput$min_Q

Qinput$variable<-as.character(format(as.Date(Qinput$Date, "%m/%d/%y"), "%m-%d"))

Q_GW<-merge(Qinput, gw_df_melt, by="variable")

model_input<-model_input[6,c(17:18)]

upstream_streamlength<-model_input$Distance.Downstream.from.Upstream/10000
  
cc6_streamlength<-model_input$Distance.Downstream.from.CC6/10000

gw_df_melt$start<-ifelse(gw_df_melt$variable=="08-30", 11795, 11956)

gw_df_melt<-gw_df_melt %>%
  mutate(SM_step = ifelse(gw_df_melt$variable=="08-30", cc6_streamlength, 
                          upstream_streamlength))

gw_df_melt$SM_cumsum<-gw_df_melt$run*gw_df_melt$SM_step

gw_df_melt$SM<-gw_df_melt$start-gw_df_melt$SM_cumsum

gw_df_melt$GW_sum

pdf("GW_Flux_SM.pdf", width = 8, height = 5)

ggplot(gw_df_melt, aes(SM, value, col=variable))+
  geom_vline(xintercept=9600, size=1, linetype="dashed")+
  #geom_vline(xintercept =9200, size=1, linetype="dashed")+
  geom_line(lwd=1)+scale_x_continuous(trans = "reverse")+
  theme_classic()+scale_color_viridis_d(direction = -1)+theme(text = element_text(size = 20, family = "Times"))+
  labs(x="Stream Meter (m)", y="Incoming GW Flux (m/s)", col="Date")

dev.off()

# gw_df_melt$month<-gsub("-.*", "", gw_df_melt$variable)
# 
# cols<-c("06"="honeydew3", "07"="goldenrod3", "08"="darkorange3")
# 
# ggplot(gw_df_melt, aes(SM, value))+
#   geom_vline(xintercept=9600, size=1, linetype="dashed")+
#   geom_vline(xintercept =9200, size=1, linetype="dashed")+
#   geom_line(aes(col=month, group=variable),lwd=1)+scale_x_continuous(trans = "reverse")+
#   theme_classic()+scale_color_manual(values = cols)+
#   theme(text = element_text(size = 20))+
#   labs(x="Stream Meter (m)", y="Incoming GW Flux (m/s)", col="Date")

ref_table<-read.csv("RadonModelInput_FINAL.csv")

width_table<-ref_table %>%
  group_by(Date) %>%
  summarise(mean(modeledWidth, na.rm=TRUE))

width_table<-width_table[complete.cases(width_table),]
width_table$Date<-date_list
width_table$seg_len<-ifelse(width_table$Date=="08-30", cc6_streamlength, 
                            upstream_streamlength)
width_table$scale_fact<-width_table$`mean(modeledWidth, na.rm = TRUE)`*width_table$seg_len
colnames(width_table)[1]<-"variable"

gw_df_melt<-merge(gw_df_melt, width_table[,c(1,4)], by=c("variable"))

gw_df_melt$scaled_value<-gw_df_melt$value*gw_df_melt$scale_fact

gw_df_melt<-gw_df_melt %>%
  group_by(variable) %>%
  mutate(GW_sum = cumsum(scaled_value))

gw_df_melt$type<-ifelse(gw_df_melt$SM > 9600, "fracture", "alluvial")

GW_cont<-gw_df_melt %>%
  dplyr::group_by(variable, type) %>%
  summarise(GW_vol=max(GW_sum))

GW_prop<-as.data.frame(t(dcast(GW_cont, formula = type ~ variable)))

colnames(GW_prop)<-c("tot", "fracture")

GW_prop<-GW_prop[-1,]

GW_prop<-as.data.frame(sapply(GW_prop[,c(1,2)], as.numeric))

GW_prop$alluvial<-GW_prop$tot-GW_prop$fracture

GW_prop$variable<-date_list

GW_cont<-merge(GW_prop, Qinput, by="variable")

GW_cont$tot_prop<-GW_cont$tot/GW_cont$diff

GW_cont$Alluvial<-GW_cont$alluvial/GW_cont$diff

GW_cont$Fracture<-GW_cont$fracture/GW_cont$diff

colnames(GW_cont)[1]<-"month"

GW_cont_melt<-melt(GW_cont, id.vars=c("month", "tot", "fracture", "alluvial", "Date", "min_Q", "max_Q", "diff", "tot_prop"))

GW_cont_melt$variable<-factor(GW_cont_melt$variable, levels = c("Fracture", "Alluvial"))

GW_cont$alluvial_gwprop<-GW_cont$alluvial/GW_cont$tot

GW_cont$fracture_gwprop<-GW_cont$fracture/GW_cont$tot

GW_cont$fracture_alluvial<-GW_cont$fracture/GW_cont$alluvial

#write.csv(GW_cont, "ModelOutputGWSummary.csv")


# p1<-ggplot(GW_cont, aes(variable,GW_prop))+geom_bar(stat = "identity", aes(col=variable, fill=variable))+
#   theme_classic()+scale_color_viridis_d(direction = -1)+scale_fill_viridis_d(direction = -1)+
#   labs(x="Date", y="Groundwater Proportion", col="Date", fill="Date")+
#   theme(text = element_text(size = 20, family = "Times"), axis.text.x = element_text(angle=45, hjust = 1))

p1<-ggplot(GW_cont_melt, aes(month,value))+geom_bar(stat = "identity", aes(col=variable, fill=variable))+
  theme_classic()+scale_color_viridis_d(direction = -1)+scale_fill_viridis_d(direction = -1)+
  labs(x="Date", y="Groundwater Proportion", col="", fill="")+
  theme(text = element_text(size = 20, family = "Times"), axis.text.x = element_text(angle=45, hjust = 1))

pdf("Cumulative_Prop_GW_discharge.pdf", width = 15, height = 6)

p2<-ggplot(gw_df_melt, aes(SM, GW_sum, col=variable))+
  geom_vline(xintercept=9600, size=1, linetype="dashed")+
  #geom_vline(xintercept =9200, size=1, linetype="dashed")+
  geom_line(lwd=1)+scale_x_continuous(trans = "reverse")+
  theme_classic()+scale_color_viridis_d(direction = -1)+theme(text = element_text(size = 20, family = "Times"))+
  labs(x="Stream Meter (m)", y="Cumulative Groundwater Discharge (cms)", col="Date")

dev.off()

pdf()

pdf("Cumulative_Prop_GW_discharge_wbarplot.pdf", width = 15, height = 6)

ggarrange(p2, NULL, p1, widths = c(0.6, 0.02, 0.45), nrow = 1)

dev.off()
# 
# q_df<- Q_csv_list %>% 
#   lapply(read.csv) %>% 
#   as.data.frame %>%
#   bind_cols
# 
# colnames(q_df)<-date_list
# 
# q_df$run<-seq(1,10000,1)
# q_df_melt<-melt(q_df, id.vars = "run")
# 
# df_merged<-merge(gw_df_melt, q_df_melt, by=c("variable", "run"))
# df_merged$gw_prop<-df_merged$GW_sum/df_merged$value.y
# 
# pdf("GW_proportion.pdf", width = 8, height = 5)
# 
# ggplot(df_merged, aes(SM, gw_prop, col=variable))+geom_line(lwd=1)+scale_x_continuous(trans = "reverse")+
#   theme_classic()+scale_color_viridis_d(direction = -1)+theme(text = element_text(size = 20))+
#   labs(x="Stream Meter (m)", y="Groundwater Proportion", col="Date")
# 
# dev.off()
# 
# 
# colnames(gw_df)<-colnames_list

# #i=1
# 
# gw_prop<-data.frame(matrix(ncol = length(colnames_list), nrow = 10000))
# 
# for (i in 1:length(colnames_list)) {
#   
#   gw_prop[,i]<-(gw_df[,colnames_list[i]]/q_df[,colnames_list[i]])*100
#   
# }
# 
# colnames(gw_prop)<-colnames_list
# 
# df<-gw_prop
# 
# gw_prop$run<-seq(1,10000,1)
# q_df$run<-seq(1,10000,1)
# q_df_melt<-melt(q_df, id.vars = "run")
# 
# gw_prop_melt<-melt(gw_prop, id.vars = "run")
# 
# gw_df$run<-seq(1,10000,1)
# gw_df_melt<-melt(gw_df, id.vars = "run")
# 
# gw_df_melt <- gw_df_melt %>%
#   group_by(variable) %>%
#   mutate(roll_sum = cumsum(value))
# 
# gw_prop_melt <- gw_prop_melt %>%
#   group_by(variable) %>%
#   mutate(roll_sum = cumsum(value))
# 
# gw_prop<-as.data.frame(gw_df_melt$roll_sum/q_df_melt$value)
# gw_prop<-cbind(gw_df_melt[,c(1,2)],gw_prop)
# names(gw_prop)[3]<-"GWProp"
# 
# ggplot(gw_prop, aes(run, GWProp, col=variable))+geom_line(lwd=1)+theme_classic()+
#   scale_x_continuous()+xlab("Stream Meter")+ylab("GW Proportion")+
#   theme(text = element_text(size = 20))+scale_color_viridis_d(direction = -1)+labs(col="Date")
# 
# ggplot()+geom_line(gw_df_melt, mapping = aes(run, roll_sum, col=variable), lty="dashed")+
#   geom_line(q_df_melt, mapping = aes(run, value, col=variable))+theme_classic()
# 
