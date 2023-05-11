##plot 10 best AIC
setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

GW_csv_list<-list.files(path = ".", pattern = "best10AIC.csv")
GW_csv_list<-GW_csv_list[grep("GW", GW_csv_list)]

# Convert nasty filename to nice column name
translate_files_to_columns <- function(csv_filename, idx) {
  #find the first '_' after the string exported
  # or just the first character after 'exported_'
  gw_date<-sapply(strsplit(GW_csv_list, 'exported_'), `[`, 2)
  date<-substr(gw_date, 1, 6)
  return(date)
  
}

gw_df <- GW_csv_list %>% 
  lapply(read.csv) %>% 
  as.data.frame %>%
  bind_cols

colnames(gw_df)<-paste0(translate_files_to_columns(GW_csv_list), "_Run_", seq(1, ncol(gw_df)))

gw_df$run<-seq(1,10000,1)

gw_df_melt<-melt(gw_df, id.vars = "run")

gw_df_melt$date<-substr(gw_df_melt$variable, 1, 6)

model_input<-read.csv("RadonModelInput_FINAL.csv")

model_input<-model_input[6,c(17:18)]

upstream_streamlength<-model_input$Distance.Downstream.from.Upstream/10000

cc6_streamlength<-model_input$Distance.Downstream.from.CC6/10000

gw_df_melt$start<-ifelse(gw_df_melt$date=="083021", 11795, 11956)

gw_df_melt<-gw_df_melt %>%
  mutate(SM_step = ifelse(gw_df_melt$date=="083021", cc6_streamlength, 
                          upstream_streamlength))

gw_df_melt$SM_cumsum<-gw_df_melt$run*gw_df_melt$SM_step

gw_df_melt$SM<-gw_df_melt$start-gw_df_melt$SM_cumsum

pdf("AIC_top10_GWflux.pdf", width = 16, height = 10)

ggplot(gw_df_melt, aes(SM, value))+geom_line()+facet_wrap(~date)+theme_bw()+
  labs(y="Incoming GW Flux (m/s)", x="Stream Meter (m)")+
  scale_x_continuous(trans = "reverse")+theme(text = element_text(size = 20))
dev.off()



