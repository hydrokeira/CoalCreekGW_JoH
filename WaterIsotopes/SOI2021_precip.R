##seasonal origin index##
##using WY 2021 precipitation##
#install.packages("schoolmath")
require(EflowStats)
require(ggplot2)
require(schoolmath)
require(dplyr)
require(zoo)
require(lubridate)
require(data.table)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek")
iso<-read.csv("MTCBPrecip.csv")
amount<-read.csv("Butte_Snotel_Precip.csv")

iso<-iso[complete.cases(iso$Processed.Delta.D),]

#format date
iso$Date<-as.Date(iso$Date, "%Y-%m-%d")
amount$Date<-as.Date(amount$Date, "%m/%d/%y")

#get water year
iso$WY<-get_waterYear(iso$Date)
amount$WY<-get_waterYear(amount$Date)

#subset to WY 2015
iso<-subset(iso, iso$WY=="2021")
amount<-subset(amount, amount$WY=="2021")

amount<-amount[,c(1:5)]

colnames(amount)<-c("Date", "SWE", "Precip_Accumulation", "Precip_Increment", "Precip_Snow_Adj")

#ggplot(amount)+geom_line(aes(Date, SWE))+geom_line(aes(Date, Precip_Accumulation))

#take average of days where duplicate samples taken
#iso<-aggregate(iso, by=list(iso$date, iso$WY, iso$Type), FUN = mean)

#select colomns and rename
iso<-iso[,c(2,3,5)]
names(iso)<-c("Date","d2H", "dO")

#merge amount and iso
iso_tot<-merge(amount, iso, by=c("Date"), all.x = TRUE)

#get list of numbers where there is isotope data, include 0 for next step
iso_dates_nums<-c(0, which(!is.na(iso_tot$dO)))

iso_dates<-iso_tot[c(iso_dates_nums),]$Date

#open list for sums
sum_list<-list()

#sum each group of amounts between isotope data measurements, append to list
for (i in 1:length(iso_dates_nums)) {
  
  sum_these<-iso_tot[c((iso_dates_nums[i]+1):iso_dates_nums[i+1]),]
  
  amount_sum<-sum(sum_these$Precip_Increment)
  
  sum_list[[i]]<-amount_sum
  
}

#create df from list, cbind to isotope data
sum_df<-as.data.frame(unlist(sum_list))
iso_sums<-cbind(iso_dates, sum_df)
names(iso_sums)<-c("Date", "sum")
iso_sums_tot<-merge(iso_sums, iso, by="Date")

iso_sums_tot<-iso_sums_tot[-1,]

#separate out rain and snow for SOI
rain<-subset(iso_sums_tot, iso_sums_tot$Date > as.Date("2021-05-07"))
snow<-subset(iso_sums_tot, iso_sums_tot$Date < as.Date("2021-05-07"))

mean(rain$d2H)
mean(snow$d2H)

iso_sums_tot_crop<-bind_rows(rain, snow)

iso_sums_tot_crop$LMWL<-iso_sums_tot_crop$dO*7.4+2.68


#create list of dfs for list
dfs<-list(iso_sums_tot_crop, rain, snow)

#open list for weighted averages
WA_list<-list()

#calculate WA for overall, rain, and snow for 2015 WY
for (i in 1:length(dfs)) {
  
  df<-dfs[[i]]
  
  df$percent<-df$sum/sum(df$sum)
  
  df$WA<-df$dO*df$percent
  
  WA_list[[i]]<-sum(df$WA)
  
}

weighted.average <- function(x, w){
  ## Sum of the weights 
  sum.w <- sum(w, na.rm = T)
  ## Sum of the weighted $x_i$ 
  xw <- sum(w*x, na.rm = T)
  
  ## Return the weighted average 
  return(xw/sum.w)
}

#get standard error mean, x is data, w is weights, both are vectors
weighted.se.mean <- function(x, w, na.rm = T){
  ## Remove NAs 
  if (na.rm) {
    i <- !is.na(x)
    w <- w[i]
    x <- x[i]
  }
  
  ## Calculate effective N and correction factor
  n_eff <- (sum(w))^2/(sum(w^2))
  correction = n_eff/(n_eff-1)
  
  ## Get weighted variance 
  numerator = sum(w*(x-weighted.average(x,w))^2)
  denominator = sum(w)
  
  ## get weighted standard error of the mean 
  se_x = sqrt((correction * (numerator/denominator))/n_eff)
  return(se_x)
}

rain$weight<-rain$sum/(sum(rain$sum))

snow$weight<-snow$sum/(sum(snow$sum))

weighted.se.mean(rain$d2H, rain$weight, na.rm = T)

weighted.se.mean(snow$d2H, snow$weight, na.rm = T)

#create touple of types of WA and WA number
WA_touple<-list(c("total", WA_list[1]), c("rain", WA_list[2]), c("snow", WA_list[3]))

##start code for SOI calculation#

#create function to get SOI
get_SOI<-function(stream_iso, total, winter, summer){
  
  SOI<-list()
  
  for (i in 1:length(stream_iso)) {
    
    if(stream_iso[i] > total){
      
      SOI[[i]] <- (stream_iso[i]-total)/(summer - total)
      
    } else {
      
      SOI[[i]] <- (stream_iso[i]-total)/(total - winter)
      
    }
    
  }
  
  SOI_df<-as.data.frame(unlist(SOI))
  
  return(SOI_df)
  
}

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

remove_these_dates<-c(as.Date("2021-09-14"), as.Date("2021-09-28"),
                      as.Date("2021-08-17"), as.Date("2021-08-09"))

remove_these_samples<-c("CC9", "CC9.5", "ELK")

CC_iso<-read.csv("Isotope_CleanNames.csv")
CC_iso<-CC_iso[,c(4:7)]
names(CC_iso)<-c("Name", "Date", "d18O", "d2H")
CC_iso$Date<-as.Date(CC_iso$Date, "%m/%d/%y")

SM<-read.csv("RadonStreamMeter.csv")
names(SM)[1]<-"Name"

CC_iso<-merge(CC_iso, SM, by="Name", all.x = TRUE)

CC_iso<-subset(CC_iso, CC_iso$Stream.Meter > 9000)

CC_iso<-CC_iso[!c(CC_iso$Name %like% "CC7."),]

CC_iso<-CC_iso[!c(CC_iso$Name %in% remove_these_samples),]

CC_iso<-CC_iso[!c(CC_iso$Date %in% remove_these_dates),]


#get spring rows and remove them
#springs<-grep("SPRING", CC_iso$Name)
#CC_iso<-CC_iso[-c(springs),]

stream_iso_vector<-CC_iso$d18O

SOI<-get_SOI(stream_iso_vector, WA_list[[1]], WA_list[[3]], WA_list[[2]])

CC_iso<-cbind(CC_iso, SOI)
names(CC_iso)[12]<-"SOI"
CC_iso$Date<-as.Date(CC_iso$Date)
CC_iso$pretty_date<-format(as.Date(CC_iso$Date), "%m-%d")

SOI_date<-CC_iso %>%
  group_by(Date) %>%
  summarise(median(SOI))

write.csv(SOI_date, "SOI_date.csv")

# CC_iso %>%
#   group_by(Name) %>%
#   summarise(sd(d18O))

pdf("SeasonalOriginIndex.pdf", width = 8.5, height = 7)

ggplot(CC_iso, aes(pretty_date,SOI))+
  geom_hline(yintercept = 0, col="gray")+
  geom_boxplot(outlier.size = 2, outlier.shape = NA)+
  geom_jitter(shape=1, size=2, col="grey47")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20, family="Times"))+labs(x="Date")+ylab("Seasonal Origin Index")+
  theme(plot.margin = margin(l=40))
  

dev.off()

# ##precipitation data - 5d ma
# setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Data")
# 
# #read in SNOTEL data
# WYprecip<-read.delim("BUTTEsnotelUpdated.txt", sep = ",", skip = 63)
# WYprecip<-WYprecip[,c(1,3)]
# names(WYprecip)<-c("Date", "Precip")
# WYprecip<-subset(WYprecip, WYprecip$Date > as.Date("2021-05-26") & WYprecip$Date < as.Date("2021-10-21"))
# 
# #calculate amount of daily precip from cummulative precip
# amount<-diff(WYprecip$Precip, lag=1)
# 
# #take 5 day moving average
# rollavg<-rollmean(amount, k=5)
# 
# #find difference between the WY df and rollavg list
# WYremove<-length(WYprecip$Precip)-length(rollavg)
# 
# #remove that amount
# WYprecip<-WYprecip[-c(1:WYremove),]
# 
# #find difference between the amount list and rollavg list
# amountremove<-length(amount)-length(rollavg)
# 
# #remove that amount
# amount<-amount[-c(1:amountremove)]
# 
# #merge together
# WYprecip<-cbind(WYprecip, amount, rollavg)
# 
# #convert in to mm
# WYprecip$amount<-WYprecip$amount*25.4
# WYprecip$rollavg<-WYprecip$rollavg*25.4
# 
# #set date to as date
# WYprecip$Date<-as.Date(WYprecip$Date)
# 
# #remove rows that are before/after sampling period
# WYprecip<-WYprecip[!(WYprecip$rollavg < 0),]
# 
# CC_iso<-merge(CC_iso, WYprecip, by="Date")
# 
# pdf("SOI_precip.pdf", height = 10, width = 8)
# 
# ggplot(CC_iso)+geom_boxplot(aes(as.character(Date),SOI))+
#   geom_point(aes(as.character(Date), (rollavg)/10), col="blue", shape=9, size=4)+
#   geom_hline(yintercept = 0)+
#   scale_y_continuous(sec.axis = sec_axis(trans = ~.*10, name = "5dma Precipitation (mm)"))+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(x="Date")
# 
# dev.off()
# 
# LT<-as.data.frame(table(CC_iso$Name))
# LT<-subset(LT, LT$Freq>3)
# 
# CC_iso_LT<-CC_iso[which(CC_iso$Name %in% LT$Var1),]
# CC_iso_LT$Name<-factor(CC_iso_LT$Name, levels =c("UPSTREAM", "CC6", "CC7", "CC8", "DOWNSTREAM", "UPSTREAMELK", "COAL15",
#                        "UPSTREAMCC11", "CC11", "COAL12"))
# 
# ggplot(CC_iso_LT)+geom_point(aes(as.character(Date),SOI))+
#   geom_point(aes(as.character(Date), (rollavg)/10), col="blue", shape=9, size=4)+
#   geom_hline(yintercept = 0)+facet_wrap(~Name)+
#   scale_y_continuous(sec.axis = sec_axis(trans = ~.*10, name = "5dma Precipitation (mm)"))+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(x="Date")
# 
# 
# 
# 
