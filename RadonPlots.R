require(schoolmath)
require(zoo)
require(EGRET)
require(scales)
require(ggpmisc)
require(ggpubr)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

#read in radon data
radon<-read.csv("Radon_CleanName.csv")
names(radon)[2]<-"Site"
radon<-radon[,-1]

#format radon data - set date as date, remove site number column, remove sites with missing data
radon$Date<-as.Date(radon$Date, "%m/%d/%y")
radon<-radon[complete.cases(radon),]

#radon_tot<-merge(radon, stream_meter, by="Site")

#radon_tot<-subset(radon_tot, radon_tot$Stream.Meter > 9000)

#create data frame of table of counts
radon_counts<-as.data.frame(table(radon$Site))

#name columns
colnames(radon_counts)<-c("Site", "Times Sampled")

#read in stream meter data
stream_meter<-read.csv("RadonStreamMeter.csv")
stream_meter<-stream_meter[,c(1,2)]

#merge table of counts with stream meter
radon_counts<-merge(radon_counts, stream_meter, by="Site")

#write csv
#write.csv(radon_counts, "RadonCounts.csv")

#create list of sites that have been sampled > 3 times "long term sites"
long_term<-radon_counts[radon_counts$`Times Sampled` > 3, ]

#extract names from this list
long_term_sites<-long_term$Site

#calculate mean
radon_mean<-aggregate(radon, by=list(radon$Site), FUN=mean)
radon_mean<-radon_mean[,c(1,4)]
names(radon_mean)<-c("Site", "Mean")

#calculate SD
radon_sd<-aggregate(radon, by=list(radon$Site), FUN=sd)
radon_sd<-radon_sd[,c(1,4)]
names(radon_sd)<-c("Site", "SD")

#merge counts, mean, and sd data into one file
radon_tot<-merge(radon_counts, radon_mean, by="Site")
radon_tot<-merge(radon_tot, radon_sd, by="Site")

#round numbers to 2 decimal places
radon_tot[,c(4,5)]<-lapply(radon_tot[,c(4,5)],round, 2)

#write csv
write.csv(radon_tot, "RadonTot.csv")

radon_tot$min<-radon_tot$Mean-radon_tot$SD
radon_tot$max<-radon_tot$Mean+radon_tot$SD

#merge radon and stream meter
radon_merge<-merge(radon, stream_meter, by="Site")

#renames
radon<-radon_merge

#extract long term sites
radon_LT<-radon[radon$Site %in% long_term_sites, ]

#set new wd
setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Data")

#read in SNOTEL data
WYprecip<-read.delim("BUTTEsnotelUpdated.txt", sep = ",", skip = 63)
WYprecip<-WYprecip[,c(1,3)]
names(WYprecip)<-c("Date", "Precip")
WYprecip<-subset(WYprecip, WYprecip$Date > as.Date("2021-05-26") & WYprecip$Date < as.Date("2021-10-21"))

#calculate amount of daily precip from cummulative precip
amount<-diff(WYprecip$Precip, lag=1)

#take 5 day moving average
rollavg<-rollmean(amount, k=5)

#find difference between the WY df and rollavg list
WYremove<-length(WYprecip$Precip)-length(rollavg)

#remove that amount
WYprecip<-WYprecip[-c(1:WYremove),]

#find difference between the amount list and rollavg list
amountremove<-length(amount)-length(rollavg)

#remove that amount
amount<-amount[-c(1:amountremove)]

#merge together
WYprecip<-cbind(WYprecip, amount, rollavg)

#convert in to mm
WYprecip$amount<-WYprecip$amount*25.4
WYprecip$rollavg<-WYprecip$rollavg*25.4

#set date to as date
WYprecip$Date<-as.Date(WYprecip$Date)

#remove rows that are before/after sampling period
WYprecip<-WYprecip[-c(122:127),]

#read in summer discharge from NWIS
#summerQ<-readNWISDaily("09111250", "00060", "2021-05-26", "2021-10-21")

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek")

summerQ<-read.csv("CCregressedQ.csv")

summerQ$date<-as.Date(summerQ$date, "%m/%d/%y")

summerQ<-subset(summerQ, summerQ$date > "2021-05-25" & summerQ$date < "2021-10-22")

#take 5dma of Q
Qrollavg<-rollmean(summerQ$CC_m3_s, k=5, na.rm=T)

Qremove<-length(summerQ$CC_m3_s)-length(Qrollavg)

summerQ<-summerQ[-c(1:Qremove),]

summerQ<-cbind(summerQ, Qrollavg)

summerQ<-summerQ[complete.cases(summerQ$Q),]

radon_UP<-subset(radon_LT, radon_LT$Stream.Meter > 9100)
radon_DOWN<-subset(radon_LT, radon_LT$Stream.Meter < 9100)

radon_UP$Site<-factor(radon_UP$Site, levels = c("UPSTREAM","CC6", "CC7", "CC8", "DOWNSTREAM",
                                                "UPSTREAMELK", "COAL15"))

radon_UP$min<-radon_UP$time.corrected.Rn-radon_UP$Uncertainty
radon_UP$max<-radon_UP$time.corrected.Rn+radon_UP$Uncertainty

radon_UP$loc<-ifelse(radon_UP$Site %in% c("UPSTREAMELK", "COAL15"), "alluvial", "fracture")

lines<-c("alluvial"="dashed", "fracture"="solid")

pdf("Radon_Precip_SM_upperCC.pdf", width = 14, height = 7)

#radon and precip
p2<-ggplot()+geom_bar(WYprecip, mapping=aes(Date, amount), stat = "identity", fill="grey")+
  geom_line(radon_UP, mapping=aes(Date, time.corrected.Rn, col=Site, lty=loc), lwd=1)+
  geom_point(radon_UP, mapping=aes(Date, time.corrected.Rn, col=Site), size=2)+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*1, name=expression(paste(" "^{222},"Rn (piC/L)"))), 
                     name = "Daily Precipitation (mm)")+theme_classic()+
  geom_errorbar(radon_UP, mapping = aes(x=Date, ymin=min, ymax=max, col=Site, width=0))+
  theme(text = element_text(size = 20, family = "Times"))+scale_color_viridis_d(option = "viridis", direction = -1)+
  scale_linetype_manual(values = lines)+labs(lty="Associated Feature")

p2

tiff("Radon_Isotope_Precip_plot.tiff", units = "in", width = 12, height = 9, res = 300)

ggarrange(p2, p1, ncol = 1)

dev.off()

# ggplot()+geom_bar(WYprecip, mapping=aes(Date, rollavg), stat = "identity", fill="grey")+
#   geom_line(radon_UP, mapping=aes(Date, time.corrected.Rn, col=Site))+
#   geom_point(radon_UP, mapping=aes(Date, time.corrected.Rn, col=Site))+
#   geom_line(summerQ, mapping = aes(date, Qrollavg))+
#   scale_y_continuous(sec.axis = sec_axis(trans = ~.*1, name="Radon (pC/L)"), 
#                      name = "5 Day Moving Average Precipitation (mm) / 
#                      Discharge (cms)")+theme_classic()
# 
# ggplot()+geom_bar(WYprecip, mapping=aes(Date, rollavg), stat = "identity", fill="grey")+
#   geom_line(radon_DOWN, mapping=aes(Date, time.corrected.Rn, col=Site))+
#   geom_point(radon_DOWN, mapping=aes(Date, time.corrected.Rn, col=Site))+
#   scale_y_continuous(sec.axis = sec_axis(trans = ~.*1, name="Radon (pC/L)"), 
#                      name = "5 Day Moving Average Precipitation (mm) / 
#                      Discharge (cms)")+
#   geom_line(summerQ, mapping = aes(date, Qrollavg))+
#   theme_classic()+theme(axis.title.x = element_blank(),
#                         axis.text.x=element_blank(),
#                         axis.ticks.x=element_blank())

pdf("RadonConcentrationDate.pdf", width = 12, height = 8)

ggarrange(p2,p1, ncol=1)

dev.off()

ggplot(radon_merge, aes(rollavg, time.corrected.Rn))+geom_point(aes(col=Site))+
  geom_smooth(aes(col=Site), method = "lm", se=FALSE)+labs(x="5 Day Moving Average Precipitation (mm)",
                                                           y="Radon (pC/L)")+
  theme(text = element_text(size=20))+
  stat_poly_eq(aes(col=Site),position = "identity", formula = y~x, geom = "text_npc",
               output.type = "text", label.x="right")+ylim(c(0,35))


stream_meter$y<--.1
#remove dam names
stream_meter<-stream_meter[-c(21,25,15),]

pal<-colorRampPalette(c("deep sky blue", "orange","maroon1"))

LT_names<-stream_meter[stream_meter$Site %in% long_term_sites,]

remove_these_dates<-c(as.Date("2021-09-14"), as.Date("2021-09-28"),
                      as.Date("2021-08-17"), as.Date("2021-08-09"))

radon_UP_dates<-radon_UP[!(radon_UP$Date %in% remove_these_dates),]

pdf("AllRadon_wNames.pdf", width = 11, height = 5)

ggplot()+
  geom_line(radon, mapping=aes(Stream.Meter, time.corrected.Rn,col=as.character(Date)))+
  geom_point(radon, mapping=aes(Stream.Meter, time.corrected.Rn, col=as.character(Date)))+
  xlab("Stream Meter")+ylab("222Rn (piC/L)")+
  labs(colour="Date")+scale_x_continuous(trans = "reverse", breaks = pretty_breaks(n=6))+
  theme(text=element_text(size=20))+
  geom_text(stream_meter, mapping=aes(x=as.numeric(Stream.Meter), y=y, label=Site),
            hjust=0, vjust=0, stat = "identity", check_overlap = FALSE, angle=90)+
  scale_color_manual(values=pal(16))+theme_classic()

pdf("LongTerm_RadonUp_SM.pdf", width = 10, height = 7)

ggplot()+
  geom_line(radon_UP_dates, mapping=aes(Stream.Meter, time.corrected.Rn,col=as.character(Date)), lwd = 1)+
  geom_point(radon_UP_dates, mapping=aes(Stream.Meter, time.corrected.Rn, col=as.character(Date)), size = 2)+
  xlab("Stream Meter")+ylab("222Rn (piC/L)")+
  labs(colour="Date")+scale_x_continuous(trans = "reverse", limits = c(12000,9000),
                                         breaks = pretty_breaks(n=6))+
  theme(text=element_text(size=20))+
  geom_text(LT_names, mapping=aes(x=as.numeric(Stream.Meter), y=y, label=Site),
            hjust=0, vjust=0, stat = "identity", check_overlap = FALSE, angle=90)+
  theme_classic()+theme(text = element_text(size = 20))+
  scale_color_viridis_d(option = "viridis", direction = -1)

dev.off()

radon_UP<-subset(radon_LT, radon_LT$Stream.Meter > 9000)
stream_meter_UP<-subset(stream_meter, stream_meter$Stream.Meter > 9000)

radon_UP<-radon_UP[radon_UP$Site %in% long_term_sites,]

stream_meter_UP<-stream_meter_UP[stream_meter_UP$Site %in% long_term_sites,]

pdf("LTRadonBoxplot_wNames.pdf", width = 9, height = 5)

ggplot()+
  geom_boxplot(radon_UP, mapping=aes(Stream.Meter, time.corrected.Rn, group=Site),
               outlier.shape = NA)+
  geom_jitter(radon_UP, mapping=aes(Stream.Meter, time.corrected.Rn), size=2, shape=1, col="grey47")+
  xlab("Stream Meter")+ylab("222Rn (piC/L)")+
  scale_x_continuous(trans = "reverse", breaks = pretty_breaks(n=6))+
  geom_hline(yintercept = 0, color="light grey")+
  geom_text(stream_meter_UP, mapping=aes(x=as.numeric(Stream.Meter), y=-3.5, label=Site),
            hjust=.5, vjust=-2.8, stat = "identity", check_overlap = FALSE, angle=45, family="Times")+
            #position = position_jitter(height = 2))+
  #geom_errorbar(radon_UP, mapping = aes(x=Stream.Meter, ymin=min, ymax=max),
                #width=50)+
  theme_classic()+theme(text=element_text(size=20, family = "Times"))+xlim(c(12150,9000))

dev.off()

