require(ggplot2)
require(schoolmath)
require(zoo)
require(EGRET)
require(scales)
require(ggpmisc)
require(ggpubr)
require(EflowStats)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

iso<-read.csv("Isotope_CleanNames.csv")
iso$Date<-as.Date(iso$Date, "%m/%d/%y")

SM<-read.csv("RadonStreamMeter.csv")
names(SM)[1]<-"Name"

springs_list<-grep("SPRING", iso$Name)

springs<-iso[c(springs_list),]

springs<-springs[,c(4:7)]

colnames(springs)<-c("Name", "Date", "d18O", "d2H")

#iso<-iso[-c(springs_list),]

#iso<-iso[complete.cases(iso),]

# setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Data")
# 
# mainQ<-read.csv("CoalCreekDischarge2021.csv")
# mainQ<-mainQ[,c(3,5)]
# mainQ$date<-as.Date(mainQ$date, "%m/%d/%y")
# names(mainQ)[1]<-"Date"
# 
iso<-merge(iso, SM, by="Name", all.x = TRUE)

iso<-subset(iso, iso$Stream.Meter > 9000)

iso<-iso[c(1,5,6,7,11,12,14)]
#iso<-merge(iso, mainQ, by="Date", all.x = TRUE)
#iso$CC_m3_s<-ifelse(is.na(iso$CC_m3_s), "0.09", iso$CC_m3_s)

names(iso)<-c("Name", "Date", "d18O", "d2H","d18Oprec", "d2Hprec", "SM")

iso$dexcess<-iso$d2H-8*iso$d18O

iso$LMWL<-iso$d18O*8+10

# iso$class<-ifelse(iso$SM > 9100, "fracture", "BMA")

counts<-as.data.frame(table(iso$Name))

long_term<-subset(counts, counts$Freq > 3)
#| Var1=="ELK")

long_term_names<-long_term$Var1

iso_LT<-iso[iso$Name %in% long_term_names,]

iso_LT$type<-"stream"

springs$type<-"spring"

iso_all<-bind_rows(iso_LT, springs)

write.csv(iso_all, "IsotopeData.csv")

write.csv(iso_all, "Iso_All.csv")

names(iso_LT)[2]<-"Date"

iso_LT %>%
  group_by(Name) %>%
  summarise(mean(d18O))

pdf("IsoBoxplot_withNames.pdf", , width = 9, height = 5)

ggplot(iso_LT, aes(SM, d18O))+geom_boxplot(aes(group=Name), outlier.size = 2)+
  geom_jitter(size=2, shape=1)+
  xlab("Stream Meter")+ylab("d18O (permil)")+
  scale_x_continuous(trans = "reverse", breaks = pretty_breaks(n=6))+
  geom_text(stream_meter_UP, mapping=aes(x=as.numeric(Stream.Meter), y=-17.5, label=Site),
            hjust=.5, vjust=-2.8, stat = "identity", check_overlap = FALSE, angle=45)+
  #position = position_jitter(height = 2))+
  #geom_errorbar(radon_UP, mapping = aes(x=Stream.Meter, ymin=min, ymax=max),
  #width=50)+
  theme_classic()+theme(text=element_text(size=20))+xlim(c(12150,9000))

dev.off()


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
WYprecip<-WYprecip[!(WYprecip$rollavg < 0),]

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

iso_LT$Name<-factor(iso_LT$Name, levels = c("UPSTREAM", "CC6", "CC7", "CC8", "DOWNSTREAM",
                                            "UPSTREAMELK", "COAL15"))

iso_LT$uncert_min<-iso_LT$d18O-0.1
iso_LT$uncert_max<-iso_LT$d18O+0.1

#iso_LT$loc<-ifelse(iso_LT$Name %in% c("UPSTREAMELK", "COAL15"), "alluvial", "fracture")

#colnames(iso_LT)[2]<-"Date"

pdf("Isotope_Precip_Plot.pdf", width = 12, height = 7)

p1<-ggplot()+
  geom_bar(WYprecip, mapping=aes(Date, amount/4), stat = "identity", fill="grey")+
  #geom_line(summerQ, mapping = aes(date, Qrollavg))+
  geom_line(iso_LT, mapping=aes(Date, d18O+18, col=Name), lwd=1)+
  geom_point(iso_LT, mapping=aes(Date, d18O+18, col=Name))+
  geom_errorbar(iso_LT, mapping = aes(Date, d18O+18, ymin=uncert_min+18, max=uncert_max+18, col=Name),
                width=0)+
  scale_y_continuous(labels = ~(.*4), sec.axis = sec_axis(trans = ~(.*1)-18, name=expression(paste(delta^{18}, "O (\u2030)"))), 
                     name = "Daily Precipitation (mm)")+
  theme_classic()+theme(text = element_text(size=20, family = "Times"))+labs(col="Site")+
  scale_color_viridis_d(direction = -1)+labs(tag = "b")

p1
dev.off()

p1<-ggplot()+
  geom_bar(WYprecip, mapping=aes(Date, amount), stat = "identity", fill="grey")+
  theme_classic()+theme(text = element_text(size=20, family = "Times"))+labs(y="Precipitation (mm)")

p2<-ggplot()+
  geom_line(iso_LT, mapping=aes(Date, d18O, col=Name), lwd=1)+
  geom_point(iso_LT, mapping=aes(Date, d18O, col=Name))+
  geom_errorbar(iso_LT, mapping = aes(Date, d18O, ymin=uncert_min, max=uncert_max, col=Name),
                width=0)+
  theme_classic()+theme(text = element_text(size=20, family = "Times"), legend.position = "null")+labs(col="Site")+
  scale_color_viridis_d(direction = -1)

ggarrange(p1, p2, ncol = 1)

dev.off()

pdf("Isotope_Precip_Plot_Dexcess.pdf", width = 15, height = 7)

#dexcess plot
ggplot()+
  geom_bar(WYprecip, mapping=aes(Date, rollavg*1), stat = "identity", fill="grey")+
  geom_line(summerQ, mapping = aes(date, Qrollavg))+
  geom_line(iso_LT, mapping=aes(Date, dexcess/2, col=Name))+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*2, name="d-excess (permil)"), 
                     name = "5 Day Moving Average Precipitation (mm) / 
                     Discharge (cms)")+
  theme(text = element_text(size=20))+theme_classic()

dev.off()

pdf("Isotope_SM_Plot_Dexcess.pdf", width = 13, height = 7)

ggplot()+geom_point(iso_LT, mapping=aes(SM, dexcess, col=format(as.Date(Date), "%m")), size=2)+theme_classic()+
  labs(x="Stream Meter", y="d-excess (per mil)", col="Month")+scale_x_continuous(trans = "reverse")

dev.off()

getwd()
plot.new()
