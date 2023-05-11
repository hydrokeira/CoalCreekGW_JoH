#install.packages(c("gplots", "plotrix"))

library(gplots) # for plotCI
library(plotrix) # for axis.break
require(dplyr)
require(lubridate)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

isoall<-read.csv("Iso_All.csv")

precip_weightedavg<-data.frame(c("rain", "snow"), 
                               c("-6.95", "-19.56"), 
                               c("-44.64", "-143.49"))

colnames(precip_weightedavg)<-c("type", "d18O", "d2H")

precip_weightedavg$d18O<-as.numeric(precip_weightedavg$d18O)
precip_weightedavg$d2H<-as.numeric(precip_weightedavg$d2H)
precip_weightedavg$LMWL<-precip_weightedavg$d18O*7.4+2.68

iso_LT<-subset(isoall, isoall$type=="stream")

springs<-subset(isoall, isoall$type=="spring")

iso_LT$month<-as.character(month(iso_LT$Date))

myColors <- ifelse(iso_LT$month=="6" , "deepskyblue", 
                   ifelse(iso_LT$month=="7", "seagreen3",
                        ifelse(iso_LT$month=="8", "goldenrod",
                            ifelse(iso_LT$month=="9", "sienna2",
                                 "darkred"))))

pdf("DualIsotopeCoalCreek.pdf", width = 8.5, height = 7.5)

jpeg("DualIsotopeCoalCreek.jpeg", width = 800, height = 700)

tiff("DualIsotopeCoalCreek.tiff", units = "in", width = 8.5, height = 7.5, res = 300)

par(mar=c(5,5,5,5))

gplots::plotCI(
  x=iso_LT$d18O, # note that since there's a plotrix::plotCI and a gplots::plotCI, we must be explicit
  y=iso_LT$d2H,
  pch=19, # 24 = filled triangle pointing up
  cex=1, # symbol size multiplier,
  col=myColors,
  xlim=c(-21,-10), # x axis limits
  xaxp=c(-21,-10,12), # x-min tick mark, x-max tick mark, number of intervals between tick marks
  xaxs="i", # makes the axis just fit the range, rather than extending 4% beyond
  xlab=expression(paste(delta^{18}, "O (\u2030)")), # x axis label
  ylim=c(-155, -75), # NOTE AXIS FROM 0.4 TO 1.0
  yaxp=c(-155, -75, 14),
  yaxs="i", # makes the axis just fit the range, rather than extending 4% beyond
  axes=FALSE,
  ylab=expression(paste(delta^{2}, "H (\u2030)")),
  las=1, # axis labels horizontal (default is 0 for always parallel to axis)
  font.lab=2, # 1 plain, 2 bold, 3 italic, 4 bold italic, 5 symbol
  family="Times"
)

lines(x=iso_sums_tot_crop$dO,
      y=iso_sums_tot_crop$LMWL,
      axes=FALSE)

gplots::plotCI(
  x=iso_LT$d18O, # note that since there's a plotrix::plotCI and a gplots::plotCI, we must be explicit
  y=iso_LT$d2H,
  pch=19, # 24 = filled triangle pointing up
  cex=1, # symbol size multiplier,
  col=myColors,
  xlim=c(-21,-10), # x axis limits
  xaxp=c(-21,-10,12), # x-min tick mark, x-max tick mark, number of intervals between tick marks
  xaxs="i", # makes the axis just fit the range, rather than extending 4% beyond
  xlab=expression(paste(delta^{18}, "O (\u2030)")), # x axis label
  ylim=c(-155, -75), # NOTE AXIS FROM 0.4 TO 1.0
  yaxp=c(-155, -75, 14),
  yaxs="i", # makes the axis just fit the range, rather than extending 4% beyond
  axes=FALSE,
  ylab=expression(paste(delta^{2}, "H (\u2030)")),
  las=1, # axis labels horizontal (default is 0 for always parallel to axis)
  font.lab=2, # 1 plain, 2 bold, 3 italic, 4 bold italic, 5 symbol
  family="Times",
  add=TRUE
)

box() # draw four sides to the plot area
axis(
  side=1, # X axis
  at=c(-21, -20,-19,-18, -17, -16, -15, -14, -13, -12, -11, -10),
  labels=c("-21", "-20","-19", "-18","-17","-16", "-15","-14","-8", "-7", "-6", "-5"),
  family="Times"
)
axis(side=1,at=0,labels="surgery", tick=FALSE) # additional X axis label without a tick mark
axis(
  side=2, # Y axis
  las=1, # text horizontal
  at=c(-155, -150, -145, -140,-135,-130, -125, -120, -115, -110, -105, -100, -90, -85, -80, -75), # position of labels
  labels=c("-155", "-150","-145", "-140", "-135","-130","-125","-120","-115","-110",
           "-105","-100", "-50", "-45", "-40", "-35"), # NOTE ACTUAL LABELS: THE "0.0" IS FALSELY POSITIONED
  family="Times"
)
axis.break(axis=2,breakpos=-95,style="slash") # break the left Y axis
axis.break(axis=4,breakpos=-95,style="slash") # break the right Y axis
axis.break(axis=1,breakpos=-13.5,style="slash") # break the bottom X axis
axis.break(axis=3,breakpos=-13.5,style="slash") # break the top X axis
gplots::plotCI(
  x=springs$d18O,
  y=springs$d2H,
  pch=21, # symbol type 21 = filled circle
  pt.bg="white",
  cex=1.5,
  lty=1,
  add=TRUE # ADD this plot to the previous one
)

gplots::plotCI(
  x=-11.95,
  y=-84.64,
  uiw = 6.9046,
  err="y",
  pch=15, # symbol type 21 = filled circle
  pt.bg="black",
  cex=1.5,
  lty=1,
  add=TRUE # ADD this plot to the previous one
)

gplots::plotCI(
  x=-11.95,
  y=-84.64,
  uiw = 0.85,
  err="x",
  pch=15, # symbol type 21 = filled circle
  pt.bg="black",
  cex=1.5,
  lty=1,
  add=TRUE # ADD this plot to the previous one
)

gplots::plotCI(
  x=precip_weightedavg$d18O[2],
  y=precip_weightedavg$d2H[2],
  uiw=9.46,
  err="y",
  pch=17, # symbol type 21 = filled circle
  pt.bg="black",
  cex=1.5,
  lty=1,
  add=TRUE # ADD this plot to the previous one
)

gplots::plotCI(
  x=precip_weightedavg$d18O[2],
  y=precip_weightedavg$d2H[2],
  uiw=1.20,
  err="x",
  pch=17, # symbol type 21 = filled circle
  pt.bg="black",
  cex=1.5,
  lty=1,
  add=TRUE # ADD this plot to the previous one
)


op<-par(family="Times")

legend(
  x=-12.7,
  y=-116,
  box.lty=0, # line type to surround the legend box (0 for none)
  title = "Sample Type",
  legend=c("rain","snow","springs", "stream"), # sequence of text for the legend
  pch=c(15,17,21, 19), # sequence of point types for the legend
  pt.bg=c("black","black","white", "black"), # sequence of fill colours for the points
  pt.cex=c(1.5,1.5,1.5, 1.5) # sequence of line types for the legend
)

legend(
  x=-12.7,
  y=-133,
  box.lty=0, # line type to surround the legend box (0 for none)
  title = "Collection Month",
  legend=c("June","July","August", "September", "October"), # sequence of text for the legend
  pch=c(21,21,21,21,21), # sequence of point types for the legend
  pt.bg=c("deepskyblue","seagreen3","goldenrod", "sienna2", "darkred"), # sequence of fill colours for the points
  pt.cex=c(1.5,1.5,1.5,1.5) # sequence of line types for the legend
)

dev.off()

