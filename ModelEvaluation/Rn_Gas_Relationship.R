setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/Coal_Creek/Summer")

MC<-read.csv("MC_min_max_values.csv")

MC$flow<-c(0.28, 0.29, 0.26, 0.12, 0.2, 0.08)

pdf("Rn_Gas_relationship.pdf", height = 5.5, width = 7)

ggplot(MC, aes(med_rn, med_gas))+geom_point(aes(col=flow), size=4)+theme_classic()+
  labs(x="Median Radon Concentration (piC/L)", y="Median Gas Exchange (m/d)", col="Discharge")+
  theme(text = element_text(size=20))

dev.off()
