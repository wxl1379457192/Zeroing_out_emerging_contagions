library("rstan")
library("bayesplot")
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
source('code/main_code.R')
source('code/NPI_code.R')
dir ="dataset"
vacW = read.csv(paste0(dir,"/vaccinations_weighted.csv"),stringsAsFactors = F)
vac = read.csv(paste0(dir,"/vaccinations.csv"),stringsAsFactors = F)
vaccination.rate = province_vaccination(vacW,vac)

library(dplyr)
library(reshape2)
vacc = vaccination.rate%>%group_by(date)%>%summarise(Total_vaccination_rate = sum(total_vaccination_rate),
                                                     Fully_vaccination_rate = sum(fully_vaccination_rate),
                                                     Boost_dose_rate = sum(boost_dose_rate))
vacc = melt(vacc,id.vars=c("date"))
p = ggplot(vacc)+geom_col(aes(x=date,y=value,fill=variable),size=0.01,position = "identity",alpha = 0.35)+
  scale_y_continuous(limits = c(0,50),
               expand = expansion(mult = c(0,0)),
               breaks = seq(0,50,15))+
  scale_x_date(expand = expansion(mult = c(0,0)))+
  scale_fill_manual(values = c("Total_vaccination_rate"="#91D1C2CC",
                               "Fully_vaccination_rate"="#DC0000CC",
                               "Boost_dose_rate"="#3C5488CC"),
                    labels=c("Total_vaccination_rate"="Total vaccination rate",
                             "Fully_vaccination_rate"="Fully vaccination rate",
                             "Boost_dose_rate" = "Boost dose rate"))+
  theme(legend.position = c(0.15,0.8),
        legend.background = element_rect(
          fill = "transparent", # Ìî³äÉ«
          colour = "transparent", # ¿òÏßÉ«
          size = 1.5),
        legend.title=element_blank(),
        panel.grid = element_blank(),
        axis.line=element_line(color="grey90",size=0.5),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size =unit(9,"pt")),
        axis.title.x= element_blank(),
        axis.text.x = element_text(color="black",size = unit(9,"pt")),
        axis.title.y= element_blank(),
        axis.text.y = element_text(color="black",size = unit(9,"pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))

outdir = "Figures"
ggsave(paste0(outdir,"/Vaccination_FigureS2.pdf"),p, width=14, height=8, units="cm", scale=1)

  
  