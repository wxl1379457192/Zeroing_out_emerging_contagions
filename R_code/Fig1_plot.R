library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(reshape2)
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
indir = "pre-5dataset"
casesISO = read.csv("dataset/New_cases_all.csv",stringsAsFactors = F) 
dat = read.csv(paste0(indir,"/Rt&smooth_NPI.csv"),stringsAsFactors = F)
dat = dat[,c("Date","city","New_cases_smooth")]
dat$Date <-as.Date(dat$Date) 
dat = unique(dat)
cases = casesISO[,c("Date","city","New_cases")]
cases$Date <-as.Date(cases$Date) 
dat = merge(dat,cases,by=c("Date","city"),all=T)
dat$New_cases[is.na(dat$New_cases)] = 0
dat$New_cases_smooth[is.na(dat$New_cases_smooth)] = 0
dat = dat%>%group_by(Date)%>%summarise(New_cases = sum(New_cases))
dat = melt(dat,id.vars=c("Date"))

p1 = ggplot()+geom_col(data=dat,aes(x=Date,y=value,fill=variable,alpha=variable),
                  position = "identity")+
  scale_fill_manual(values = c("New_cases"="#a3cd5b"),
                    labels=c("New_cases"="Daily reported new cases"))+
  scale_alpha_manual(values = c("New_cases"=0.4),
                    labels=c("New_cases"="Daily reported new cases"))+
  scale_y_sqrt(limits = c(0,31000),
               expand = expansion(mult = c(0,0)),
               breaks = c(100,3000,10000,20000))+
  scale_x_date(expand = expansion(mult = c(0,0)))+
  theme_grey()+
  annotate("text", x=as.Date("2020-08-20"),y=8000,label="SARS-CoV-2",size =unit(4,"pt"))+
  geom_segment(aes(x = as.Date("2020-04-01"), y = 10000, xend = as.Date("2021-01-17"), yend = 10000), 
               colour = "#024700", linetype = 7,size=0.4, 
               arrow= arrow(ends = "both",length = unit(0.03, "npc"))) +
  
  annotate("text", x=as.Date("2021-03-20"),y=25000,label="Alpha",size =unit(4,"pt"))+
  geom_segment(aes(x = as.Date("2021-01-18"), y = 28000, xend = as.Date("2021-05-20"), yend = 28000), 
             colour = "#024700", linetype = 7,size=0.4, 
             arrow= arrow(ends = "both",length = unit(0.03, "npc"))) +
  geom_vline(xintercept = as.Date("2021-01-17"),lty="dashed")+
  annotate("text", x=as.Date("2021-08-30"),y=25000,label="Delta",size =unit(4,"pt"))+
  geom_segment(aes(x = as.Date("2021-05-22"), y = 28000, xend = as.Date("2021-12-08"), yend = 28000), 
               colour = "#024700", linetype = 7,size=0.4, 
               arrow= arrow(ends = "both",length = unit(0.03, "npc"))) +
  geom_vline(xintercept = as.Date("2021-05-21"),lty="dashed")+
  geom_vline(xintercept = as.Date("2021-12-09"),lty="dashed")+
  annotate("text", x=as.Date("2022-03-10"),y=25000,label="Omicron",size =unit(4,"pt"))+
  geom_segment(aes(x = as.Date("2021-12-10"), y = 28000, xend = as.Date("2022-05-30"), yend = 28000), 
               colour = "#024700", linetype = 7,size=0.4, 
               arrow= arrow(ends = "both",length = unit(0.03, "npc"))) +
  theme(legend.position = c(0.15,0.8),
        legend.background = element_rect(
          fill = "transparent", # 填充色
          colour = "transparent", # 框线色
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

library(circlize)
library(tidyverse)
library(aplot)
all = read.csv(paste0(indir,"/Rt&smooth_NPI.csv"),stringsAsFactors = F)
all = all[which(all$Mean.R.>0),]
all$Control[which(all$Control>1)] = 1
NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")
auxname = c("Temp","fully_vaccination_rate")
library(dplyr)
all$Date = as.Date(all$Date)
all$province = substring(all$citycode,1,1)
ns = all[,c("VG","province",NPIname)]%>%group_by(VG,province)%>%
  summarise(L = mean(Lockdown),
            BPC = mean(Business_Premises_Closure),
            PTC = mean(Public_Transportation_Closure),
            GR = mean(Gathering_restriction),
            WC = mean(Workplace_Closure),
            SC = mean(School_Closure),
            MM = mean(Medicine_Management),
            MS = mean(Mass_screening),
            FM = mean(Facial_Mask),
            C = mean(Control))
library(reshape2)
ns = melt(ns,id.vars=c("VG","province"))
#my.col <-colorRampPalette(brewer.pal(11, "BrBG"))(140)
ns$group = "SD"
ns$group[which(ns$variable=="MS"|ns$variable=="MM")] = "PCR"
ns$group[which(ns$variable=="FM")] = "FM"
ns$group[which(ns$variable=="C")] = "C"
ns$group = factor(ns$group,ordered=TRUE,levels = c("SD","PCR","FM","C"))

fx=ns %>% mutate(Y = "NPIs")%>%
  ggplot()+
  geom_tile(aes(x=variable,y=Y,fill=group)) + theme_void()+
  labs(fill = "NPIs")+
  scale_fill_manual(values = c("SD"="#9B9CC9",
                               "PCR"="#B296C7",
                               "FM" ="#B3B3B3",
                               "C" = "grey50"))+
  theme(legend.position ="top")

ns$VG[which(ns$VG=="original&alpha")]="Pre-Delta"
ns$VG[which(ns$VG=="omicron")]="Omicron"
ns$VG[which(ns$VG=="delta")]="Delta"

ns$VG = factor(ns$VG,levels = c("Omicron","Delta","Pre-Delta"))
ns = ns[order(ns$VG),]

ns = do.call(rbind,lapply(seq(1,length(unique(ns$VG))),function(v){
  vall = subset(ns,ns$VG==unique(ns$VG)[[v]])
  vall$Region = paste0(vall$province,v)
  return(vall)
}))

ns$Region = factor(ns$Region,levels = c("61","62","63","51","52","53",
                                        "41","42","43","31","32","33",
                                        "21","22","23","11","12","13"))
nsy = ns[,c("Region","VG")]
fy= unique(nsy)%>% mutate(X = "Variant")%>%
  ggplot()+
  geom_tile(aes(y=Region,x=X,fill=VG),alpha=0.5) +
  theme_void()+
  scale_fill_manual(values = c("Pre-Delta" = "#99BAAC",
                               "Delta" = "#A49B90",
                               "Omicron" = "#257187"))+
  labs(fill = "")+ theme(legend.position ="bottom")


p = ggplot(ns,aes(y=Region,x=variable))+
  geom_tile(aes(fill=value))+
  scale_fill_gradient2(low="#0000ff", high="#365a87")+
  scale_y_discrete(position="right",) +
  theme_grey()+
  theme(legend.position ="bottom",
        legend.background = element_rect(
          fill = "transparent", # 填充色
          colour = "transparent", # 框线色
          size = 1.5),
        legend.title=element_blank(),
        axis.line=element_line(color="grey90",size=0.5),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size =unit(9,"pt")),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text.x = element_text(color="black",hjust = 1,
                                   vjust=0,angle = 90,size = unit(9,"pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.grid.major=element_blank(),
        axis.text.y= element_text(color="black", size =unit(9,"pt")))
library(aplot)
p2 = p %>%
  insert_top(fx, height = .02) %>%
  insert_left(fy, width = .12)


library(ggalt)
library(ggsci)
g = read.csv(paste0(indir,"/NPI&cases_smooth&Rt.csv"),stringsAsFactors = F)
g = g[,c("Date","city","original_start","VG","Mean.R.","I")]#,"Quantile.0.05.R","Quantile.0.95.R","variant"
g$Date = as.Date(g$Date)

g = do.call(rbind,lapply(split(g,g$city),function(c){
  c = do.call(rbind,lapply(split(c,c$original_start),function(c1){
    #if(unique(c1$VG=="original&alpha")){
    # c1 = c1[which(c1$I>=3),]
    #}
    c1$days = seq(1,nrow(c1),1)
    return(c1)
  }))
  return(c)
}))
g$VG = factor(g$VG,levels = c("omicron","delta","original&alpha"))
g = g[order(g$VG),]
g = do.call(rbind,lapply(seq(1,length(unique(g$VG))),function(v){
  vall = subset(g,g$VG==unique(g$VG)[[v]])
  vall$city = paste0(v,vall$city)
  return(vall)
}))

#g=melt(g,id.vars=c("Date","city"))

my.col <- c(colorRampPalette(brewer.pal(12,"Set3"))(130),colorRampPalette(brewer.pal(9, "Greys"))(27),
            colorRampPalette(brewer.pal(9,"Blues"))(12))
#my.col <-colorRampPalette(brewer.pal(11, "BrBG"))(140)
p3 = ggplot(g)+geom_line(aes(x=days,y=Mean.R.,color=paste(city,original_start),linetype=VG),
                       alpha=0.5,show.legend = FALSE,size=0.5)+
  #scale_color_viridis_c()+
  scale_linetype_manual(values = c("original&alpha" = 5,
                        "delta" = 2,
                        "omicron" = 3))+
  theme_grey()+ geom_hline(yintercept = 1,lty="dashed")+
  geom_smooth(aes(x=days,y=Mean.R.),alpha=0.4,span=0.2,color="#2D6587",
              fill="#87a0b7",method = "gam",size=0.9)+
  scale_color_manual(values=rev(my.col))+
  scale_y_sqrt(limits = c(0.1,11),breaks=c(0,1,3,6,10))+
  labs(x="Days",y="Instantaneous reproduction number (Rt)")+
  annotate("text", x=75,y=2.5,label="Rt = 1",size =unit(4,"pt"))+
  annotate("rect",xmin=70,xmax=80,ymin=2,ymax=3,
           color="#455787",fill="transparent",linetype=5) +
  geom_segment(aes(x = 75, y = 2, xend = 75, yend = 1), 
               colour = "#455787", linetype = 5,size=0.4, 
               arrow= arrow(length = unit(0.03, "npc"))) +
  scale_x_continuous(expand = expansion(mult = c(0,0)),limits=c(1,96),breaks = c(15,30,45,60,75,90))+
  theme(legend.position = c(0.15,0.8),
      legend.background = element_rect(
        fill = "transparent", # 填充色
        colour = "transparent", # 框线色
        size = 1.5),
      legend.title=element_blank(),
      panel.grid = element_blank(),
      axis.line=element_line(color="black",size=0.5),
      plot.title = element_text(color="black",hjust = 0,vjust=0,size =unit(9,"pt")),
      axis.title.x= element_text(color="black",hjust = 0.5,vjust=0,size =unit(9,"pt")),
      axis.text.x = element_text(color="black",size = unit(9,"pt")),
      axis.title.y= element_text(color="black",hjust = 0.5,vjust=0,size =unit(9,"pt")),
      axis.text.y = element_text(color="black",size = unit(9,"pt")),
      plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))

outdir = "FiguresV1129"
ggsave(paste0(outdir,"/fig1c.pdf"),p2, width=12, height=13, units="cm", scale=1)
ggsave(paste0(outdir,"/fig1a.pdf"),p1, width=16, height=8, units="cm", scale=1)
ggsave(paste0(outdir,"/fig1b.pdf"),p3, width=16, height=8, units="cm", scale=1)

f = read.csv(paste0(indir,"/Rt&smooth_NPI.csv"),stringsAsFactors = F)
fn = read.csv(paste0(indir,"/cityname.csv"),stringsAsFactors = F)
colnames(fn) = c("city","City")
f = merge(f,fn,by="city")

f = f[,c("Date","City","original_start","Label_q","VG",
         "Facial_Mask",
         "Mass_screening","Medicine_Management",
         "Control",
         "Lockdown","School_Closure","Workplace_Closure","Gathering_restriction",
         "Public_Transportation_Closure","Business_Premises_Closure")]

f = f%>%group_by(Date)%>%summarise(FM = mean(Facial_Mask),MC = mean(Mass_screening),
                               MM = mean(Medicine_Management),CC = mean(Control),
                               Lock = mean(Lockdown), SC = mean(School_Closure),
                               WC = mean(Workplace_Closure), GR = mean(Gathering_restriction),
                               PTC = mean(Public_Transportation_Closure),BPC = mean(Business_Premises_Closure)
                          )
f1 = melt(f,id="Date")
f1$Date = as.Date(f1$Date)
ggplot(f1)+geom_raster(aes(x=Date,y = variable,fill = value))+
  scale_x_date(expand = expansion(mult = c(0,0)),date_labels = "%Y %m")



