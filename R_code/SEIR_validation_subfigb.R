library(ggplot2)
library(openxlsx)
setwd("E:/zeroCOVID_NPI/Version0504")
SEIR =  read.xlsx("SEIR_validation/Validation_SEIR.xlsx", sheet = 1)

SEIR$VG[which(SEIR$VG=="delta")]="Delta"
SEIR$VG[which(SEIR$VG=="omicron")]="Omicron"
SEIR$VG[which(SEIR$VG=="original&alpha")]="Pre-Delta"
SEIR$VG = factor(SEIR$VG,ordered=TRUE,levels = c("Pre-Delta","Delta","Omicron"))
SEIR$Label_q = paste("Group",SEIR$Label_q)
library(reshape)
library(ggsignif)
library(ggprism)
colnames(SEIR)[2] = "NRMSE"
SEIR$r2[which(SEIR$r2<0)] = -SEIR$r2[which(SEIR$r2<0)]
val = melt(SEIR[,c("VG","NRMSE","r2")],id=c("VG"))
g2 = ggplot(val,aes(x=VG,y=value))+
  stat_boxplot(aes(color=VG),geom="errorbar",width=0.2,size=0.5)+
  geom_boxplot(outlier.shape=NA,show.legend = T,outlier.size=0.1,aes(color=VG),
               shape=2,outlier.stroke = 0.1, outlier.alpha = 45,
               notch = F,notchwidth = 0.2,
               fill="grey95",width=0.6,cex=0.5)+
  geom_jitter(aes(color=VG),width=.2,size=0.2)+ 
  theme_classic()+
  facet_wrap( ~ variable, scales = "free_y")+
  scale_color_manual(values = c("#899ebd","#427996",
                                "#003153"))+
  theme_prism(base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif")+
  theme(legend.position = "",
        panel.grid.major = element_line(colour = "transparent"),
        panel.grid.minor = element_line(colour = "transparent"),
        legend.title=element_text(color="black",size = unit(9, "pt")),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(9, "pt")),
        axis.title.x= element_blank(),
        axis.text.x = element_text(color="black",hjust=0.5,size = unit(9, "pt")),
        axis.title.y= element_blank(),
        axis.text.y= element_text(color="black",hjust=1,size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
ggsave("Figures_Reviewed/figSEIR_valb.pdf",g2, width=15, height=8, units="cm")
  
library(tidyverse)
val%>%group_by(variable,VG)%>%summarise(mean(value))












