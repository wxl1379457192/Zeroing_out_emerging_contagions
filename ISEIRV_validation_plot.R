library(ggplot2)
library(openxlsx)
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
SEIR =  read.xlsx("SEIR_1220_validation/Validation_SEIR.xlsx", sheet = 1)
SEIR$VG[which(SEIR$VG==1)]="Delta"
SEIR$VG[which(SEIR$VG==2)]="Omicron"
SEIR$VG[which(SEIR$VG==3)]="Pre-Delta"
SEIR$VG = factor(SEIR$VG,ordered=TRUE,levels = c("Pre-Delta","Delta","Omicron"))
SEIR$Label_q = paste("Group",SEIR$Label_q)
model.lm<-lm(formula = Real_cases ~ Predicted_cases, data = SEIR)
summary(model.lm)
l <- list(r2 = format(summary(model.lm)$r.squared, digits = 4),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))

eq <- substitute(italic(R)^2~"="~r2~","~italic(P-value)~"="~p, l)


g = ggplot(data=SEIR)+geom_point(aes(x=Predicted_cases,y=Real_cases,shape = as.factor(Label_q), 
                                 color = as.factor(VG),fill = as.factor(VG)),alpha=0.4,size=1.5)+
  coord_trans(x = "log",y="log")+scale_y_continuous(breaks=c(250,2500,25000,250000))+
  scale_x_continuous(breaks=c(250,2500,25000,250000))+
  scale_color_manual(values=c("#500105","grey50","#0b5c9e"),name = "Variant",
                     guide = guide_legend(ncol = 2, byrow = TRUE))+
  scale_fill_manual(values=c("#500105","grey50","#0b5c9e"),name = "Variant",
                    guide = guide_legend(ncol = 2, byrow = TRUE))+
  scale_shape_manual(values=c(20,17,22,21),name = "Group",
                     guide = guide_legend(ncol = 2, byrow = TRUE))+
  geom_text(aes(x = 600, y = 500000, label = as.character(as.expression(eq))), parse = TRUE,size = unit(3, "pt"))+
  theme_test()+labs(x="Simulated cases",y="Observed cases")+
  #guides(col = )+
  theme(legend.position = "right",
        panel.grid.major = element_line(colour = "transparent"),
        panel.grid.minor = element_line(colour = "transparent"),
        legend.title=element_text(color="black",size = unit(9, "pt")),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(9, "pt")),
        axis.title.x= element_text(color="black",size = unit(9, "pt")),
        axis.text.x = element_text(color="black",hjust=0.5,size = unit(9, "pt")),
        axis.title.y= element_text(color="black",size = unit(9, "pt")),
        axis.text.y= element_text(color="black",hjust=1,size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
outdir = "Figures"
ggsave(paste0(outdir,"/Figure_ISEIRV_validation_point.pdf"),g, width=15, height=9, units="cm")
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
  geom_jitter(aes(color=VG),width=.2,size=0.5)+ 
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
ggsave(paste0(outdir,"/Figure_ISEIRV_validation_boxplot.pdf"),g2, width=15, height=8, units="cm")
  
library(tidyverse)
val%>%group_by(variable,VG)%>%summarise(mean(value))












