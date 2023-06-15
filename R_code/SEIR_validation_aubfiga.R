library(ggplot2)
library(openxlsx)
library(readxl)
setwd("E:/zeroCOVID_NPI/Version0504")

dir = "SEIR_validation/SEIR_simulation"
filelist  = list.files(dir, pattern=paste0("*.xlsx"))
valnum  = do.call(rbind,lapply(filelist,function(file){
  f = read.xlsx(file.path(dir,file))
  valcum = data.frame(VG = unique(f$VG), citycode = unique(f$citycode),
                      original_start = unique(f$original_start),
                      Observed_cases = sum(f$Observed_cases[(nrow(f)-unique(f$validation_length)):nrow(f)]),
                      Predicted_cases = sum(f$Predicted_cases[(nrow(f)-unique(f$validation_length)):nrow(f)]),
                      Label_q = unique(f$Label_q))
  #f = f[(nrow(f)-unique(f$validation_length)):nrow(f),]
  print(paste(unique(f$citycode),unique(f$original_start)))
  return(valcum)
}))

valnum = valnum[which(valnum$Observed_cases<=1000),]
model.lm<-lm(formula = Predicted_cases~Observed_cases, data = valnum)

summary(model.lm)
l <- list(r2 = format(summary(model.lm)$r.squared, digits = 2),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))
eq <- substitute(italic(R)^2~"="~r2~","~italic(P-value)~"="~p, l)

g = ggplot(valnum)+geom_point(aes(x=Observed_cases,y = Predicted_cases,
                          shape = as.factor(Label_q), 
                          color = as.factor(VG),fill = as.factor(VG)),alpha=0.4,size=1.5)+
  scale_x_continuous(breaks=seq(0,200,50),limits = c(0,200))+
  scale_y_continuous(breaks=seq(0,200,50),limits = c(0,200))+
  #scale_y_continuous(breaks=c(250,2500,25000,250000))+
  #scale_x_continuous(breaks=c(250,2500,25000,250000))+
  scale_color_manual(values=c("#500105","grey50","#0b5c9e"),name = "Variant",
                     guide = guide_legend(ncol = 2, byrow = TRUE))+
  scale_fill_manual(values=c("#500105","grey50","#0b5c9e"),name = "Variant",
                    guide = guide_legend(ncol = 2, byrow = TRUE))+
  scale_shape_manual(values=c(20,17,22,21),name = "Group",
                     guide = guide_legend(ncol = 2, byrow = TRUE))+
  #geom_text(aes(x = 600, y = 500000, label = as.character(as.expression(eq))), parse = TRUE,size = unit(3, "pt"))+
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
ggsave("Figures_Reviewed/figSEIR_vala_testdata.pdf",g, width=15, height=9, units="cm")

vale  = do.call(rbind,lapply(filelist,function(file){
  f = read.xlsx(file.path(dir,file))
  f = f[(nrow(f)-unique(f$validation_length)):nrow(f),]
  model.lm<-lm(formula = Predicted_cases~Observed_cases, data = f)
  valcum = data.frame(VG = unique(f$VG), citycode = unique(f$citycode),
                      original_start = unique(f$original_start),
                      r2 = format(summary(model.lm)$r.squared, digits = 2),
                      p = format(summary(model.lm)$coefficients[2,4], digits = 4),
                      label_q = unique(f$Label_q))
  
  print(paste(unique(f$citycode),unique(f$original_start)))
  return(valcum)
}))

SEIR%>%group_by(VG)%>%summarise(mean(r2))




