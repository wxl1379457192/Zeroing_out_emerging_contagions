setwd("E:/zeroCOVID_NPI/Version0504")
source('code/main_code.R')
#########################绘制NPI变化热力图，检查数据及平滑过程是否有误##########################
library(reshape2)
library(cowplot)
library(ggplot2)

outdir = "dataset"
all.cor = read.csv(paste0(outdir,"/alldata&NPI_smoothed.csv"),stringsAsFactors = F)


pdf(paste0(outdir,"/plot_NPI.pdf"))
NPInames = c("Lockdown","Facial_Mask",
             "Business_Premises_Closure","Public_Transportation_Closure","Gathering_restriction",
             "Workplace_Closure","School_Closure","Logistics_Management","Medicine_Management",
             "Mass_screening")

for(p in split(all.cor,all.cor$citycode)){
  p$Date=as.Date(p$Date)
  for(pc in split(p,p$original_start)){
    for(pd in split(pc,pc$wave)){
      k1 = pd[,c(NPInames,"Date")]
      k1 = melt(k1,id="Date")
      p1 = ggplot(k1)+ 
        geom_tile(aes(x=Date,y=variable,fill=value))+
        labs(title = paste(unique(pd$citycode[which(is.na(pd$citycode)==F)]),unique(pd$original_start),unique(pd$wave)))+
        theme(legend.position = "bottom",legend.title=element_blank(),
              plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
              axis.title.y= element_blank(),
              axis.text.y = element_text(color="black",size = 10),
              plot.margin=unit(c(2,0.1,2,0.1),"cm"),
              panel.background=element_rect(fill = "transparent",color="black"))
      print(p1)
      print(paste(unique(pd$citycode[which(is.na(pd$citycode)==F)]),unique(pd$original_start)))
    }
  }
}
dev.off()

#################################################################################################
#for(i in NPInames){
#  all.cor[,i] =  all.cor[,i]/max(all.cor[,i])
#  all.cor[,i][which(all.cor[,i]<0)] = 0
#}

N = all.cor[,c(NPInames,"pop_density","Temp","Hum")]
N = N[complete.cases(N),]

library("ggcor")

COR<-function(N){
  p1<-quickcor(N, cor.test = TRUE)+geom_star(data = get_data(type = "upper" ,show.diag = FALSE))+
    geom_mark(data = get_data(type = "lower", show.diag = FALSE),size=2.5)+
    geom_abline(size=0.5,slope=-1,intercept =ncol(N)+1)+
    scale_fill_gradient2(midpoint = 0, low = "#00498d", mid = "white", high = "#74230a",space="Lab")+
    theme(legend.position = "right",
          legend.text = element_text(color="black",size=unit(9, "pt")),
          legend.title = element_text(color="black",size=unit(9, "pt")),
          legend.key.height=unit(10,'mm'),
          legend.key.width=unit(2,'mm'),
          axis.text = element_text(color="black",size=unit(9, "pt")))
  return(p1)
  #ggsave("picture/cor.pdf",p1,width=200,height=200,units="mm",device = cairo_pdf)
}
COR(N)


