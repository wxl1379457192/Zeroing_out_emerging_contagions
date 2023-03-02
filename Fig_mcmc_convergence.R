library(ggplot2)
library(rstan)
library("bayesplot")
library("gridExtra")
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
path = "1220Version_IndividualNPI"
outpath = "FiguresV1129"
efdata<-list()
list<-list.files(dir,pattern="*.rds")
conresult<-do.call(rbind,lapply(list,function(l){
  f<-readRDS(paste0(path,"/",l))
  rhat<-as.data.frame(rhat(f))
  Rneff<-as.data.frame(neff_ratio(f))
  dd<-do.call(cbind,list(rhat,Rneff))
  colnames(dd)<-c("rhat","R_ESS")
  return(dd)
}))
valplot<-list()
valplot[[1]]<-ggplot(conresult,aes(rhat))+
  geom_histogram(bins=100,colour="white",fill="#2D6587",size=0.01)+
  labs(x=NULL,y=NULL,title="a")+
  scale_x_continuous(expand=c(0,0),limits=c(0.995,1.008))+
  scale_y_continuous(expand=c(0,0))+geom_vline(xintercept = 1,linetype ="dotdash",size=1)+
  theme(axis.line.y= element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size = unit(9, "pt")),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title=element_text(size = unit(9, "pt"),hjust=0),
        panel.grid=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        panel.background=element_rect(fill = "transparent",colour = NA),
        plot.margin=unit(c(0.1,0.2,0.1,0.2),"cm"),
        axis.text=element_text(size = unit(9, "pt")))
valplot[[2]]<-ggplot(conresult,aes(R_ESS))+
  geom_histogram(bins=100,colour="white",fill="#2D6587",size=0.01)+
  labs(x=NULL,y=NULL,title="b")+scale_x_continuous(expand=c(0.02,0.02),limits=c(0,1.2))+
  geom_vline(xintercept = 1,linetype ="dotdash",size=1)+
  scale_y_continuous(expand=(c(0,0)))+
  theme(axis.line.y= element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.x =element_text(size = unit(9, "pt")),
        plot.title=element_text(size =unit(9, "pt"),hjust=0),
        panel.grid=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        panel.background=element_rect(fill = "transparent",colour = NA),
        plot.margin=unit(c(0.1,0.2,0.1,0.2),"cm"),
        axis.text=element_text(size = unit(9, "pt")))
plot<-grid.arrange(arrangeGrob(grobs =valplot,ncol =2))
ggsave(paste0(outpath,"/mcmc_convergence.pdf"),plot,units="mm",width=150,height=40,device = cairo_pdf)
############################################################################################
