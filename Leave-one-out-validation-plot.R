#sensitivity analysis : Bayesian#
library("rstan")
library("bayesplot")
library("ggalt")
library("ggplot2")
library("scales")
library("data.table")
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
Sys.setlocale("LC_TIME","English")
R2<-function(x,y){
  xm<-mean(x)
  ssres<-sum((x-xm)^2)
  ssreg<-sum((y-x)^2)
  return(ssreg/(ssres+ssreg))
}
CI95lowfun<-function(data){
  if(sd(data)>0){
    k<-t.test(data)$conf.int[1]
  }else{k<-0}
  return(k)
}
CI95highfun<-function(data){
  if(sd(data)>0){
    k<-t.test(data)$conf.int[2]
  }else{k<-0}
  return(k)
}

alldir = "1220Version_Bayesian_model_validation"
conlist = list.files(alldir,pattern="*.rds")
NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")
ObNPI = c("Static.management",
          "Facial_Mask",
          "PCR.screening",
          "Contact.control") 
auxname = c("Temp","fully_vaccination_rate")

valdir = "1214Version_IndividualNPI_overall"
pd = read.csv(paste0(valdir,"/testdata_IndNPI.csv"),stringsAsFactors = F)
pd[,"Temp"] = (pd[,"Temp"]-min(pd[,"Temp"]))/(max(pd[,"Temp"])-min(pd[,"Temp"]))
outpath = alldir
par = do.call(rbind,lapply(conlist,function(c){
  citycode = substring(strsplit(c,"_")[[1]][3],1,6)
  ordate = substring(strsplit(c,"_")[[1]][3],8,17)
  data = pd[which(pd$citycode==citycode),]
  data = subset(data,data$original_start==ordate)
  data$Date=as.Date(data$Date)
  data = data[order(data$Date),]
  fit = readRDS(paste0(alldir,"/",c))

  out = rstan::extract(fit)
  r0=mean(out$R0)
  
  al=as.data.frame(out$alpha)
  colnames(al)=NPIname
  for (i in NPIname){
    data[,i]= -quantile(al[,i])[3]*data[,i]
  }
  be=as.data.frame(out$beta)
  colnames(be)=auxname
  for (j in auxname){
    data[,j]= -quantile(be[,j])[3]*data[,j]
  }
  al = apply(data[,NPIname],1,sum)
  be = apply(data[,auxname],1,sum)
  Rt = exp(al+be)*r0
  data$pre = Rt
  
  #CI05
  data1 = pd[which(pd$citycode==citycode),]
  data1 = subset(data1,data1$original_start==ordate)
  data1$Date=as.Date(data1$Date)
  data1 = data1[order(data1$Date),]
  al=as.data.frame(out$alpha)
  colnames(al)=NPIname
  for (i in NPIname){
    data1[,i]= -quantile(al[,i])[2]*data1[,i]
  }
  be=as.data.frame(out$beta)
  colnames(be)=auxname
  for (j in auxname){
    data1[,j]= -quantile(be[,j])[2]*data1[,j]
  }
  al = apply(data1[,NPIname],1,sum)
  be = apply(data1[,auxname],1,sum)
  Rt = exp(al+be)*r0
  data$y_CI05 = Rt
  
  #CI95
  data2 = pd[which(pd$citycode==citycode),]
  data2 = subset(data2,data2$original_start==ordate)
  data2$Date=as.Date(data2$Date)
  data2 = data2[order(data1$Date),]
  al=as.data.frame(out$alpha)
  colnames(al)=NPIname
  for (i in NPIname){
    data2[,i]= -quantile(al[,i])[4]*data2[,i]
  }
  be=as.data.frame(out$beta)
  colnames(be)=auxname
  for (j in auxname){
    data2[,j]= -quantile(be[,j])[4]*data2[,j]
  }
  al = apply(data2[,NPIname],1,sum)
  be = apply(data2[,auxname],1,sum)
  Rt = exp(al+be)*r0
  data$y_CI95 = Rt
  meltr = melt(data[,c("Date","Mean.R.","pre")],id="Date")
  rtfig<-ggplot()+
    geom_ribbon(data=data,aes(x=as.Date(Date),ymin=Quantile.0.05.R.,ymax=Quantile.0.95.R.),fill="#ae7d74",alpha=0.3)+
    geom_ribbon(data=data,aes(x=as.Date(Date),ymin=y_CI05,ymax=y_CI95),fill="#8491B4B2",alpha=0.3)+
    geom_point(data=meltr,aes(x=as.Date(Date),y=value,color=variable),alpha=0.5,size=0.3)+
    geom_xspline(data=meltr,aes(x=as.Date(Date),y=value,color=variable),size=0.2)+
    scale_x_date(labels=date_format("%b %d\n%Y"),expand=c(0.01,0.01))+
    scale_color_manual(values=c("#dc0000cc","#365a87"))+
    labs(y =expression(R[t]),x = NULL)+
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(10, "pt")),
          axis.title.x= element_blank(),
          axis.text.x = element_text(color="black",hjust=0.5,size = unit(9,"pt")),
          axis.title.y= element_text(color="black",size = unit(9,"pt")),
          axis.text.y= element_text(color="black",hjust=1,size = unit(9,"pt")),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",colour ="black"))
  ggsave(paste0(outpath,"/",gsub(".rds",".pdf",c)),rtfig,units="mm",width=70,height=70,device = cairo_pdf)
  return(data)
}))

pofig = ggplot(par,aes(x=Mean.R.,y=pre))+
  geom_point(alpha=0.3,size=0.6,color="#008080")+
  labs(y =expression("Prediction of "*R[t]),x = expression(R[t]))+
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(10, "pt")),
        axis.text.x = element_text(color="black",hjust=0.5,size = unit(10,"pt")),
        axis.title.y= element_text(color="black",vjust=0,size = unit(10,"pt")),
        axis.title.x= element_text(color="black",vjust=0,size = unit(10,"pt")),
        axis.text.y= element_text(color="black",hjust=1,size = unit(10,"pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = "black"))
ggsave(paste0(outpath,"/ALL_Point.pdf"),pofig,units="mm",width=140,height=140,device = cairo_pdf)
as.data.frame(t(caret::postResample(par$Mean.R., par$pre)))
val = do.call(rbind,lapply(split(par,par$citycode),function(dc){
  k = do.call(rbind,lapply(split(dc,dc$original_start),function(data){
    validation<-as.data.frame(t(caret::postResample(data$Mean.R., data$pre)))
    validation$R2N<-R2(data$Mean.R., data$pre)
    validation$or = unique(data$original_start)
    return(validation)
  }))
  k$citycode = unique(dc$citycode)
  return(k)
}))

write.csv(val,paste0(outpath,"/validation_all.csv"),row.names = F)