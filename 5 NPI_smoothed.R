setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
source('code/main_code.R')
#####################################
outdir = "1220dataset"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}
#############################################
indir="1025dataset"
Rt.NPI = read.csv(paste0(indir,"/NPI&cases_smooth&Rt&env&Vac.csv"), stringsAsFactors = F)

group = read.csv(paste0(indir,"/alldata_group_Qstatistics.csv"),stringsAsFactors = F)
Cases.control = read.csv(paste0("dataset/New_cases_all.csv"),stringsAsFactors = F)
#处理under。control数据缺失
Control.data = data.frame()
#Cases.control = do.call(rbind,lapply(split(Cases.control,Cases.control$city),function(c){
for(c in split(Cases.control,Cases.control$city)){
  c = c[,c(1:10)]
  c$Date = as.Date(c$Date)
  cg = subset(group,group$city==unique(c$city))
  if(nrow(cg)>0){
    cg$transmission_start = as.Date(cg$transmission_start)
    cg$end = as.Date(cg$end)
    call = do.call(rbind,lapply(seq(1,nrow(cg)),function(j){
      cgw = cg[j,]
      c1 = subset(c,c$Date>=cgw$transmission_start&c$Date<=cgw$end)
      f1 = data.frame(Date = seq.Date(from = cgw$transmission_start, to=cgw$end, by = "day"))
      c1 = merge(f1,c1,id="Date",all.x=TRUE)
      c1$city = unique(c1$city[which(is.na(c1$city)==FALSE)])
      c1$province = unique(c1$province[which(is.na(c1$province)==FALSE)])
      c1$ID = unique(c1$ID[which(is.na(c1$ID)==FALSE)])
      c1$New_cases[is.na(c1$New_cases)] = 0
      c1$New_asymptomatic[is.na(c1$New_asymptomatic)] = 0
      c1$New_confirmed[is.na(c1$New_confirmed)] = 0
      c1$Confirmed_from_asymptomatic[is.na(c1$Confirmed_from_asymptomatic)] = 0
      
      c1$Cases_under_control = as.numeric(c1$Cases_under_control)
      c1$确诊病例 = as.numeric(c1$确诊病例)
      c1$Control = c1$Cases_under_control/c1$New_cases
      
      if(length(is.na(c1$Cases_under_control))>0){
        k = which(c1$Control==1)[1]
        #if(is.na(k)==F){
        #  if(length(which(is.na(c1$Control[k:nrow(c1)]))==T)>0){
        #    c1$Control[which(is.na(c1$Control[k:nrow(c1)])==T)] = 1
        #  }
        #}
        c1$Control[is.na(c1$Control)] = c1$确诊病例[is.na(c1$Control)]/c1$New_confirmed[is.na(c1$Control)]
        c1$Control[which(is.na(c1$Control)&c1$Date<as.Date(cgw$plateau_start))] = 0
        c1$Control[which(is.na(c1$Control)&c1$Date>=as.Date(cgw$decline_start))] =0.8
        
        if(length(is.na(c1$Cases_under_control))>0){
          if(length(c1$Control[which(c1$Date>=as.Date(cgw$plateau_start)&
                                     c1$Date<as.Date(cgw$decline_start)&is.na(c1$Control)==F)])>0){
            max = c1$Control[which(c1$Date>=as.Date(cgw$plateau_start)&
                                     c1$Date<as.Date(cgw$decline_start)&is.na(c1$Control)==F)][1]
            if(max>0){
              c1$Control[which(is.na(c1$Control)&
                                 c1$Date>=as.Date(cgw$plateau_start)&
                                 c1$Date<as.Date(cgw$decline_start))] = seq(0,max,
                                                                            max/length(c1$Control[which(is.na(c1$Control)&
                                                                                                          c1$Date>=as.Date(cgw$plateau_start)&
                                                                                                          c1$Date<as.Date(cgw$decline_start))]))
            }else{
              c1$Control[which(is.na(c1$Control)&
                                 c1$Date>=as.Date(cgw$plateau_start)&
                                 c1$Date<as.Date(cgw$decline_start))] = 0 
            }
          }else{
            c1$Control[which(is.na(c1$Control)&
                               c1$Date>=as.Date(cgw$plateau_start)&
                               c1$Date<as.Date(cgw$decline_start))] = seq(0,1,
                                                                          1/length(c1$Control[which(is.na(c1$Control)&
                                                                                                      c1$Date>=as.Date(cgw$plateau_start)&
                                                                                                      c1$Date<as.Date(cgw$decline_start))]))
          }
        }
        
      }
      return(c1)
    }))
    #return(call)
    smo = as.numeric(stats::filter(call$Control,rep(1/5,5),sides=1))
    call$Control=smo
    if(length(is.na(call$Control))>0){
      call$Control[is.na(call$Control)]=0
    }
    call$Control[which(call$Control>1)]=1
    Control.data = rbind(Control.data,call)
  }
  #else{return(NULL)}
}
 
#}))

#outdir = "0916Version_withoutLOG_Label"
#outdir = "0916Version_withoutLOG_Label_VTran"

Cases.control = Control.data[,c(1,4,11)]
colnames(Cases.control)[2] = "citycode"
Cases.control$Date = as.Date(Cases.control$Date)
Rt.NPI$Date = as.Date(Rt.NPI$Date)
Rt.NPI = merge(Rt.NPI,Cases.control,by=c("citycode","Date"),Rt.NPI.x=TRUE)
Rt.NPI = do.call(rbind,lapply(split(Rt.NPI,Rt.NPI$citycode),function(k){
  l = length(k$Control[is.na(k$Control)])
  if(l == nrow(k)){
    k = NULL
  }else{
    k$Control[is.na(k$Control)] = 0
  }
  return(k)
}))
#write.csv(Rt.NPI,paste0(outdir,"/Rt&smooth_NPI_withControl.csv"),row.names = F)
Rt.NPI = unique(Rt.NPI)
NPInames = c("Lockdown","Facial_Mask",
             "Business_Premises_Closure","Public_Transportation_Closure","Gathering_restriction",
             "Workplace_Closure","School_Closure","Logistics_Management","Medicine_Management",
             "Mass_screening")


all.cor = do.call(rbind,lapply(split(Rt.NPI,Rt.NPI$VG),function(v){
  cor.variant = do.call(rbind,lapply(NPInames,function(n){
    cortest = do.call(rbind,lapply(seq(0,7),function(s){
      smooth.npi = smooth_NPI(n,s,v)
      merge = merge(v,smooth.npi,by=c("Date","city"))
      merge = merge[which(merge$Mean.R.!=0),]
      merge$NPI =  merge$NPI+rnorm(nrow(merge),sd=1e-6)
      core = cor.test(x=merge$Mean.R.,y=merge$NPI,method = "spearman",conf.level = 0.95)
      cor.result = data.frame(NPIname = n, step = s, variant = unique(v$VG),
                              rho = round(as.numeric(core$estimate),2), p.value = as.numeric(core$p.value))
      print(paste(unique(v$variant),"-",n,"-",s,"has been processed"))
      return(cor.result)
    }))
    return(cortest)
  }))
  return(cor.variant)
}))


write.csv(all.cor,paste0(outdir,"/Spearman rank correlation.csv"),row.names = F)
library(dplyr)
best = all.cor%>%group_by(variant,NPIname)%>%summarise(Best.step = step[which(rho==min(rho[is.na(rho)==F]))[1]],
                                                       rho = rho[which(rho==min(rho[is.na(rho)==F]))[1]])
#best$Best.step[is.na(best$Best.step)] = 7
write.csv(best,paste0(outdir,"/Spearman rank correlation_best_step.csv"),row.names = F)

all.cor = do.call(rbind,lapply(split(best,best$variant),function(k){
  NPI = subset(Rt.NPI,Rt.NPI$VG==unique(k$variant))
  n = NPInames[1]
  s = k$Best.step[which(k$NPIname==n)]
  new.npi = smooth_NPI(n,s,NPI)
  colnames(new.npi)[1] = n
  for(n in NPInames[2:length(NPInames)]){
    s = k$Best.step[which(k$NPIname==n)]
    smooth.npi = smooth_NPI(n,s,NPI)
    colnames(smooth.npi)[1] = n
    smooth.npi = unique(smooth.npi)
    new.npi = merge(smooth.npi,new.npi,by=c("Date","city","original_start","wave"))
    new.npi = unique(new.npi)
    print(n)
  }
  NPI = NPI[,-c(4:22)]
  newdata = merge(NPI,new.npi,by=c("Date","city","original_start","wave"))
  print(unique(k$variant))
  return(newdata)
}))
write.csv(all.cor,paste0(outdir,"/Rt&smooth_NPI.csv"),row.names = F)

#########################绘制NPI变化热力图，检查数据及平滑过程是否有误##########################
library(reshape2)
library(cowplot)
library(ggplot2)
pdf(paste0(outdir,"/plot_NPI_1025.pdf"))
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
