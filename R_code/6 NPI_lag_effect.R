setwd("E:/zeroCOVID_NPI/Version0504")
rm(list = ls())
source('code/main_code.R')
#####################################
outdir = "dataset"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}
#############################################
indir="pre-dataset"
Rt.NPI = read.csv(paste0(indir,"/NPI&cases_smooth&Rt&env&Vac.csv"), stringsAsFactors = F)

group = read.csv(paste0(indir,"/alldata_group_Qstatistics.csv"),stringsAsFactors = F)
Cases.control = read.csv(paste0("dataset/New_cases_all.csv"),stringsAsFactors = F)


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
      c1$ȷ�ﲡ�� = as.numeric(c1$ȷ�ﲡ��)
      c1$Control = c1$Cases_under_control/c1$New_cases
      
      if(length(is.na(c1$Cases_under_control))>0){
        k = which(c1$Control==1)[1]
        #if(is.na(k)==F){
        #  if(length(which(is.na(c1$Control[k:nrow(c1)]))==T)>0){
        #    c1$Control[which(is.na(c1$Control[k:nrow(c1)])==T)] = 1
        #  }
        #}
        c1$Control[is.na(c1$Control)] = c1$ȷ�ﲡ��[is.na(c1$Control)]/c1$New_confirmed[is.na(c1$Control)]
        c1$Control[which(is.na(c1$Control)&c1$Date<as.Date(cgw$plateau_start))] = 0
        c1$Control[which(is.na(c1$Control)&c1$Date>=as.Date(cgw$decline_start))] = 0.8
        
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

Cases.control = Control.data[,c(1,4,11)]
colnames(Cases.control)[2] = "citycode"
Cases.control$Date = as.Date(Cases.control$Date)
Rt.NPI$Date = as.Date(Rt.NPI$Date)
Rt.NPI = merge(Rt.NPI,Cases.control,by=c("citycode","Date"),Rt.NPI.x=TRUE)
Rt.NPI = do.call(rbind,lapply(split(Rt.NPI,Rt.NPI$citycode),function(k){
  k = unique(k)
  l = length(k$Control[is.na(k$Control)])
  if(l == nrow(k)){
    k = NULL
  }else{
    k$Control[is.na(k$Control)] = 0
  }
  return(k)
}))
Rt.NPI = do.call(rbind,lapply(split(Rt.NPI,Rt.NPI$city),function(c){
  c = do.call(rbind,lapply(split(c,c$original_start),function(co){
    co = do.call(rbind,lapply(split(co,co$wave),function(cow){
      datelist = unique(cow$Date)
      cow = do.call(rbind,lapply(datelist,function(d){
        k=subset(cow,cow$Date==d)
        if(nrow(k)>1){
          k = k[order(k$Lockdown),]
          k = k[1,]
        }
        return(k)
      }))
      return(cow)
    }))
    return(co)
  }))
  return(c)
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

