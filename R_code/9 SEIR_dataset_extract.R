library("rstan")
library("bayesplot")
setwd("E:/zeroCOVID_NPI/Version0504")
source('code/NPI_code.R')

indir = "Reviewed_IndividualNPI_stanV5_pooled_onebeta"
ordir = "dataset"
dir = "SEIR"
if (dir.exists(dir)==F){dir.create(dir)
}else{print("This path has been exists")}
cities = read.csv("dataset/New_cases_ID.csv",stringsAsFactors = F)

getresult = function(indir,outdir,NPIname,conname){
  conlist = list.files(indir,pattern="*.rds")
  all = read.csv(paste0(indir,"/testdata_IndNPI.csv"),stringsAsFactors = F)
  all = all[which(all$Mean.R.>0),]
  all = merge(all,cities,by="city")
  write.csv(all,paste0(outdir,"/dataset.csv"),row.names = F)
  f = readRDS(paste0(indir,"/",conlist))
  out = rstan::extract(f)
  al=as.data.frame(out$alpha)
  variant = c("original&alpha","delta","omicron")
  syn2 = do.call(rbind,lapply(seq(1,length(variant)),function(x){
    v = all[which(all$VG==variant[x]),]
    alv = al[,((x-1)*length(NPIname)+1):(x*length(NPIname))] 
    colnames(alv)=NPIname

    for (i in NPIname){
      alv[,i] = (1-exp(-alv[,i]))*100
    }
    
    synv=mcmc_intervals_data(alv,prob = .5,prob_outer= .95,point_est="mean")
    synv$strength = 0
    for (p in unique(synv$parameter)){
      synv$strength[which(synv$parameter==p)] = mean(v[,p]) 
    }
    synv$variant = unique(v$VG)
    
    return(synv)
  }))
  
  be = as.data.frame(out$beta)
  colnames(be) = conname
  for (i in conname){
    be[,i] = (1-exp(-be[,i]))*100
  }
  synb=mcmc_intervals_data(be,prob = .5,prob_outer= .95,point_est="mean")
  synb$strength = 0
  for (p in unique(synb$parameter)){
    synb$strength[which(synb$parameter==p)] = mean(all[,p]) 
  }
  synb$variant = "all"
  return(do.call(rbind,list(syn2,synb)))
}
NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")
auxname = c("Temp","fully_vaccination_rate","pop_density")
l = getresult(indir,dir,NPIname,auxname)
write.csv(l,paste0(dir,"/NPIefficacy.csv"),row.names = F)