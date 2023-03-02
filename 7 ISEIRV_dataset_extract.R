library("rstan")
library("bayesplot")
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
source('code/NPI_code.R')

indir = "1214Version_IndividualNPI"
ordir = "dataset"
dir = "SEIR_1214"
if (dir.exists(dir)==F){dir.create(dir)
}else{print("This path has been exists")}

getresult = function(indir,outdir,NPIname,vacname){
  conlist = list.files(indir,pattern="*.rds")
  testdata = read.csv(paste0(indir, "/testdata_IndNPI.csv"),stringsAsFactors = F)
  write.csv(testdata,paste0(outdir,"/dataset.csv"),row.names = F)
  result  = do.call(rbind,lapply(seq(1,length(conlist)),function(j){
    conf = conlist[j]
    f = readRDS(paste0(indir,"/",conf))
    variant = strsplit(conf,"_")[[1]][1]
    
    v = testdata[which(testdata$VG==variant),]
    out = rstan::extract(f)
    al=as.data.frame(out$alpha)
    colnames(al)=NPIname
    for (i in NPIname){
      al[,i] = (1-exp(-al[,i]))*100
      #if(max(v[,i])!=mean(v[,i])){
      #  al[,i] = (1-exp(-al[,i]))*100
      #}else{al[,i] = 0}
    }
    pard=mcmc_intervals_data(al,prob = .5,prob_outer= .95,point_est="mean")
    pard$variant = variant
    pard$strength = 0
    for (p in unique(pard$parameter)){
      pard$strength[which(pard$parameter==p)] = mean(v[,p]) 
    }
    
    be = as.data.frame(out$beta)
    colnames(be) = vacname 
    vac = v[,vacname]
    
    if (length(vacname)>1){
      for (j in vacname){
        if(max(vac[,j])>0){
          be[,j]=(1-exp(-be[,j]))*100
        }else{
          be[,j]=0
        }
      }
    }else{
      if(max(vac)>0){
        be[,vacname]=(1-exp(-be[,vacname]))*100
      }else{
        be[,vacname]=0
      }
    }
    
    parv=mcmc_intervals_data(be,prob = .5,prob_outer= .95,point_est="mean")
    parv$variant = variant
    #parv$strength = mean(vac)
    parv$strength = 0
    for (p in unique(parv$parameter)){
      parv$strength[which(parv$parameter==p)] = mean(v[,p]) 
    }
    #r0 = as.data.frame(out$R0)
    #colnames(r0) = "R0"
    #r0 = mcmc_intervals_data(r0,prob = .5,prob_outer= .95,point_est="mean")
    #r0$variant = variant
    par = do.call(rbind,list(pard,parv))
    return(par)
  }))
  return(result)
}
NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")

auxname = c("Temp","fully_vaccination_rate")
l = getresult(indir,dir,NPIname,auxname)
write.csv(l,paste0(dir,"/NPIefficacy.csv"),row.names = F)