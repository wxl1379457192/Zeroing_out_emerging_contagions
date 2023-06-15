library("rstan")
library("bayesplot")
library(dplyr)
setwd("E:/zeroCOVID_NPI/Version0504")
source('code/NPI_code.R')

indir = "dataset"
##########################################################################
##################Individual NPI##########################
outdir = "Reviewed_Bayesian_model_validation"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}


all = read.csv(paste0(indir,"/Rt&smooth_NPI.csv"),stringsAsFactors = F)
all = all[which(all$Mean.R.>0),]
all$Control[which(all$Control>1)] = 1
NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")
auxname = c("Temp","fully_vaccination_rate","pop_density")



mcon=stan_model('stanmodel/stanmodel_V5_pooled_variant_onebeta.stan')

outfile = data.frame()
ID = unique(paste0(all$citycode,"-",all$original_start))
all$ID = paste0(all$citycode,"-",all$original_start)

for(i in seq(1,length(ID))){
  c1 = all[which(all$ID!=ID[i]),]
  c1[,"Temp"] = (c1[,"Temp"]-min(c1[,"Temp"]))/(max(c1[,"Temp"])-min(c1[,"Temp"]))
  c1$VG = factor(c1$VG,levels = c("original&alpha","delta","omicron"))
  
  
  al = subset(c1,c1$variant=="alpha")
  or =  subset(c1,c1$variant=="original")
  alw = length(unique(paste0(al$country,"_",al$original_start,"_",al$wave)))
  orw = length(unique(paste0(or$country,"_",or$original_start,"_",or$wave)))
  r01 = 3.32*orw/(alw+orw)+4.28*alw/(alw+orw)
  
  r0 = as.vector(c(r01,4.9,9.5))
  output = paste0(outdir, "/validation_without_",ID[i],".rds")
  if (file.exists(output)==F){
    v = c1[complete.cases(c1),]
    v = do.call(rbind,lapply(split(v,v$VG),function(k){
      k = do.call(rbind,lapply(split(k,k$city),function(k1){
        k1 = arrange(k1,Date)
        return(k1)
      }))
      return(k)
    }))
    orw = nrow(v[which(v$VG =="original&alpha"),])
    dew = nrow(v[which(v$VG=="delta"),])
    omw = nrow(v[which(v$VG=="omicron"),])
    
    X1 = as.matrix(v[,NPIname])
    X2 = as.matrix(v[,auxname])
    Rt = as.vector(v$Mean.R.)
    Rt_sd = as.vector(v$Std.R.)
    dataset<-list(m = 3,lor = orw, lde = orw+dew, lom = nrow(v),
                  n =length(Rt), k1=ncol(X1), k2=ncol(X2),
                  X1=X1, X2=X2 , Rt=Rt, Rt_sd = Rt_sd, r0=r0)  
    rstan_options(auto_write = TRUE)
    set.seed(0419)
    
    fit<-rstan::sampling(object = mcon,data=dataset,iter=5000,warmup=2500,
                         chains=5,thin=1,control = list(adapt_delta = 0.95, max_treedepth =15))
    saveRDS(fit, file = output)
    
    write.csv(c1,paste0(outdir,"/validation_dataset_",ID[i],".csv"),row.names = F)
  }else{
    print("This file has been exists")}
  
}


