library("rstan")
library("bayesplot")
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
source('code/NPI_code.R')

indir = "1220dataset"
##########################################################################
##################Individual NPI##########################
outdir = "1220Version_Bayesian_model_validation"
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
auxname = c("Temp","fully_vaccination_rate")

mcon=stan_model('stanmodel/stanmodel_V3.stan')
outfile = data.frame()
all = subset(all,all$variant=="omicron")
ID = unique(paste0(all$citycode,"-",all$original_start))
all$ID = paste0(all$citycode,"-",all$original_start)

for(i in seq(1,length(ID))){
  c1 = all[which(all$ID!=ID[i]),]
  c1[,"Temp"] = (c1[,"Temp"]-min(c1[,"Temp"]))/(max(c1[,"Temp"])-min(c1[,"Temp"]))
  k = unique(c1[,c("city","original_start","variant")])
  
  r0=9.50
  output = paste0(outdir, "/validation_without_",ID[i],".rds")
  if (file.exists(output)==F){
    v = c1[complete.cases(c1),]
    X1 = as.matrix(v[,NPIname])
    X2 = as.matrix(v[,auxname])
    Rt = as.vector(v$Mean.R.)
    
    dataset<-list(m=length(Rt), k1=ncol(X1), k2=ncol(X2),
                  X1=X1, X2=X2 , Rt=Rt, r0=as.numeric(r0))  
    rstan_options(auto_write = TRUE)
    set.seed(0419)
    
    fit<-rstan::sampling(object = mcon,data=dataset,iter=5000,warmup=2500,
                         chains=5,thin=1,control = list(adapt_delta = 0.95, max_treedepth =15))
    saveRDS(fit, file = output)
    write.csv(c1,paste0(outdir,"/validation_dataset_",ID[i],".csv"),row.names = F)
  }else{
    print("This file has been exists")}
  
}


