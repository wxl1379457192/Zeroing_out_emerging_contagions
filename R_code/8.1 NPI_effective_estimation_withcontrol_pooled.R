library(rstan)
library(bayesplot)
library(dplyr)
rm(list = ls())
setwd("E:/zeroCOVID_NPI/Version0504")
source('code/NPI_code.R')

indir = "dataset"
##########################################################################
##################Individual NPI##########################
outdir = "Reviewed_IndividualNPI_stanV5_pooled_onebeta"
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
#mcon=stan_model('stanmodel/stanmodel_V5_pooled_variant.stan')
mcon=stan_model('stanmodel/stanmodel_V5_pooled_variant_onebeta.stan')
all[,"Temp"] = (all[,"Temp"]-min(all[,"Temp"]))/(max(all[,"Temp"])-min(all[,"Temp"]))
all$VG = factor(all$VG,levels = c("original&alpha","delta","omicron"))



al = subset(all,all$variant=="alpha")
or =  subset(all,all$variant=="original")
alw = length(unique(paste0(al$country,"_",al$original_start,"_",al$wave)))
orw = length(unique(paste0(or$country,"_",or$original_start,"_",or$wave)))
r01 = 3.32*orw/(alw+orw)+4.28*alw/(alw+orw)

r0 = as.vector(c(r01,4.9,9.5))
output = paste0(outdir, "/eachvariant.rds")
v = all[complete.cases(all),]
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

write.csv(all,paste0(outdir,"/testdata_IndNPI.csv"),row.names = F)
#####################Overall effect##################

outdir = "Reviewed_overallNPI_stanV5_overall"
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
mcon=stan_model('stanmodel/stanmodel_V5.stan')
outfile = data.frame()
all[,"Temp"] = (all[,"Temp"]-min(all[,"Temp"]))/(max(all[,"Temp"])-min(all[,"Temp"]))
k = unique(all[,c("city","original_start","variant")])
orw = length(which(k$variant=="original"))
alw = length(which(k$variant=="alpha"))
dew = length(which(k$variant=="delta"|k$variant=="delta&omicron"))
omw = length(which(k$variant=="omicron"))
r0 = 3.32*orw/nrow(k)+4.28*alw/nrow(k)+4.9*dew/nrow(k)+9.5*omw/nrow(k)
output = paste0(outdir, "/overall.rds")

v = all[complete.cases(all),]
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

write.csv(v,paste0(outdir,"/testdata_IndNPI.csv"),row.names = F)
