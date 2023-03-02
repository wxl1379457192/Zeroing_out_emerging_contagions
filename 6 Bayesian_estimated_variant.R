library("rstan")
library("bayesplot")
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
source('code/NPI_code.R')

indir = "1220dataset"
##########################################################################
##################Individual NPI##########################
outdir = "1220Version_IndividualNPI"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}
labelname = "Label_q"

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
for(c1 in split(all,all$VG)){
  c1[,"Temp"] = (c1[,"Temp"]-min(c1[,"Temp"]))/(max(c1[,"Temp"])-min(c1[,"Temp"]))
  testdata = domcmcV3(mcon, c1, NPIname, auxname, outdir)
  outfile  = rbind(outfile,testdata)
}

write.csv(outfile,paste0(outdir,"/testdata_IndNPI.csv"),row.names = F)
##########################################################################
outdir = "1220Version_IndividualNPI_overall"
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
