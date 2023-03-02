
domcmcV2 <- function(mcon, all, NPIname, vacname, outdir,labelname ="label"){
  i = unique(all$VG)
  #r0: orginal = 3.32 From: Estimate of the Basic Reproduction Number for COVID-19: A Systematic Review and Meta-analysis
  #alpha= 4.02
  #delta: 4.9 From: Transmission dynamics and epidemiological characteristics of SARS-CoV-2 Delta variant infections in Guangdong, China, May to June 2021
  #omericn = 9.5
  if (i=="omicron"){
    r0=9.50
    output = paste0(outdir, "/", i,"_Group",unique(all[,labelname]),"_all.rds")
  }else if(i=="delta"){
    r0=4.9
    output = paste0(outdir, "/", i,"_Group",unique(all[,labelname]),"_all.rds")
  }else{
    al = subset(all,all$variant=="alpha")
    or =  subset(all,all$variant=="original")
    alw = length(unique(paste0(al$country,"_",al$original_start,"_",al$wave)))
    orw = length(unique(paste0(or$country,"_",or$original_start,"_",or$wave)))
    r0 = 3.32*orw/(alw+orw)+4.28*alw/(alw+orw)
    output = paste0(outdir, "/original&alpha_Group",unique(all[,labelname]),"_all.rds")
  }
  
  v = all[complete.cases(all),]
  X1 = as.matrix(v[,NPIname])
  X2 = as.matrix(v[,vacname])
  Rt = as.vector(v$Mean.R.)
  
  dataset<-list(m=length(Rt), k1=ncol(X1), k2=ncol(X2),
                X1=X1, X2=X2 , Rt=Rt, r0=as.numeric(r0))  
  rstan_options(auto_write = TRUE)
  set.seed(0419)
  
  fit<-rstan::sampling(object = mcon,data=dataset,iter=5000,warmup=2500,
                       chains=5,thin=1,control = list(adapt_delta = 0.95, max_treedepth =15))
  saveRDS(fit, file = output)
  
  return(v)
}
domcmcV3 <- function(mcon, all, NPIname, vacname, outdir){
  i = unique(all$VG)
  #r0: orginal = 3.32 From: Estimate of the Basic Reproduction Number for COVID-19: A Systematic Review and Meta-analysis
  #alpha= 4.02
  #delta: 4.9 From: Transmission dynamics and epidemiological characteristics of SARS-CoV-2 Delta variant infections in Guangdong, China, May to June 2021
  #omericn = 9.5
  if (i=="omicron"){
    r0=9.50
    output = paste0(outdir, "/", i,"_all.rds")
  }else if(i=="delta"){
    r0=4.9
    output = paste0(outdir, "/", i,"_all.rds")
  }else{
    al = subset(all,all$variant=="alpha")
    or =  subset(all,all$variant=="original")
    alw = length(unique(paste0(al$country,"_",al$original_start,"_",al$wave)))
    orw = length(unique(paste0(or$country,"_",or$original_start,"_",or$wave)))
    r0 = 3.32*orw/(alw+orw)+4.28*alw/(alw+orw)
    output = paste0(outdir, "/original&alpha","_all.rds")
  }
  
  v = all[complete.cases(all),]
  X1 = as.matrix(v[,NPIname])
  X2 = as.matrix(v[,vacname])
  Rt = as.vector(v$Mean.R.)
  
  dataset<-list(m=length(Rt), k1=ncol(X1), k2=ncol(X2),
                X1=X1, X2=X2 , Rt=Rt, r0=as.numeric(r0))  
  rstan_options(auto_write = TRUE)
  set.seed(0419)
  
  fit<-rstan::sampling(object = mcon, data=dataset, iter=5000, warmup=2500,
                       chains=5, thin=1,control = list(adapt_delta = 0.95, max_treedepth =15))
  saveRDS(fit, file = output)
  
  return(v)
}
domcmc_group <- function(mcon, all, NPIname, vacname, outdir,labelname ="label"){
  #r0: orginal = 3.32 From: Estimate of the Basic Reproduction Number for COVID-19: A Systematic Review and Meta-analysis
  #alpha= 4.02
  #delta: 4.9 From: Transmission dynamics and epidemiological characteristics of SARS-CoV-2 Delta variant infections in Guangdong, China, May to June 2021
  #omericn = 9.5
  al = subset(all,all$variant=="alpha")
  or =  subset(all,all$variant=="original")
  om = subset(all,all$variant == "omicron")
  de =  subset(all,all$variant == "delta")
  alw = length(unique(paste0(al$country,"_",al$original_start,"_",al$wave)))
  orw = length(unique(paste0(or$country,"_",or$original_start,"_",or$wave)))
  omw = length(unique(paste0(om$country,"_",om$original_start,"_",om$wave)))
  dew = length(unique(paste0(de$country,"_",de$original_start,"_",de$wave)))
  
  r0 = 3.32*orw/(alw+orw+omw+dew)+4.28*alw/(alw+orw+omw+dew)+
    9.50*omw/(alw+orw+omw+dew)+4.90*dew/(alw+orw+omw+dew)
  output = paste0(outdir, "/Group",unique(all[,labelname]),"_all.rds")
  
  v = all[complete.cases(all),]
  X1 = as.matrix(v[,NPIname])
  X2 = as.matrix(v[,vacname])
  Rt = as.vector(v$Mean.R.)
  
  dataset<-list(m=length(Rt), k1=ncol(X1), k2=ncol(X2),
                X1=X1, X2=X2 , Rt=Rt, r0=as.numeric(r0))  
  rstan_options(auto_write = TRUE)
  set.seed(0419)
  
  fit<-rstan::sampling(object = mcon,data=dataset,iter=5000,warmup=2500,
                       chains=5,thin=1,control = list(adapt_delta = 0.95, max_treedepth =15))
  saveRDS(fit, file = output)
  
  return(v)
}
getresult = function(dir,all,
                     NPIname = c("Static_management","Gathering_restriction",
                                 "Logistics_Management","Medicine_Management",
                                 "Mass_screening","Facial_Mask"),
                     vacname = c("Practical.vaccination"),labelname="label"
){
  conlist = list.files(dir,pattern="*.rds")
  testdata=all[complete.cases(all),]
  #testdata = read.csv(paste0(dir, "/testdata.csv"),stringsAsFactors = F)
  result  = do.call(rbind,lapply(seq(1,length(conlist)),function(j){
    conf = conlist[j]
    f = readRDS(paste0(dir,"/",conf))
    variant = strsplit(conf,"_")[[1]][1]
    label =  substring(strsplit(conf,"_")[[1]][2],6,7)
    
    #if(variant=="original&alpha"){
    #  v = testdata[which(testdata$variant=="original"|testdata$variant=="alpha"),]
    #}else{v = testdata[which(testdata$variant==variant),]}
    v = testdata[which(testdata$VG==variant),]
    v = v[which(v[,labelname]==label),]
    #for(k in vacname){
    #  v[,k] = v[,k]/max(v[,k])
    #}
    
    out = rstan::extract(f)
    al=as.data.frame(out$alpha)
    colnames(al)=NPIname
    for (i in NPIname){
      al[,i] = (1-exp(-al[,i]*mean(v[,i])))*100
      #if(max(v[,i])!=mean(v[,i])){
      #  al[,i] = (1-exp(-al[,i]))*100
      #}else{al[,i] = 0}
    }
    pard=mcmc_intervals_data(al,prob = .5,prob_outer= .95,point_est="mean")
    pard$variant = variant
    pard$label = label
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
          be[,j]=(1-exp(-be[,j]*mean(vac[,j])))*100
        }else{
          be[,j]=0
        }
      }
    }else{
      if(max(vac)>0){
        be[,vacname]=(1-exp(-be[,vacname]*mean(vac)))*100
      }else{
        be[,vacname]=0
      }
    }
    
    parv=mcmc_intervals_data(be,prob = .5,prob_outer= .95,point_est="mean")
    parv$variant = variant
    parv$label =label
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

getpre = function(dir,all,
                     NPIname = c("Static_management","Gathering_restriction",
                                 "Logistics_Management","Medicine_Management",
                                 "Mass_screening","Facial_Mask"),
                     vacname = c("Practical.vaccination"),labelname="label"
){
  conlist = list.files(dir,pattern="*.rds")
  testdata=all[complete.cases(all),]
  #testdata = read.csv(paste0(dir, "/testdata.csv"),stringsAsFactors = F)
  result  = do.call(rbind,lapply(seq(1,length(conlist)),function(j){
    conf = conlist[j]
    f = readRDS(paste0(dir,"/",conf))
    variant = strsplit(conf,"_")[[1]][1]
    label =  substring(strsplit(conf,"_")[[1]][2],6,7)
    
    #if(variant=="original&alpha"){
    #  v = testdata[which(testdata$variant=="original"|testdata$variant=="alpha"),]
    #}else{v = testdata[which(testdata$variant==variant),]}
    v = testdata[which(testdata$VG==variant),]
    v = v[which(v[,labelname]==label),]
    for(k in vacname){
      v[,k] = v[,k]/max(v[,k])
    }
    
    out = rstan::extract(f)
    
    y = as.data.frame(out$ypre)
  
    colnames(y)= seq(1,ncol(y))
   
    pard=mcmc_intervals_data(y,prob = .5,prob_outer= .95,point_est="mean")
    
    v$ypre = pard$m
    v$ypre_CI5 = pard$ll
    v$ypre_CI25 = pard$l
    v$ypre_CI75 = pard$h
    v$ypre_CI95 = pard$hh
    
    return(v)
  }))
  return(result)
}