setwd("E:/zeroCOVID_NPI/Version0504")
library(data.table) 
library(EpiEstim) # see https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
library(ggplot2)
library(ggthemes)
library(dplyr)
source('code/main_code.R')
Sys.setlocale("LC_TIME","English")

dir="1025dataset"
##################################################
p.all = read.csv(paste0(dir,"/NPIdataset_processed.csv"),stringsAsFactors = F)
city = split(p.all,p.all$citycode)
all=data.frame()
for(k in 1:length(city)){
  c= city[[k]]
  clist = unique(c$original_start)
  
  f = do.call(rbind,lapply(clist,function(j){
    cw = subset(c,c$original_start==j)
    cw$Date = as.Date(cw$Date)
    casedata_dates = data.frame(dates= seq.Date(from = min(cw$Date)-7, to=max(cw$Date)+6, by = "day"))
    casedata_dates$I = 0
    for (k in seq(1,nrow(cw))){
      d = cw$Date[k]
      casedata_dates$I[which(casedata_dates$dates==d)] = cw$New_cases[which(cw$Date==d)]
    }
    casedata_dates = casedata_dates[order(casedata_dates$dates),]
    casedata_dates$I = smooth_Gaussian(casedata_dates$I) 
    
    cw$New_cases_smooth = casedata_dates$I[c(which(casedata_dates$dates==min(cw$Date)):which(casedata_dates$dates==max(cw$Date)))]
    cw$growth_rate = 0
    for(i in seq(2,nrow(cw))){
      if(cw$New_cases_smooth[i-1]!=0){
        cw$growth_rate[i] = cw$New_cases_smooth[i]/cw$New_cases_smooth[i-1]
      }else{
        cw$growth_rate[i] = 0
      }
      cw$growth_rate = round(cw$growth_rate,2)
    }
    
    return(cw)
  }))
  all = rbind(all,f)
}
#write.csv(all,paste0(dir,"/NPI&cases_smooth.csv"),row.names = F)
###########################################################################
out.ID = read.csv("dataset/alldata_group_Qstatistics.csv", stringsAsFactors = F)
#casesISO = read.csv(paste0(dir,"/alldata&NPI_smoothed.csv"),stringsAsFactors = F)
casesISO = all

or = 5.8  
sda = 3.2

al = 3.18
sda = 4.36

del = 2.3
sdd = 3.4

om = 2.9
sdo = 1.6

outdir = dir
#if (dir.exists(outdir )==F){dir.create(outdir)
#}else{print("This path has been exists")}

alldata = data.frame()

#pdf(paste0(outdir,"/",unique(case$³ÇÊĞ),"ID",pandc$outbreak_ID,"-",f,".pdf"))
pdf(paste0(outdir,"/plot_all_countriesV2.pdf"))

for (i in 1:nrow(out.ID)){
  pand = out.ID[i,]
  cases = subset(casesISO,casesISO$citycode==pand$citycode)
  cases$Date=as.Date(cases$Date)
  pand$Outbreak_date = as.Date(pand$Outbreak_date)
  cases = cases[which(cases$original_start==pand$Outbreak_date),]
  cases = cases[which(cases$Date>=as.Date(pand$transmission_start)&cases$Date<=as.Date(pand$end)),]

  T <- nrow(cases)+5
  startday = which(cases$New_cases_smooth>0)[1]
  t_start <- seq(startday+2, T-6) # starting at 3 as conditional on the past observations
  t_end <- t_start + 6 # adding 4 to get 6-day windows as bounds included in window
  cases$I = cases$New_cases_smooth
  estimatedI = c(cases$I,rep(0,5))
  if(pand$variant == "omicron"){
    res_parametric_si <- estimate_R(estimatedI,
                                    method="parametric_si",
                                    config = make_config(list(
                                      t_start = t_start,
                                      t_end = t_end,
                                      mean_si = om, # the mean of serial interval - please change it according to covid Variant ref:https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taac037/6545354
                                      std_si = sdo))) # the sd of serial interval - please change it according to covid Variant
    
    }else if(pand$variant == "alpha"){
     res_parametric_si <- estimate_R(estimatedI, 
                                     method="parametric_si",
                                     config = make_config(list(
                                       t_start = t_start,
                                       t_end = t_end,
                                       mean_si = al, # the mean of serial interval - please change it according to covid Variant ref:https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taac037/6545354
                                       std_si = sda))) # the sd of serial interval - please change it according to covid Variant
   
     }else if(pand$variant == "delta"){
     res_parametric_si <- estimate_R(estimatedI, 
                                     method="parametric_si",
                                     config = make_config(list(
                                       t_start = t_start,
                                       t_end = t_end,
                                       mean_si = del, # the mean of serial interval - please change it according to covid Variant
                                       std_si = sdd))) # the sd of serial interval - please change it according to covid Variant
    }else if(pand$variant == "delta&omicron"){
      res_parametric_si <- estimate_R(estimatedI, 
                                      method="parametric_si",
                                      config = make_config(list(
                                        t_start = t_start,
                                        t_end = t_end,
                                        mean_si = (del+om)/2, # the mean of serial interval - please change it according to covid Variant
                                        std_si = (sdd+sdo)/2))) # the sd of serial interval - please change it according to covid Variant
    }else{
      res_parametric_si <- estimate_R(estimatedI,
                                      method="parametric_si",
                                      config = make_config(list(
                                        t_start = t_start,
                                        t_end = t_end,
                                        mean_si = or, # the mean of serial interval - please change it according to covid Variant
                                        std_si = sda)))
    }
  
  daily_R_data = data.frame(res_parametric_si$R)
  daily_R_data = data.frame(Date = res_parametric_si$dates[res_parametric_si$R$t_start]+cases$Date[which(cases$New_cases_smooth>0)[1]]
                            ,daily_R_data[,3:11])
  daily_R_data = daily_R_data[order(daily_R_data$Date),]
  
  daily_R_data$variant = pand$variant
  #daily_R_data$label = pand$label
  #daily_R_data$label_Tran = pand$label_Tran
  #daily_R_data$label_Env = pand$label_Env
  daily_R_data$GDP = pand$GDP
  daily_R_data$pop_density = pand$pop_density
  daily_R_data$pop = pand$pop
  daily_R_data$Doctor.number = pand$Doctor.number
  daily_R_data$wave = pand$wave
  daily_R_data$VG = pand$VG
  daily_R_data$Label_q = pand$Group
  daily_R_data = merge(cases,daily_R_data,by="Date",all.y=TRUE)
  s1 = which(daily_R_data$Quantile.0.025.R.>1)[1]
  e1 = which(daily_R_data$Quantile.0.975.R.<1)[length(which(daily_R_data$Quantile.0.975.R.<1))]
  if(is.na(s1)==F){
    if(length(e1)>0){
      daily_R_data = daily_R_data[c(s1:e1),]
    }else{
      daily_R_data = daily_R_data[c(s1:nrow(daily_R_data)),]
    }
  }
  plotdata = daily_R_data[which(daily_R_data$Mean.R.>0),]
  daily_R_data[is.na(daily_R_data)] = 0
  alldata = rbind(alldata,daily_R_data)
  g1<-ggplot(plotdata)+geom_col(aes(x=Date,y=I),color="grey80",alpha=0.9)+
    geom_vline(aes(xintercept = as.Date(pand$Outbreak_date)),linetype ="dotdash", color= "black")+
    labs(x=NULL,y="Incidence",title=paste0("Epidemic curve of cityID: ",pand$citycode,
                                           " Pandemic start:",pand$Outbreak_date))+
    theme(legend.position = "right",
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x = element_text(color="black",size = 10),
          axis.title.y= element_text(color="black",size = 10),
          axis.text.y = element_text(color="black",size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",color="black"))
  
  ##plot estimate Rt
  g2<-ggplot(plotdata)+
    geom_vline(aes(xintercept = as.Date(pand$Outbreak_date)),linetype ="dotdash", color= "black")+
    geom_line(aes(x=Date,y=Mean.R.,color="Mean"))+
    geom_ribbon(aes(x=Date, ymin = Quantile.0.025.R., ymax = Quantile.0.975.R.,fill="95% CI"),
                alpha = 0.3)+
    labs(x=NULL,y="Rt",title="Estimated Rt")+
    scale_colour_manual("",values="#14a1a7")+
    scale_fill_manual("",values="#c6e7e6")+
    geom_hline(yintercept = 1,linetype ="dotdash",size=1)+
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x = element_text(color="black",size = 10),
          axis.title.y= element_text(color="black",size = 10),
          axis.text.y = element_text(color="black",size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",color="black"))
  
  gridExtra::grid.arrange(grobs = list(g1, g2))
  print(paste(i,pand$citycode," -pan:",pand$Outbreak_date,"has been processed"))
  
            
}

dev.off()
alldata=unique(alldata)
write.csv(alldata,paste0(dir,"/NPI&cases_smooth&Rt.csv"),row.names = F)
