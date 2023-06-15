setwd("E:/zeroCOVID_NPI/Version0504")
library(data.table) 
library(EpiEstim) # see https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
library(ggplot2)
library(ggthemes)
library(dplyr)
library(EpiNow2)
source('code/main_code.R')
Sys.setlocale("LC_TIME","English")
#cases = merge_cases_data(dir="E:/zeroCOVID_NPI/Rt_0709/codeV0803/dataset/")
#casesISO = do.call(rbind,lapply(split(cases,cases$city),function(c){
#  d1 = do.call(rbind,lapply(split(c,c$Date),function(d){
#    if(nrow(d)>1){
#      d = d[which(d$Confirmed_from_asymptomatic>0),]
      #d = d[which(d$New_cases==d$New_asymptomatic+d$New_confirmed-d$Confirmed_from_asymptomatic),]
      #d = unique(d)
#      print(paste(d$city,d$Date))
#    }
#    return(d)
#  }))
#  return(d1)
#}))
casesISO = read.csv("dataset/New_cases_all.csv",stringsAsFactors = F)                                                                     
ID = unique(casesISO[,c("city","province")])
#write.csv(casesISO,"dataset/New_cases_all.csv",row.names = F)
#write.csv(ID,"dataset/New_cases_ID.csv",row.names = F)
out.ID = read.csv("dataset/Outbreak_information.csv",stringsAsFactors = F)[,c(1:12)]


# #################Rt estimation based on EpiNow2 package###########
# out.ID = read.csv("dataset/alldata_group_Qstatistics.csv", stringsAsFactors = F)
# Rt_estimation = data.frame()
# for(i in unique(out.ID$city)){
#   U = out.ID[which(out.ID$city==i),]
#   cases.city = casesISO[which(casesISO$city==i),]
#   cases.city$Date = as.Date(cases.city$Date)
#   U$Outbreak_date = as.Date(U$Outbreak_date)
#   U$transmission_start = as.Date(U$transmission_start)
#   U$end = as.Date(U$end)
#   U$label = 0
#   U = arrange(U, Outbreak_date)
#   if(nrow(U)>1){for(x in seq(2,nrow(U))){
#     if(U$transmission_start[x] == U$transmission_start[x-1]){
#       U$label[x-1] = 1
#     }
#   }}
#   U = U[which(U$label==0),]
#   for (x in seq(1,nrow(U))){
#     ID = U[x,] 
#     if (ID$end-ID$transmission_start>3){
#       start = ifelse(ID$transmission_start>=ID$Outbreak_date,ID$transmission_start,ID$Outbreak_date)
#       case = cases.city[which(cases.city$Date>=start&
#                               cases.city$Date<=ID$end+7),]
#       case$Date = as.Date(case$Date)
#       case = arrange(case, Date)
#       if (nrow(case)>0&ID$variant!=""){
#         case <- case[,c('Date', 'New_cases')]
#         casesISO2 <- data.frame(Date = seq(min(case$Date)-7,
#                                            max(case$Date),by = '1 day'))
#         casesISO2  = merge(case,casesISO2,by = "Date",all =T)
#         casesISO2$New_cases[is.na(casesISO2$New_cases)] = 0
#         
#         if(ID$variant == "original"){
#           generation_time <- generation_time_opts(mean = 5.7, sd = 3.8)
#           incubation_period <- data.table(mean = 1.78,sd = 0.52,
#                                           mean_sd = 0.1,sd_sd =0.1,
#                                           max= 10,
#                                           fixed = TRUE,dist = "lognormal")
#         
#         }else if(ID$variant == "alpha"){
#           generation_time <- generation_time_opts(mean = 4.7, sd = 3.3)
#           incubation_period <- data.table(mean = 1.50,sd = 0.46,
#                                           mean_sd = 0.1,sd_sd =0.1, max= 10,
#                                           fixed = TRUE,dist = "lognormal")
#         }else if(ID$variant == "delta"){
#           generation_time <- generation_time_opts(mean = 2.9, sd = 3.0)
#           incubation_period <- data.table(mean = 1.25,sd = 0.34,
#                                           mean_sd = 0.1,sd_sd =0.1, max= 10, 
#                                           fixed = TRUE,dist = "lognormal")
#         }else if(ID$variant == "omicron"){
#           generation_time <- generation_time_opts(mean = 2.36, sd = 0.59)
#           incubation_period <- data.table(mean = 1.02,sd = 0.45,
#                                           mean_sd = 0.1,sd_sd =0.1, max= 10,
#                                           fixed = TRUE,dist = "lognormal")
#         }else{
#           generation_time <- generation_time_opts(mean = (2.9+2.36)/2, sd = (3.0+0.59)/2)
#           incubation_period <- data.table(mean = (1.25+1.02)/2,sd = (0.34+0.45)/2,
#                                           mean_sd = 0.1,sd_sd =0.1, max= 10,
#                                           fixed = TRUE,dist = "lognormal")
#         }
#         reporting_delay <- list(mean = convert_to_logmean(2, 1),mean_sd = 0.1,
#                                 sd = convert_to_logsd(2, 1),sd_sd = 0.1,max = 10)
#         colnames(casesISO2) = c("date","confirm")
#         estimates <-  estimate_infections(
#           reported_cases = casesISO2,
#           generation_time = generation_time,
#           delays = delay_opts(incubation_period, reporting_delay),
#           rt = rt_opts(prior = list(mean = 2, sd = 0.5)),
#           backcalc = backcalc_opts(),
#           stan =  stan_opts(control = list(adapt_delta = 0.99)),
#           # obs = obs_opts(scale = list(mean = 0.8, sd = 0.05)),
#           horizon = 0
#         )
#         estimate_R = summary(estimates, type = "parameters", params = "R")
#         
#         casesISO2 <- merge(casesISO2, estimate_R, by= 'date', all=T)
#         
#         casesISO2 = casesISO2[!is.na(casesISO2$variable)]
#         casesISO2 = casesISO2[,c(1,2,6:14)]
#         casesISO2$city = i
#         casesISO2$citycode = ID$citycode
#         casesISO2$original_start = ID$Outbreak_date
#         casesISO2$wave = ID$wave
#         Rt_estimation = rbind(Rt_estimation,casesISO2)
#       }
#     }
#     print(paste(i,ID$start,"has successed"))
#   }
# }
# write.csv(Rt_estimation,"dataset/Rt_estimated_VEpiNow2_V2.csv",row.names = F)
#################Rt estimation based on EpiEstim package###########
repeat.cases = do.call(rbind,lapply(seq(1,50),function(j){
  cases.delay = do.call(rbind,lapply(unique(out.ID$city),function(i){
    cases.city = casesISO[which(casesISO$city==i),]
    cases.city$Date = as.Date(cases.city$Date)
    U  = out.ID[which(out.ID$city==i),]
    U$start = as.Date(U$start)
    U$end = as.Date(U$end)
    f = do.call(rbind,lapply(seq(1,nrow(U)),function(x){
      ID = U[x,] 
      if (ID$end-ID$start>3){
        case = cases.city[which(cases.city$Date>=ID$start&cases.city$Date<=ID$end),]
        mass.start = case$Date[which(case$Mass_screening>=2)][1]
        case = case[which(case$New_cases!=0),]
        if (nrow(case)>0&ID$variant!=""){
          case.time <- case[,c('Date', 'New_cases')]
          case.time$New_cases[is.na(case.time$New_cases)] = 0
          case.time <- case.time[rep(seq(nrow(case.time)), case.time$New_cases),]
          
          if(ID$variant == "original"){
            #case.time$exposepd <- round(rlnorm(nrow(case.time), 1.788, 0.52),0)#https://www.nature.com/articles/s41598-021-91834-8/tables/1
            #case.time$exposepd <- round(rlnorm(nrow(case.time), log(5.2), log(1.6)),0)
            case.time$exposepd <- round(rlnorm(nrow(case.time), 1.67, 0.28),0)
          }else if(ID$variant == "alpha"){
            case.time$exposepd <- round(rlnorm(nrow(case.time), 1.15, 0.41),0)
          }else if(ID$variant == "delta"){
            case.time$exposepd <- round(rlnorm(nrow(case.time), 0.91, 0.42),0)
            
          }else if(ID$variant == "omicron"){
            case.time$exposepd <- round(rlnorm(nrow(case.time), 0.72, 0.48),0)
          }else{
            case.time$exposepd <- round(rlnorm(nrow(case.time), (0.91+0.72)/2, (0.42+0.48)/2),0)
          }
          case.time$expose_report <- case.time$Date - case.time$exposepd 
          if(min(case.time$Date)<=ID$start){
            before.mass = nrow(case.time[which(case.time$Date <= ID$start),])
            case.time$onset_report <- 0 
            case.time$onset_report[which(case.time$Date <= ID$start)] <- round(rlnorm(before.mass, 0.82, 0.84),0)#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9225637/
            case.time$onset_report[which(case.time$Date > ID$start)] <- round(rbinom(nrow(case.time)-before.mass, 1, 0.3),0)
          }else{
            case.time$onset_report[which(case.time$Date > ID$start)] <- round(rbinom(nrow(case.time), 1, 0.3),0)
          }
          case.time$expose_report <- case.time$expose_report - case.time$onset_report 
          #case.time$expose_report <- case.time$Date
          case.time$New_cases <- 1
          case.time <- data.table(case.time)
          case.time <- case.time[, list(New_cases = sum(New_cases)), by='expose_report']
          casesISO2 <- data.frame(Date = seq(min(case.time$expose_report),
                                             max(case.time$expose_report),by = '1 day'))
          casesISO2 <- merge(casesISO2, case.time, by.x= 'Date', by.y='expose_report', all.x=T)
          casesISO2$New_cases[is.na(casesISO2$New_cases)] <- 0
          casesISO2$city = i
          casesISO2$citycode = ID$citycode
          casesISO2$original_start = ID$start
          return(casesISO2)
        }else{return(NULL)}
      }else(return(NULL))
    }))
    return(f)
  }))
  dd = cases.delay%>%group_by(Date,city,citycode,original_start)%>%summarise(New_cases=sum(New_cases))
  dd$WID = j
  return(dd)
  #return(cases.delay)
}))
cases.delay = repeat.cases%>%group_by(Date,city,
                                      citycode,original_start)%>%summarise(New_cases=round(sum(New_cases)/50,0))
                                                             #Total_cases= sum(New_cases))
cases.delay = unique(cases.delay)
################plot NPI################
library(reshape2)
library(cowplot)
pdf(paste0("plot_NPI.pdf"))
NPInames = c("Lockdown","Facial_Mask",
             "Business_Premises_Closure","Public_Transportation_Closure","Gathering_restriction",
             "Workplace_Closure","School_Closure","Logistics_Management","Medicine_Management",
             "Mass_screening")
dir = "dataset"
p.all = read.csv(paste0(dir,"/NPI_dataset.csv"),stringsAsFactors = F)

NPI.all = data.frame()
for(p in split(p.all,p.all$citycode)){
  colnames(p)[1] = "Date"
  p$Date=as.Date(p$Date)
  p.delay = subset(cases.delay,cases.delay$city==unique(p$city))
  for(pd in split(p.delay,p.delay$original_start)){
    if(sum(pd$New_cases)>=60){
      pd$Date = as.Date(pd$Date)
      pall = merge(p,pd[,c(1,4,5)],id="Date",all.y=T)
      pall = pall[order(pall$Date),]
      if(length(pall$Lockdown[is.na(pall$Lockdown)])!=nrow(pall)){
        pall[is.na(pall$Lockdown),c(2:22)] = pall[which(is.na(pall$Lockdown)==F)[1],c(2:22)]
        pall$New_cases[is.na(pall$New_cases)] = 0
        k1 = pall[,c(NPInames,"Date")]
        k1 = melt(k1,id="Date")
        p1 = ggplot(k1)+ 
          geom_tile(aes(x=Date,y=variable,fill=value))+
          labs(title = paste(unique(pall$citycode[which(is.na(pall$citycode)==F)]),unique(pall$original_start)))+
          theme(legend.position = "bottom",legend.title=element_blank(),
                plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
                axis.title.y= element_blank(),
                axis.text.y = element_text(color="black",size = 10),
                plot.margin=unit(c(2,0.1,2,0.1),"cm"),
                panel.background=element_rect(fill = "transparent",color="black"))
        print(p1)
        NPI.all = rbind(NPI.all,pall)
        print(paste(unique(pall$citycode[which(is.na(pall$citycode)==F)]),unique(pall$original_start)))
      }
      
    }
  }
}
dev.off()

write.csv(NPI.all,"dataset/NPIdataset.csv",row.names = F)

