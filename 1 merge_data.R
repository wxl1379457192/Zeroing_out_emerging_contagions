setwd("D:/ÂÛÎÄÍ¶¸å/Zero_COVID/20230131/")
library(data.table) 
library(EpiEstim) # see https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
library(ggplot2)
library(ggthemes)
library(dplyr)
source('code/main_code.R')
Sys.setlocale("LC_TIME","English")
casesISO = read.csv("dataset/New_cases_all.csv",stringsAsFactors = F)                                                                     
ID = unique(casesISO[,c("city","province")])
#write.csv(casesISO,"dataset/New_cases_all.csv",row.names = F)
#write.csv(ID,"dataset/New_cases_ID.csv",row.names = F)
out.ID = read.csv("dataset/Outbreak_information_all.csv",stringsAsFactors = F)[,c(1:12)]
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
p.all = read.csv(paste0(dir,"/NPI_all.csv"),stringsAsFactors = F)

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

