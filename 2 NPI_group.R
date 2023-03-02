setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
library(data.table) 
library(EpiEstim) # see https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
library(ggplot2)
library(ggthemes)
library(dplyr)
source('code/main_code.R')
Sys.setlocale("LC_TIME","English")
###########################
################plot NPI################
library(reshape2)
library(cowplot)

pdf(paste0("plot_NPI_1025.pdf"))
NPInames = c("Lockdown","Facial_Mask",
             "Business_Premises_Closure","Public_Transportation_Closure","Gathering_restriction",
             "Workplace_Closure","School_Closure","Logistics_Management","Medicine_Management",
             "Mass_screening")

dir = "dataset"
outdir = "1025dataset"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}

p.all = read.csv(paste0(dir,"/NPIdataset.csv"),stringsAsFactors = F)

p.all$Business_Premises_Closure = p.all$Business_Premises_Closure_intensity/4+p.all$Business_Premises_Closure-1
p.all$Public_Transportation_Closure = p.all$Public_Transportation_Closure_intensity/4+p.all$Public_Transportation_Closure-1
p.all$Gathering_restriction = p.all$Gathering_restriction_intensity/4+p.all$Gathering_restriction-1
p.all$Workplace_Closure = p.all$Workplace_Closure_intensity/4+p.all$Workplace_Closure-1
p.all$School_Closure = p.all$School_Closure_intensity/4+p.all$School_Closure-1
p.all$Logistics_Management = p.all$Logistics_Management_intensity/4+p.all$Logistics_Management-1
p.all$Medicine_Management = p.all$Medicine_Management_intensity/4+p.all$Medicine_Management-1
p.all$Mass_screening = p.all$Mass_screening*p.all$Mass_screening_intensity/7
for(i in NPInames){
  p.all[,i] =  p.all[,i]/max(p.all[,i][is.na(p.all[,i])==F])
  p.all[,i][which(p.all[,i]<0)] = 0
}

for(p in split(p.all,p.all$citycode)){
  colnames(p)[1] = "Date"
  p$Date=as.Date(p$Date)
  for(pd in split(p,p$original_start)){
    if(sum(pd$New_cases)>=60){
      pd$Date = as.Date(pd$Date)
      k1 = pd[,c(NPInames,"Date")]
      k1 = melt(k1,id="Date")
      p2 = ggplot(pd)+geom_line(aes(x=Date,y=New_cases))+
        labs(title = paste(unique(pd$citycode[which(is.na(pd$citycode)==F)]),unique(pd$original_start)))+
        theme(legend.position = "bottom",legend.title=element_blank(),
              plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
              axis.title.y= element_blank(),
              axis.title.x= element_blank(),
              axis.text.y = element_text(color="black",size = 10),
              plot.margin=unit(c(0.5,0.1,0.1,4.5),"cm"),
              panel.background=element_rect(fill = "transparent",color="black"))
      p1 = ggplot(k1)+ 
        geom_tile(aes(x=Date,y=variable,fill=value))+
        labs(title = paste(unique(pd$citycode[which(is.na(pd$citycode)==F)]),unique(pd$original_start)))+
        theme(legend.position = "bottom",legend.title=element_blank(),
              plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
              axis.title.y= element_blank(),
              axis.text.y = element_text(color="black",size = 10),
              plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
              panel.background=element_rect(fill = "transparent",color="black"))
   
      gridExtra::grid.arrange(grobs = list(p2, p1))
      print(paste(unique(pd$citycode[which(is.na(pd$citycode)==F)]),unique(pd$original_start)))
    }
  }
}
dev.off()

write.csv(p.all,paste0(outdir,"/NPIdataset_processed.csv"),row.names = F)
#################################
#dir = "dataset"
#p.all = read.csv(paste0(dir,"/NPIdataset_processed.csv"),stringsAsFactors = F)
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
    casedata_dates = unique(casedata_dates)
    for (g in seq(1,nrow(cw))){
      d = cw$Date[g]
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
   # cw$smooth = cw$New_cases_smooth
    #cw$smooth[which(cw$New_cases_smooth<=2)] = 0
    cw$diff_cases = c(0,diff(cw$New_cases_smooth))
    
    transmission = cw[which(cw$growth_rate>=1.3),]
    transmission$diff = c(100,diff(transmission$Date))
    transmission_start = transmission$Date[which(transmission$diff>1&transmission$diff_cases>0)]
    if(length(transmission_start)>1){
      for (t in 2:length(transmission_start)){
        if(transmission_start[t]-transmission_start[t-1]<=7){
          transmission_start[t] = transmission_start[t-1]
        }
      }
      transmission_start = unique(transmission_start)
    }

    plateau = cw[which(cw$growth_rate<1.3&cw$growth_rate>=0.9),]
    plateau$diff = c(100,diff(plateau$Date))
    plateau_start = plateau$Date[which(plateau$diff>1&plateau$diff_cases>0)]
    #if(length(plateau_start)>1){
    #  for (t in 2:length(plateau_start)){
    #    if(plateau_start[t]-plateau_start[t-1]<=7){
    #      plateau_start[t] = plateau_start[t-1]
    #    }
    #  }
    #  plateau_start = unique(plateau_start)
    #}
    
    decline = cw[which(cw$growth_rate<0.9&cw$growth_rate>0),]
    decline$diff = c(100,diff(decline$Date))
    decline_start = decline$Date[which(decline$diff>1&decline$diff_cases<0)]
    if(length(decline_start)>1){
      for (t in 2:length(decline_start)){
        if(decline_start[t]-decline_start[t-1]<=7){
          decline_start[t] = decline_start[t-1]
        }
      }
      decline_start = unique(decline_start)
    }
   
    end = cw[which(cw$diff_cases<0),]
    end$diff = c(diff(end$Date),100)
    end_date = end$Date[which(end$diff>2)]
    if(length(end_date)>1){
      for (t in 2:length(end_date)){
        if(end_date[t]-end_date[t-1]<=7){
          end_date[t] = end_date[t-1]
        }
      }
      end_date = unique(end_date)
    }
    
    k1 = data.frame()

    if(length(end_date)>=length(transmission_start)){
      if(length(transmission_start)>1){
        for(h in 1:(length(transmission_start)-1)){
          if(transmission_start[h+1]<end_date[h]){
            transmission_start[h+1] = NA
          }
        }
        transmission_start = transmission_start[which(is.na(transmission_start)==F)]
        if(max(transmission_start)>max(end_date)){
          transmission_start = transmission_start[which(transmission_start!=max(transmission_start))]
        }
      }
     
      for(i in 1:length(transmission_start)){
        k = data.frame(citycode = unique(cw$citycode),
                       city = unique(cw$city),
                       Outbreak_date = unique(cw$original_start),
                       wave = i,
                       transmission_start = as.Date(ifelse(i<length(transmission_start),
                                                   ifelse(transmission_start[i+1]>end_date[i],transmission_start[i],transmission_start[i-1])
                                                   ,transmission_start[i])),
                       plateau_start = plateau_start[which(plateau_start>transmission_start[i])][1],
                       decline_start = decline_start[which(decline_start>transmission_start[i])][1],
                       end = as.Date(ifelse(i<length(transmission_start),
                                            ifelse(end_date[i]<=transmission_start[i+1],
                                                   max(end_date[which(end_date<transmission_start[i+1]&end_date>transmission_start[i])]),
                                                   end_date[i+1]),
                                            max(end_date))))
        #k = k[complete.cases(k),]
        k$last = as.numeric(k$end-k$transmission_start)
        k$total_cases = sum(cw$New_cases[which(cw$Date>=k$transmission_start&cw$Date<=k$end)])
        k$max_cases = max(cw$New_cases[which(cw$Date>=k$transmission_start&cw$Date<=k$end)])
        cww = cw[which(cw$Date>=k$transmission_start&cw$Date<=k$end),]
        k$Outbreak_date = as.Date(k$Outbreak_date)
        
        cww$Lockdown.diff = c(0,diff(cww$Lockdown))
        k$Lockdown.ts = ifelse(max(cww$Lockdown.diff)>0,as.numeric(min(cww$Date[which(cww$Lockdown.diff>0)])-k$transmission_start),0)
      #k$Lockdown.ts = ifelse(max(cww$Lockdown)>=0.4,as.numeric(min(cww$Date[which(cww$Lockdown>0.4)])-k$transmission_start),0)
        k$Lockdown.os = ifelse(max(cww$Lockdown.diff)>0,as.numeric(min(cww$Date[which(cww$Lockdown.diff>0)])-k$Outbreak_date),0)
        k$Lockdown.intensity = ifelse(max(cww$Lockdown.diff)>0,cww$Lockdown[which(cww$Date==min(cww$Date[which(cww$Lockdown.diff>0)]))],max(cww$Lockdown))
        
        cww$MS.diff = c(0,diff(cww$Mass_screening))
        k$MS.ts = ifelse(max(cww$MS.diff)>0,as.numeric(min(cww$Date[which(cww$MS.diff>0)])-k$transmission_start),0)
        k$MS.os = ifelse(max(cww$MS.diff)>0,as.numeric(min(cww$Date[which(cww$MS.diff>0)])-k$Outbreak_date),0)
        #k$MS.ts = ifelse(max(cww$Mass_screening)>=4/35,as.numeric(min(cww$Date[which(cww$Mass_screening>=4/35)])-k$transmission_start),0)
        #k$MS.os = ifelse(max(cww$Mass_screening)>=4/35,as.numeric(min(cww$Date[which(cww$Mass_screening>=4/35)])-k$Outbreak_date),0)
        k$MS.intensity = ifelse(max(cww$MS.diff)>0,cww$Mass_screening[which(cww$Date==min(cww$Date[which(cww$MS.diff>0)]))],
                                max(cww$Mass_screening))
        k1 = rbind(k1,k)
      }
    }else{
      for(i in 1:length(end_date)){
        k = data.frame(citycode = unique(cw$citycode),
                       city = unique(cw$city),
                       Outbreak_date = unique(cw$original_start),
                       wave = i,
                       transmission_start = min(transmission_start[which(transmission_start<end_date[i])]), 
                       plateau_start = min(plateau_start[which(plateau_start<end_date[i])]),
                       decline_start = min(decline_start[which(decline_start<end_date[i])]),
                       end = end_date[i])
        k$last = as.numeric(k$end-k$transmission_start)
        k$total_cases = sum(cw$New_cases[which(cw$Date>=k$transmission_start&cw$Date<=k$end)])
        k$max_cases = max(cw$New_cases[which(cw$Date>=k$transmission_start&cw$Date<=k$end)])
        cww = cw[which(cw$Date>=k$transmission_start&cw$Date<=k$end),]
        k$Outbreak_date = as.Date(k$Outbreak_date)
        
        cww$Lockdown.diff = c(0,diff(cww$Lockdown))
        k$Lockdown.ts = ifelse(max(cww$Lockdown.diff)>0,as.numeric(min(cww$Date[which(cww$Lockdown.diff>0)])-k$transmission_start),0)
        #k$Lockdown.ts = ifelse(max(cww$Lockdown)>=0.4,as.numeric(min(cww$Date[which(cww$Lockdown>0.4)])-k$transmission_start),0)
        k$Lockdown.os = ifelse(max(cww$Lockdown.diff)>0,as.numeric(min(cww$Date[which(cww$Lockdown.diff>0)])-k$Outbreak_date),0)
        k$Lockdown.intensity = ifelse(max(cww$Lockdown.diff)>0,cww$Lockdown[which(cww$Date==min(cww$Date[which(cww$Lockdown.diff>0)]))],max(cww$Lockdown))
        
        cww$MS.diff = c(0,diff(cww$Mass_screening))
        k$MS.ts = ifelse(max(cww$MS.diff)>0,as.numeric(min(cww$Date[which(cww$MS.diff>0)])-k$transmission_start),0)
        k$MS.os = ifelse(max(cww$MS.diff)>0,as.numeric(min(cww$Date[which(cww$MS.diff>0)])-k$Outbreak_date),0)
        #k$MS.ts = ifelse(max(cww$Mass_screening)>=4/35,as.numeric(min(cww$Date[which(cww$Mass_screening>=4/35)])-k$transmission_start),0)
        #k$MS.os = ifelse(max(cww$Mass_screening)>=4/35,as.numeric(min(cww$Date[which(cww$Mass_screening>=4/35)])-k$Outbreak_date),0)
        k$MS.intensity = ifelse(max(cww$MS.diff)>0,cww$Mass_screening[which(cww$Date==min(cww$Date[which(cww$MS.diff>0)]))],
                                max(cww$Mass_screening))
        k1 = rbind(k1,k)
      }
    }
    
    return(k1)
  }))
  all = rbind(all,f)
}

all = all[which(all$total_cases>=60),]
variant = read.csv(paste0(dir,"/Outbreak_information_all.csv"),stringsAsFactors = F)
pop = read.csv(paste0(dir,"/city_pop.csv"),stringsAsFactors = F)
all = do.call(rbind,lapply(split(all,all$citycode),function(f){
  ID = variant[which(variant$citycode==unique(f$citycode)),]
  ID$start = as.Date(ID$start)
  p = pop[which(pop$ID==unique(f$citycode)),]
  f$variant = 0
  for(i in 1:nrow(f)){
    f$variant[i] = ID$variant[which(ID$start==f$Outbreak_date[i])]
  }
  f$pop_density = p$pop/p$Land.Area
  f$GDP = p$GDP
  f$pop = p$pop
  f$Doctor.number = p$Number.of.Licensed..Assistant..Doctors..person.
  
  return(f)
}))

all$Outbreak_date = as.Date(all$Outbreak_date)
all$transmission.period = as.numeric(all$Outbreak_date-all$transmission_start)
all$transmission.to.plateau = as.numeric(all$plateau_start-all$transmission_start)
all$plateau.to.decline = as.numeric(all$decline_start-all$plateau_start)

all$lock = all$Lockdown.ts+all$transmission_start
all$lock.time = ifelse(all$lock>=all$transmission_start&all$lock<all$plateau_start,1,
                       ifelse(all$lock>=all$plateau_start&all$lock<all$decline_start,2,
                              ifelse(all$lock>=all$decline_start&all$lock<all$end,3,4)))

all$decline.to.end = as.numeric(all$plateau_start-all$end)
all$mass.screening = all$MS.ts+all$transmission_start
all$mass.screening.time = ifelse(all$mass.screening>=all$transmission_start&all$mass.screening<all$plateau_start,1,
                       ifelse(all$mass.screening>=all$plateau_start&all$mass.screening<all$decline_start,2,
                              ifelse(all$mass.screening>=all$decline_start&all$mass.screening<all$end,3,4)))
all$VG = "delta"
all$VG[which(all$variant=="original"|all$variant=="alpha")] = "original&alpha"
all$VG[which(all$variant=="omicron")] = "omicron"

write.csv(all,paste0(outdir,"/pandemic_shape.csv"),row.names = F)

##################Group with geodetector###########
library(dplyr)
library(geodetector)
library(ggplot2)
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
dir = "1025dataset"
group = read.csv(paste0(dir,"/pandemic_shape.csv"),stringsAsFactors = F)
groupnew = strata_data(group,dir)
write.csv(groupnew,paste0(dir,"/alldata_group_Qstatistics.csv"),row.names = F)
