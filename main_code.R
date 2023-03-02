library(openxlsx)
library(zoo)

cases_adjusted = function(cases_data){
  cases = read.csv(cases_data,stringsAsFactors = F)
  ###########按防范区要求，修订NPI###############
  cases$Business_Premises_Closure_intensity[which(cases$Business_Premises_Closure==0)] = 4
  cases$Business_Premises_Closure_intensity[which(cases$Lockdown>=2&cases$Business_Premises_Closure<=1&
                                                    cases$Date>="2021-09-01")] = 2
  cases$Business_Premises_Closure[which(cases$Lockdown>=2&cases$Business_Premises_Closure<=1&
                                          cases$Date>="2021-09-01")] = 4
  
  cases$Public_Transportation_Closure_intensity[which(cases$Public_Transportation_Closure==0)] = 4
  
  cases$School_Closure_intensity[which(cases$School_Closure==0)] = 4
  cases$School_Closure_intensity[which(cases$Lockdown<=3&cases$School_Closure==0&cases$Lockdown>=2&
                                         cases$Date>="2021-09-01")] = 2
  cases$School_Closure_intensity[which(cases$Lockdown==4&cases$School_Closure==0&
                                         cases$Date>="2021-09-01")] = 3
  cases$School_Closure_intensity[which(cases$Lockdown==5&cases$School_Closure==0&
                                         cases$Date>="2021-09-01")] = 4
  cases$School_Closure[which(cases$Lockdown>=2&cases$School_Closure==0&
                               cases$Date>="2021-09-01")] = 2
  
  cases$Mass_screening_intensity[which(cases$Lockdown>=2&cases$Mass_screening<2&
                                       cases$Mass_screening_intensity<3&cases$Date>="2021-09-01")] = 3
  cases$Mass_screening[which(cases$Lockdown>=2&cases$Mass_screening<2&
                               cases$Date>="2021-09-01")] = 2
  
  cases$Gathering_restriction_intensity[which(cases$Gathering_restriction==0)] = 4
  cases$Gathering_restriction_intensity[which(cases$Lockdown>=2&cases$Gathering_restriction==0&
                                                cases$Date>="2021-09-01")] = 2
  cases$Gathering_restriction[which(cases$Lockdown>=2&cases$Gathering_restriction==0&
                                      cases$Date>="2021-09-01")] = 3
  
  cases$Workplace_Closure_intensity[which(cases$Workplace_Closure==0)] = 4
  cases$Workplace_Closure_intensity[which(cases$Lockdown<=3&cases$Workplace_Closure==0&cases$Lockdown>=2&
                                            cases$Date>="2021-09-01")] = 2
  cases$Workplace_Closure_intensity[which(cases$Lockdown==4&cases$Workplace_Closure==0&
                                            cases$Date>="2021-09-01")] = 3
  cases$Workplace_Closure_intensity[which(cases$Lockdown==5&cases$Workplace_Closure==0&
                                            cases$Date>="2021-09-01")] = 4
  cases$Workplace_Closure[which(cases$Lockdown>=2&cases$Workplace_Closure==0&
                                  cases$Date>="2021-09-01")] = 2
  
  cases$Logistics_Management_intensity[which(cases$Logistics_Management==0)] = 4
  
  cases$Medicine_Management_intensity[which(cases$Medicine_Management==0)] = 4
  
  cases1 = do.call(rbind,lapply(split(cases,cases$citycode),function(c){
    c$Lockdown = as.numeric(c$Lockdown)
    c = c[,c(1:22,24:28)]
    c1 = do.call(rbind,lapply(split(c,c$Date),function(d){
      if(nrow(d)>1){
        d = d[which(d$New_cases==d$New_asymptomatic+d$New_confirmed-d$Confirmed_to_asymptomatic),] 
        d = unique(d)
        d = d[complete.cases(d),]
        d = d[order(d$Lockdown),]
        d = d[1,]
      }
      return(d)
    }))
    return(c1)
  }))
  write.csv(cases1,"cases&NPI_adjusted.csv",row.names = F)
  return(cases1)
}

province_vaccination = function(vacW,vac){
  vac = subset(vac,vac$location=="China")
  vacW$one_dose_weight = round(vacW$one_dose/sum(vacW$one_dose),5)
  vacW$fully_vaccination_weight = round(vacW$fully_vaccination/sum(vacW$fully_vaccination),5)
  vacW$booster.dose_weight= round(vacW$booster.dose/sum(vacW$booster.dose),5)
  vac$date = as.Date(vac$date)
  vac = vac[order(vac$date),]
  vac$people_vaccinated[1] = 0
  vac$people_vaccinated = na.locf(na.locf(vac$people_vaccinated),fromLast=TRUE)
  
  vac$people_fully_vaccinated[1] = 0 
  vac$people_fully_vaccinated = na.locf(na.locf(vac$people_fully_vaccinated),fromLast=TRUE)
  
  vac$total_boosters[1] = 0
  vac$total_boosters= na.locf(na.locf(vac$total_boosters),fromLast=TRUE)
  
  provac = do.call(rbind,lapply(split(vacW,vacW$province),function(p){
    v = data.frame(province = p$province , date = vac$date,
                   total_vaccination_rate = vac$people_vaccinated*p$one_dose_weight/(p$pop*10000),
                   fully_vaccination_rate = vac$people_fully_vaccinated*p$fully_vaccination_weight/(p$pop*10000),
                   boost_dose_rate = vac$total_boosters*p$booster.dose_weight/(p$pop*10000))
  }))
  return(provac)
}
merge_cases_data = function(dir){
  
  case1 = read.xlsx(paste0(dir,"3 Daily new cases_21Apr2022V1.xlsx"))[,c(1:10)]
  case2 = read.csv(paste0(dir,"0607.csv"),stringsAsFactors = F)[,c(1:9)]
  sel = read.csv(paste0(dir,"1 Outbreak_information_all_V1.csv"),stringsAsFactors = F)
  pop = read.csv(paste0(dir,"city_pop.csv"),stringsAsFactors = F)
  case2$'城市代码' = NA
  case1$Date = as.Date(as.character(case1$'日期'),format = "%Y%m%d")
  case2$Date = as.Date(case2$日期)
  colnames(case1) = c("city","province","ID","date","year","month","day","New_cases",
                      "New_asymptomatic","New_confirmed","Date")
  colnames(case2) = c("province","city","date","year","month","day",
                      "New_asymptomatic","New_confirmed","Confirmed_from_asymptomatic","ID","Date")
  
  case2 = do.call(rbind,lapply(split(case2,case2$province),function(c){
    c = c[,c("Date","city","province","ID",
             "New_asymptomatic","New_confirmed","Confirmed_from_asymptomatic")]
    if(c$province==as.character("北京市")|c$province==as.character("上海市")|
       c$province==as.character("天津市")|c$province==as.character("重庆市")){
      c = c%>%group_by(Date)%>%summarise(city=unique(province),province = unique(province),ID = NA,
                                         New_asymptomatic = sum(New_asymptomatic),New_confirmed = sum(New_confirmed),
                                         Confirmed_from_asymptomatic=sum(Confirmed_from_asymptomatic))
    }
    return(c)
  }))
  
  #NPI = read.csv(paste0(dir,"NPI_all_0705.csv"),stringsAsFactors = F)
  cases = do.call(rbind,lapply(split(case2,case2$city),function(c2){
    se = subset(sel,sel$city==unique(c2$city))
    #NPIc =subset(NPI,NPI$city==unique(c2$city))
    if(nrow(se)>0){
      popc = pop$pop[which(pop$ID==unique(se$citycode))]
      c1 = subset(case1,case1$city==unique(c2$city))
      c2$New_cases = c2$New_asymptomatic+c2$New_confirmed-c2$Confirmed_from_asymptomatic
      c1$Confirmed_from_asymptomatic = c1$New_confirmed+c1$New_asymptomatic-c1$New_cases
      c1 = c1[,c("city","province","ID","Date",
                 "New_cases","New_asymptomatic","New_confirmed","Confirmed_from_asymptomatic")]
      c2 = c2[,c("city","province","ID","Date",
                 "New_cases","New_asymptomatic","New_confirmed","Confirmed_from_asymptomatic")]
      cityID = sel$citycode[which(sel$city==unique(c2$city))]
      c2$ID=unique(cityID)
      
      c1 = c1[order(c1$Date),]
      c2 = c2[order(c2$Date),]
      c1 = c1[which(c1$New_cases==c1$New_asymptomatic+c1$New_confirmed-c1$Confirmed_from_asymptomatic),]
      c2 = c2[which(c2$New_cases==c2$New_asymptomatic+c2$New_confirmed-c2$Confirmed_from_asymptomatic),]
      c1 = unique(c1)
      c2 = unique(c2)
      c = do.call(rbind,list(c1,c2[which(c2$Date>max(c1$Date)),]))
      
      c[is.na(c)]<-0
      c = unique(c)
      #colnames(NPIc)[1] = "Date"
      #NPIc$Date=as.Date(NPIc$Date)
      #ALL = merge(NPIc,c,by=c("Date"),all.x=T)
      #ALL$pop = popc
      #return(ALL)
      f = do.call(rbind,lapply(seq(1,nrow(se)),function(k){
        k = se[k,]
        if(sum(c$New_cases[which(c$Date>=as.Date(k$start)&c$Date<=as.Date(k$end))])>60){
          return(c[which(c$Date>=as.Date(k$start)&c$Date<=as.Date(k$end)),])
        }else{return(NULL)}
      }))
      return(f)
    }
  }))
  
  #write.csv(cases,"all_cases.csv",row.names = F)
  return(cases)
}
smooth_Gaussian <-function(cases){
  
  
  r <- 7
  
  cases_smoothed <- cases
  
  GaussTemp <- c(0.0004,0.0022,0.0088,0.0270,0.0648,0.1210,0.1761,0.1995,0.1761,0.1210,0.0648,0.0270,0.0088,0.0022,0.0004)

  for (i in (8:(length(cases)-r))){
      cases_smoothed[i] <- sum(cases[c((i-r):(i+r))]*t(GaussTemp))
  }
  return(cases_smoothed)
}
smooth_NPI = function(NPIname,step = 7,data){
  NPI = do.call(rbind,lapply(split(data,data$city),function(w){
    npi.w = do.call(rbind,lapply(split(w,w$original_start),function(f){
      npi.ww = do.call(rbind,lapply(split(f,f$wave),function(w1){
        npi = w1[,NPIname]
        npi[is.na(npi)]=0
        #npio = npi
        diff = which(diff(npi)!=0)
        orn = npi[diff]
        ord = npi[diff+1]
        len = length(npi)
        npi1 = npi
        if(length(diff)>1){
          s = diff[1]
          npi[s:(s+step)] = seq(orn[1],ord[1],(ord[1]-orn[1])/step)
          for(i in seq(2,length(diff))){
            s = diff[i]
            if(ord[i]>orn[i]){
              if(ord[i-1]>orn[i-1]){
                npi[s:(s+step)] = npi[s:(s+step)] + seq(orn[i],ord[i],(ord[i]-orn[i])/(step)) -orn[i]
                npi[is.na(npi)] = 0
                if(length(which(npi[s:(s+step)]>max(orn[i],ord[i])))>0){
                  start = min(which(npi[s:(s+step)]>max(orn[i],ord[i])))
                  npi[s:(s+step)][start:length(npi[s:(s+step)])] = npi1[s:(s+step)][start:length(npi[s:(s+step)])]
                  npi = npi[1:len]
                }else{
                  npi = npi[1:len]
                }
              }else{
                npi[s:(s+step)] = seq(orn[i],ord[i],(ord[i]-orn[i])/(step))
                npi[is.na(npi)] = 0
                if(length(which(npi[s:(s+step)]>max(orn[i],ord[i])))>0){
                  start = min(which(npi[s:(s+step)]>max(orn[i],ord[i])))
                  npi[s:(s+step)][start:length(npi[s:(s+step)])] = npi1[s:(s+step)][start:length(npi[s:(s+step)])]
                  npi = npi[1:len]
                }else{
                  npi = npi[1:len]
                }
              }
            }else{
              npi[s:(s+step)] = seq(orn[i],ord[i],(ord[i]-orn[i])/(step))
              npi[is.na(npi)] = 0
              npi = npi[1:len]
            }
          }
        }else if(length(diff)>0){
          s = diff
          npi[s:(s+step)] = seq(orn[1],ord[1],(ord[1]-orn[1])/(step))
          npi = npi[1:len]
        }else{
          npi = npi[1:len]
        }
        
        npi = data.frame(NPI=npi,Date = w1$Date)
        npi$city = unique(w1$city)
        npi$original_start = unique(w1$original_start)
        npi$wave = unique(w1$wave)
        return(npi)
      }))
      return(npi.ww)
    }))
    return(npi.w)
  }))
  return(NPI) 
}

strata_data<-function (group,dir){
  pdf(paste0(dir,"/Starta.pdf"))
  gdata = data.frame()
  for (Dataset in split(group,group$VG)){
    Dataset$rate_cases <- Dataset$total_cases/Dataset$pop*100
    q_rate_case <- lapply(seq(round(quantile(Dataset$rate_cases)[2],-1), 
                              #round(quantile(Dataset$rate_cases)[4],-1),
                              max(round(max(Dataset$rate_cases)*1/20,-2),round(quantile(Dataset$rate_cases)[4],0)),
                              by=10), function(v) {
                                d <- data.frame(r=Dataset$rate_cases, strat=Dataset$rate_cases > v)
                                c(v, geodetector::factor_detector("r","strat",d)[[1]][1,1],
                                  geodetector::factor_detector("r","strat",d)[[1]][1,2])
                              }) %>% do.call(rbind,.)
    q_rate_last <- lapply(seq(round(quantile(Dataset$last)[2],0), 
                              #=round(quantile(Dataset$last)[4],0),
                              max(Dataset$last)-5,
                              by=5), function(v) {
                                d <- data.frame(r=Dataset$last, strat=Dataset$last > v)
                                c(v, geodetector::factor_detector("r","strat",d)[[1]][1,1],
                                  geodetector::factor_detector("r","strat",d)[[1]][1,2])
                              }) %>% do.call(rbind,.)
    q_rate_case<-as.data.frame(q_rate_case)
    q_rate_last<-as.data.frame(q_rate_last)
    
    q1<-q_rate_case$V1[which(q_rate_case$V2==max(q_rate_case$V2))]
    q2<-q_rate_last$V1[which(q_rate_last$V2==max(q_rate_last$V2))]
    Dataset$Group <- "1"
    Dataset$Group[Dataset$rate_cases <= max(q1) & Dataset$last > max(q2)] <- "2"
    Dataset$Group[Dataset$rate_cases > max(q1) & Dataset$last <= max(q2)] <- "3"
    Dataset$Group[Dataset$rate_cases > max(q1) & Dataset$last > max(q2)] <- "4"
    Dataset$q1 = max(q1)
    Dataset$Q.value_q1 = unique(q_rate_case$V2[which(q_rate_case$V2==max(q_rate_case$V2))])
    Dataset$P.value_q1 = unique(q_rate_case$V3[which(q_rate_case$V2==max(q_rate_case$V2))])
    Dataset$q2 = max(q2)
    Dataset$Q.value_q2 = unique(q_rate_last$V2[which(q_rate_last$V2==max(q_rate_last$V2))])
    Dataset$P.value_q2 = unique(q_rate_last$V3[which(q_rate_last$V2==max(q_rate_last$V2))])
    
    gdata = rbind(gdata,Dataset)
    g<-ggplot(Dataset, aes(rate_cases,last)) +
      geom_point(aes(color=Group),size=2) + xlab("Infection rate(per 100,000)") +
      scale_color_manual(values=c("#ffdb47","#8da636","#1581c7","#a52a2a"))+ ylab("Duration(Days)") +
      geom_hline(yintercept = max(q2), color="#6a6a6a") +
      geom_vline(xintercept = max(q1), color="#6a6a6a")+
      theme_bw()+
      labs(title = unique(Dataset$VG))+
      theme(legend.position = "right",
            legend.title=element_blank(),
            plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
            axis.title.x= element_text(color="black",size = 10),
            axis.text.x = element_text(color="black",size = 10),
            axis.title.y= element_text(color="black",size = 10),
            axis.text.y = element_text(color="black",size = 10),
            plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
    
    if(unique(Dataset$VG)=="omicron"){
      g = g +coord_trans(x = "log")
    }
    #print(g)
    ggsave(paste0(dir,"/",unique(Dataset$VG),"_starta.pdf"),g, width=70, height=60, units="mm",scale=2)
  }
  #dev.off()
  return(gdata)
}

