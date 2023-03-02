setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
source('code/main_code.R')
datadir = "dataset/auxdata"
filenames = list.files(datadir,pattern="csv$")   
alldatas = data.frame()

for (i in filenames){
  Hum = read.csv(paste0(datadir,"/",i),stringsAsFactors = F)
  Hum = Hum[,-c(1,ncol(Hum))]
  datelist = colnames(Hum)[-ncol(Hum)]
  startdate = as.Date(substring(strsplit(i,"series")[[1]][2],1,10))
  for (n in datelist){
    datenum = as.numeric(gsub("X","",n))
    d = as.data.frame(Hum[,c(n,"PAC")])
    colnames(d) = c("value","ID")
    d$Date = startdate+datenum
    d$variable = strsplit(i,"_time_series")[[1]][1]
    alldatas = rbind(alldatas,d)
    print(paste(i,startdate+datenum,"has been process"))
  }
}

Temp = subset(alldatas,alldatas$variable=="Temp")
colnames(Temp)[1] = "Temp"
Hum= subset(alldatas,alldatas$variable=="Humidity")
colnames(Hum)[1] = "Hum"
all = merge(Temp[,-4],Hum[,-4],id=c("ID","Date"))
#write.csv(all,"dataset/auxdata.csv",row.names = F)
###################疫苗&环境变量数据拼接，实际疫苗接种率计算########################
#all = read.csv("dataset/auxdata.csv",stringsAsFactors = F)
dir ="dataset"
indir = "1025dataset"
vacW = read.csv(paste0(dir,"/vaccinations_weighted.csv"),stringsAsFactors = F)
vac = read.csv(paste0(dir,"/vaccinations.csv"),stringsAsFactors = F)
vaccination.rate = province_vaccination(vacW,vac)
cases = read.csv(paste0(indir,"/NPI&cases_smooth&Rt.csv"),stringsAsFactors = F)
ID = read.csv(paste0(dir,"/1 Outbreak_information_all_V1.csv"),stringsAsFactors = F)
cases.envir = do.call(rbind,lapply(split(cases,cases$citycode),function(c){
  id = subset(ID,ID$citycode==unique(c$citycode))
  e = subset(all,all$ID==unique(c$citycode))
  c$Date = as.Date(c$Date)
  e = e[,2:4]
  c1 = merge(c,e,id = "Date",all.x = T)
  vc = subset(vaccination.rate,vaccination.rate$province==unique(id$province))[,-1]
  colnames(vc)[1] = "Date"
  c1 = merge(c1,vc,id = "Date",all.x = T)
  c1[is.na(c1)] = 0
  return(c1)
}))



cases.envir$'Practical.vaccination' = cases.envir$boost_dose_rate*(86+79)/200+(cases.envir$fully_vaccination_rate - 
                                                                                 cases.envir$boost_dose_rate)*(50+51)/200
write.csv(cases.envir,paste0(indir,"/NPI&cases_smooth&Rt&env&Vac.csv"),row.names = F)




