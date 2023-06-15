library("rstan")
library("bayesplot")
rm(list = ls())
setwd("E:/zeroCOVID_NPI/Version0504")
alldir = "Reviewed_overallNPI_stanV5_addpop"

pd = read.csv(paste0(alldir,"/testdata_IndNPI.csv"),stringsAsFactors = F)
NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")
auxname = c("Temp","fully_vaccination_rate","pop_density")
ObNPI = c("Static.management",
          "Facial_Mask",
          "PCR.screening",
          "Contact.control") 
conlist = list.files(alldir,pattern="*.rds")
f = readRDS(paste0(alldir,"/",conlist [[1]]))
out = rstan::extract(f)
al=as.data.frame(out$alpha)
colnames(al)=NPIname
for (i in NPIname){
  al[,i] = al[,i]*mean(pd[,i])
}
al$Static.management = al$Lockdown+al$Business_Premises_Closure+al$Public_Transportation_Closure+
  al$Gathering_restriction+al$Workplace_Closure+al$School_Closure
al$PCR.screening = al$Medicine_Management+al$Mass_screening
colnames(al)[which(colnames(al)=="Control")] = "Contact.control"
al= al[,ObNPI]
for (i in c(ObNPI)){
  al[,i] = (1-exp(-al[,i]))
}

syn=mcmc_intervals_data(al,prob = .5,prob_outer= .95,point_est="mean")
#al1 = data.frame(t(al))
#al1$var = colnames(al)
syn$variant = "Overall"

################################
dir =  "Reviewed_IndividualNPI_stanV5_pooled_onebeta"
leaveonedir = "Reviewed_Bayesian_model_validation"
conlist = list.files(dir,pattern="*.rds")
all = read.csv(paste0(dir,"/testdata_IndNPI.csv"),stringsAsFactors = F)
all = all[which(all$Mean.R.>0),]


f = readRDS(paste0(dir,"/",conlist))
out = rstan::extract(f)
al=as.data.frame(out$alpha)
variant = c("original&alpha","delta","omicron")
val1 = do.call(rbind,lapply(seq(1,3),function(x){
  v = all[which(all$VG==variant[x]),]
  alv = al[,((x-1)*length(NPIname)+1):(x*length(NPIname))] 
  colnames(alv)=NPIname
  for (i in NPIname){
    alv[,i] = alv[,i]*mean(v[,i])
  }
  alv$Static.management = alv$Lockdown+alv$Business_Premises_Closure+alv$Public_Transportation_Closure+
    alv$Gathering_restriction+alv$Workplace_Closure+alv$School_Closure
  alv$PCR.screening = alv$Medicine_Management+alv$Mass_screening
  colnames(alv)[which(colnames(alv)=="Control")] = "Contact.control"
  alv= alv[,ObNPI]
  for (i in c(ObNPI)){
    alv[,i] = (1-exp(-alv[,i]))
  }
  alv$variant = variant[x]
  return(alv)
}))

valmodel = list.files(leaveonedir,pattern="*.rds")
val2 = do.call(rbind,lapply(valmodel,function(model){
  f = readRDS(paste0(leaveonedir,"/",model))
  out = rstan::extract(f)
  al=as.data.frame(out$alpha)
  
  citycode = substring(strsplit(model,"_")[[1]][3],1,6)
  ordate = substring(strsplit(model,"_")[[1]][3],8,17)
  data = all[which(all$citycode==citycode),]
  data = subset(data,data$original_start==ordate)
  data$Date=as.Date(data$Date)
  data = data[order(data$Date),]
  
  val2 = do.call(rbind,lapply(seq(1,3),function(x){
    v = data[which(data$VG==variant[x]),]
    alv = al[,((x-1)*length(NPIname)+1):(x*length(NPIname))] 
    colnames(alv)=NPIname
    for (i in NPIname){
      alv[,i] = alv[,i]*mean(v[,i])
    }
    alv$Static.management = alv$Lockdown+alv$Business_Premises_Closure+alv$Public_Transportation_Closure+
      alv$Gathering_restriction+alv$Workplace_Closure+alv$School_Closure
    alv$PCR.screening = alv$Medicine_Management+alv$Mass_screening
    colnames(alv)[which(colnames(alv)=="Control")] = "Contact.control"
    alv= alv[,ObNPI]
    for (i in c(ObNPI)){
      alv[,i] = (1-exp(-alv[,i]))
    }
    alv$variant = variant[x]
    return(alv)
  }))
  print(model)
  return(val2)
  
}))

all = do.call(rbind,list(val1,val2))

syn2=do.call(rbind,lapply(split(all,all$variant),function(v){
  var = unique(v$variant)
  v = v[,c(1,2,3,4)]
  v = v[complete.cases(v),]
  f = mcmc_intervals_data(v,prob = .5,prob_outer= .95,point_est="mean")
  f$variant = var 
  return(f)
}))
  

##################plot######################

synall = do.call(rbind,list(syn,syn2))
synall$variant = factor(synall$variant,levels = c("omicron","delta","original&alpha","Overall"))

library(ggtext)
library(rlang)

element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL, hi.family = NULL) {
  structure(
    c(element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col, hi.family = hi.family)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}

element_grob.element_textbox_highlight <- function(element, label = "", ...) {
  if (label %in% element$hi.labels) {
    element$fill <- element$hi.fill %||% element$fill
    element$colour <- element$hi.col %||% element$colour
    element$box.colour <- element$hi.box.col %||% element$box.colour
    element$family <- element$hi.family %||% element$family
  }
  NextMethod()
}
library(ggthemes)
library(scales)
library(ggsci)
f = ggplot(synall)+
  geom_linerange(aes(xmin =ll, xmax=hh, y=variant,color= variant,size = variant == "Overall"),
                 show.legend = FALSE,alpha=0.3,
                 position=position_dodge(width = 0.7))+
  geom_hline(show.legend = TRUE,aes(yintercept = -3, color=variant,size = variant == "Overall"),alpha=0.3)+
  geom_linerange(aes(xmin =l, xmax=h, y=variant,color=variant,size = variant == "Overall"),show.legend = FALSE,
                 position=position_dodge(width = 0.8))+
  geom_pointrange(aes(x = m, y = variant, xmin=m,xmax=m,
                      size = variant == "Overall",
                      color = variant), fatten = 1.5, show.legend = TRUE,alpha=0.5,
                  position= position_dodge(width = 0.7))+
  facet_grid(parameter~.)+ 
  scale_size_manual(values = c(2,4))+
  scale_color_manual(values=c("#257187",
                              "#A49B90",
                              "#99BAAC",
                              "#B4281F"
                              #"#3773a4"#"#6b9ec7",
  ))+theme_grey()+
  scale_fill_manual(values=c("#257187",
                             "#A49B90",
                             "#99BAAC",
                             "#B4281F"))+
  scale_x_continuous(expand=c(0.01,0.01),limits=range(-0.01,0.80),labels=percent)+
  labs(x =expression("Reduction in "*R[t]),y = NULL)+
  theme(legend.position = "top",legend.key=element_rect(fill='transparent'),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = "transparent"),
        panel.grid.minor = element_line(colour = "transparent"),
        legend.title=element_blank(),
        plot.margin = margin(t = 5,  # ¶¥²¿±ßÔµ¾àÀë
                             r = 5,  # ÓÒ±ß±ßÔµ¾àÀë
                             b = 5,  # µ×²¿±ßÔµ¾àÀë
                             l = 80),
        strip.text = element_textbox_highlight(
          size = 12, face = "bold",
          fill = "white", box.color = "white", color = "gray40",
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(0, 1, 3, 1),
          hi.labels = "Summer", hi.family = "Bangers",
          hi.fill = "firebrick", hi.box.col = "firebrick", hi.col = "white"
        ))
print(f)
outdir = "Figures_Reviewed"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}
ggsave(paste0(outdir,"/fig_synergistic_effect_pooled.pdf"),f, width=16, height=16, units="cm", scale=1)

write.csv(synall,paste0(outdir,"/fig_synergistic_effect_pooled.csv"),row.names = F)

