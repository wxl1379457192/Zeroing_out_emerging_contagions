library(ggplot2)
library(openxlsx)
library(reshape)
library(ggalt)
library(ggtext)
library(scales) # to access break formatting functions
library(ggpubr)
library(ggthemes)
setwd("E:/zeroCOVID_NPI/Version0504")
dir = "SEIR_P2/"
element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL, hi.family = NULL) {
  structure(
    c(element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col, hi.family = hi.family)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}
######################S1 plot#####################
conlist = list.files(dir,pattern="^Scenario_")
indir = "dataset"
dat = read.csv(paste0(indir,"/Rt&smooth_NPI.csv"),stringsAsFactors = F)
dat$original_start = as.Date(dat$original_start)
S1 =  do.call(rbind,lapply(conlist,function(c){
  f = read.xlsx(paste0(dir,c), sheet = 1)
  f = do.call(rbind,lapply(split(f,f$Scenario),function(s){
    s$original_start = as.Date(as.character(s$original_start),format="%Y%m%d")
    s$Days = seq(1,nrow(s))
    s$Date = s$original_start+s$Days
    s$Mean = round(s$Mean,2)
    s$CI05 = round(s$CI05,2)
    s$CI95 = round(s$CI95,2)
    #s$CI05[which(s$CI05<0)]=0
    #s$CI95[which(s$CI95<0)]=0
    #s = s[which(s$Mean>=1|s$CI05>=1|s$CI95>=1),]
    s$pop = unique(dat$pop[which(dat$citycode==unique(s$citycode))])
    s$variant = unique(dat$VG[which(dat$citycode==unique(s$citycode)&dat$original_start==unique(s$original_start))])
    s$cummean = cumsum(s$Mean)
    s$cumCI05 = cumsum(s$CI05)
    s$cumCI95 = cumsum(s$CI95)
    s$Duration = nrow(s)
    s$id  = strsplit(c,".xlsx")[[1]][1]
    return(s)
  }))
  d = subset(dat,dat$citycode==unique(f$citycode))
  d = subset(d,d$original_start==unique(f$original_start))
  d1 = data.frame(Mean = d$New_cases_smooth,CI05=d$New_cases_smooth,
             CI95=d$New_cases_smooth,citycode=unique(f$citycode),
             original_start=unique(f$original_start),Scenario="real",
             Days = seq(1,nrow(d)),
             Date = as.Date(unique(f$original_start))+seq(1,nrow(d)),
             pop=unique(f$pop),variant=unique(f$variant),
             Duration = unique(f$Duration),id=unique(f$id))
  d1$cummean=cumsum(d1$Mean)
  d1$cumCI05 = cumsum(d1$CI05)
  d1$cumCI95 = cumsum(d1$CI95)
  f = do.call(rbind,list(f,d1))
  return(f)
}))
library(dplyr)
maxreal = max(S1$Duration[which(S1$Scenario=="realworld")])
S2 = subset(S1,S1$Days<=maxreal)
S2 = subset(S2,S2$Scenario!="real")
all = S2%>%group_by(Days,Scenario,variant)%>%summarise(Mean =sum(Mean),CI05 = sum(CI05),CI95 = sum(CI95))
all = do.call(rbind,lapply(split(all,all$variant),function(al){
  al = do.call(rbind,lapply(split(al,al$Scenario),function(k){
    k$cummean = cumsum(k$Mean)
    k$cumCI05 = cumsum(k$CI05)
    k$cumCI95 = cumsum(k$CI95)
    return(k)
  }))
  return(al)
}))
All= do.call(rbind,lapply(split(all,all$variant),function(al){
  al = do.call(rbind,lapply(split(al,al$Scenario),function(k){
    f = data.frame(Scenario = unique(k$Scenario),
                   variant = unique(k$variant),
                   AUC = sum(k$cummean),
                   AUC_CI05 = sum(k$cumCI05),
                   AUC_CI95 = sum(k$cumCI95))
    return(f)
  }))
  al$AUC_CI05 = al$AUC_CI05/al$AUC[which(al$Scenario=="Baseline")]
  al$AUC_CI95 = al$AUC_CI95/al$AUC[which(al$Scenario=="Baseline")]
  al$AUC = al$AUC/al$AUC[which(al$Scenario=="Baseline")]
  al = subset(al,al$Scenario!="Baseline")
  return(al)
}))
All$variant = factor(All$variant,levels = c("omicron","delta","original&alpha"))

g = ggplot(All,aes(x=Scenario,y=AUC,fill=variant))+
  geom_bar(stat = "identity",position=position_dodge(width = 0.95),alpha=0.8)+
  geom_errorbar(aes(ymin = AUC_CI05,ymax = AUC_CI95,width=.1,color = variant),
                stat = "identity",position=position_dodge(width = 0.95))+
  scale_color_manual(values=c("omicron" = "#257187",
                              "original&alpha" = "#99BAAC",
                              "delta" = "#A49B90"),
                     labels=c("omicron"="Omicron",
                              "original&alpha" = "Pre-Delta",
                              "delta" = "Delta"))+
  theme_grey()+
  scale_fill_manual(values=c("omicron" = "#257187",
                             "original&alpha" = "#99BAAC",
                             "delta" = "#A49B90"),
                    labels=c("omicron"="Omicron",
                             "original&alpha" = "Pre-Delta",
                             "delta" = "Delta"))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),limits=c(0,0.4))+
  labs(x = "",y ="Ratio of area under the curve")+
  theme(legend.position ="none",
        legend.background = element_rect(
          fill = "transparent", # 填充色
          colour = "transparent", # 框线色
          size = 1.5),
        legend.title=element_blank(),
        panel.grid = element_blank(),
        axis.line=element_line(color="grey90",size=0.5),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size =unit(9,"pt")),
        axis.title.x= element_text(color="black",size = unit(10,"pt")),
        axis.text.x = element_text(color="black",size = unit(9,"pt")),
        axis.title.y= element_text(color="black",size = unit(10,"pt")),
        axis.text.y = element_text(color="black",size = unit(9,"pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))

outdir = "Figures_Reviewed"
ggsave(paste0(outdir,"/fig2c_AUC.pdf"),g, width=16, height=7, units="cm", scale=1)


All = all
outdir = "Figures_Reviewed"
#"#257187","#A49B90","#99BAAC","#B4281F"
for (all in split(All,All$variant)){
  S = subset(S1,S1$variant==unique(all$variant))
  pop = sum(unique(S$pop))
  g = ggplot(all,aes(x=Days,y=cummean,color=Scenario))+
    geom_xspline()+#facet_wrap(.~variant,scales = "free")+
    geom_ribbon(aes(x=Days, ymin=cumCI05, ymax=cumCI95,fill=Scenario),color = "transparent",alpha=0.1)+
    labs(y="Daily new infections")+
    scale_y_log10(expand = expansion(mult = c(0,0.1)),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    scale_fill_manual(values = c("Baseline"="grey30",
                                 "realworld"="#452103",
                                 "withoutFM" = "#335F70",
                                 "withoutCT" = "#2A9D8F",
                                 "withoutSD" = "#E9C46A",
                                 "withoutPCR" = "#F4A261"),
                      labels=c("realworld"="All NPIs",
                               "Baseline"="No NPIs",
                                 "withoutCT" = "no CT",
                               "withoutFM" = "no FM",
                               "withoutSD" = "no SD",
                               "withoutPCR" = "no PCR"))+
    scale_color_manual(values = c("Baseline"="grey30",
                                  "realworld"="#452103",
                                  "withoutFM" = "#335F70",
                                  "withoutCT" = "#2A9D8F",
                                  "withoutSD" = "#E9C46A",
                                  "withoutPCR" = "#F4A261"),
                       labels=c("realworld"="All NPIs",
                                "Baseline"="No NPIs",
                                "withoutCT" = "no CT",
                                "withoutFM" = "no FM",
                                "withoutSD" = "no SD",
                                "withoutPCR" = "no PCR"))+
    guides(fill = guide_legend(nrow=3),color = guide_legend(nrow=3))+
    scale_x_continuous(expand = expansion(mult = c(0,0)),
                       limits = c(0,max(all$Days)+10))+
    theme_clean()+
    theme(legend.position ="none",
          legend.background = element_rect(
            fill = "transparent", # 填充色
            colour = "transparent", # 框线色
            size = 1.5),
          legend.title=element_blank(),
          panel.grid = element_blank(),
          axis.line=element_line(color="grey90",size=0.5),
          plot.title = element_text(color="black",hjust = 0,vjust=0,size =unit(9,"pt")),
          axis.title.x= element_text(color="black",size = unit(10,"pt")),
          axis.text.x = element_text(color="black",size = unit(9,"pt")),
          axis.title.y= element_text(color="black",size = unit(10,"pt")),
          axis.text.y = element_text(color="black",size = unit(9,"pt")),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
    g = g+annotate("text", x=max(all$Days)+4,
             y=max(all$cummean[which(all$Scenario=="withoutPCR")]),
             color = "#F4A261",
             label="no PCR",size =unit(3,"pt"))+
   # geom_segment(aes(x = max(all$Days)+1, xend = max(all$Days)+7,
      #               y =  max(all$cummean[which(all$Scenario=="withoutPCR")])-2, 
  #                   yend = max(all$cummean[which(all$Scenario=="withoutPCR")])), 
      #           colour = "#8CBCBD", linetype = 7,size=0.2, 
     #            arrow= arrow(ends = "both",length = unit(0.01, "npc"))) +
    annotate("text",x=max(all$Days)+4,
               y=max(all$cummean[which(all$Scenario=="realworld")])+2000,
               color = "#452103",
               label="all NPIs",size =unit(3,"pt"))+
    annotate("text",x=max(all$Days)+4,
               y=max(all$cummean[which(all$Scenario=="withoutCT")])+5000,
               color = "#2A9D8F",
               label="no CT",size =unit(3,"pt"))+
    annotate("text",x=max(all$Days)+4,
             y=max(all$cummean[which(all$Scenario=="withoutFM")])+5000,
             color = "#335F70",
             label="no FM",size =unit(3,"pt"))+
    annotate("text", x=max(all$Days)+4,
             y=max(all$cummean[which(all$Scenario=="withoutSD")])+5000,
             color = "#E9C46A",
             label="no SD",size =unit(3,"pt"))+
    annotate("text", x=max(all$Days)+4,
             y=max(all$cummean[which(all$Scenario=="Baseline")])+50000,
             color = "grey30",
             label="no NPIs",size =unit(3,"pt"))+
    geom_hline(aes(yintercept = pop*10000),lty="dashed",color= "#7E6148",size = 0.8)+
    labs(x = "Days",y ="Cumulative number of infections")
    if(unique(all$variant)=="original&alpha"){
      ggsave(paste0(outdir,"/fig3a_preDelta_4npi.pdf"),g, width=6, height=5.5, units="cm", scale=1)
    }else{
      ggsave(paste0(outdir,"/fig3a_",unique(all$variant),"_4npi.pdf"),g, width=6, height=5.5, units="cm", scale=1)
    }
    #dev.off()
}


#g = g + geom_line(data = S1,aes(x=Days,y=cummean,color = Scenario, linetype = id),
#                     show.legend = FALSE, size=0.05, alpha = 0.2)

pd = do.call(rbind,lapply(split(All,All$variant),function(s){
  s = do.call(rbind,lapply(split(s,s$Scenario),function(sv){
    data.frame(Infections = max(sv$cummean), 
               Infections_CI05 = max(sv$cumCI05),
               Infections_CI95 = max(sv$cumCI95),
               Scenario = unique(sv$Scenario),
               variant = unique(sv$variant)
               ) 
  }))
  
  S = subset(S1,S1$variant==unique(s$variant))
  s$pop = sum(unique(S$pop))
  #s$infection_ratio = s$Infections/(s$pop*10000)*100
  
  s$effect =(s$Infections[which(s$Scenario=="Baseline")]-s$Infections)/
    (s$Infections[which(s$Scenario=="Baseline")])*100
  
  s$effectCI05 = (s$Infections_CI05[which(s$Scenario=="Baseline")]-s$Infections_CI05)/
    (s$Infections_CI05[which(s$Scenario=="Baseline")])*100
  
  s$effectCI95 = (s$Infections_CI95[which(s$Scenario=="Baseline")]-s$Infections_CI95)/
    (s$Infections_CI95[which(s$Scenario=="Baseline")])*100
  
  #s$effect_real =(s$Infections[which(s$Scenario=="Baseline")]-s$Infections)/(s$pop*10000)*100
  #s$infection_decline = s$Infections[which(s$Scenario=="Baseline")]-s$Infections
  return(s)
}))
write.csv(pd,paste0(outdir,"/fig3data.csv"),row.names = F)


