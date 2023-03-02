library(ggplot2)
library(openxlsx)
library(reshape)
library(ggalt)
library(ggtext)
library(scales) # to access break formatting functions
library(ggpubr)
library(ggthemes)
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
dir = "SEIR_simulation_Realworld_1220_4NPI/"
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
indir = "1025dataset"
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
All = all
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
                                 "withoutCC" = "#2A9D8F",
                                 "withoutSD" = "#E9C46A",
                                 "withoutPCR" = "#F4A261"),
                      labels=c("realworld"="All NPIs",
                               "Baseline"="No NPIs",
                               "withoutCC" = "no CT",
                               "withoutFM" = "no FM",
                               "withoutSD" = "no SD",
                               "withoutPCR" = "no PCR"))+
    scale_color_manual(values = c("Baseline"="grey30",
                                  "realworld"="#452103",
                                  "withoutFM" = "#335F70",
                                  "withoutCC" = "#2A9D8F",
                                  "withoutSD" = "#E9C46A",
                                  "withoutPCR" = "#F4A261"),
                       labels=c("realworld"="All NPIs",
                                "Baseline"="No NPIs",
                                "withoutCC" = "no CT",
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
    g = g+annotate("text", x=max(all$Days)+7,
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
               y=max(all$cummean[which(all$Scenario=="withoutCC")])+5000,
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
    outdir = "Figures"
    if(unique(all$variant)=="original&alpha"){
      ggsave(paste0(outdir,"/fig3b_preDelta.pdf"),g, width=8, height=8, units="cm", scale=1)
    }else{
      ggsave(paste0(outdir,"/fig3b_",unique(all$variant),".pdf"),g, width=7, height=8, units="cm", scale=1)
    }
    #dev.off()
}


#g = g + geom_line(data = S1,aes(x=Days,y=cummean,color = Scenario, linetype = id),
#                     show.legend = FALSE, size=0.05, alpha = 0.2)

pd = do.call(rbind,lapply(split(All,All$variant),function(s){
  s = do.call(rbind,lapply(split(s,s$Scenario),function(sv){
    data.frame(Infections = max(sv$cummean), Scenario = unique(sv$Scenario),
               variant = unique(sv$variant),
               Infections_CI05 = max(sv$cumCI05),
               Infections_CI95 = max(sv$cumCI95)) 
  }))
  
  S = subset(S1,S1$variant==unique(s$variant))
  s$pop = sum(unique(S$pop))
  s$infection_ratio = s$Infections/(s$pop*10000)*100
  
  s$effect =(s$Infections[which(s$Scenario=="Baseline")]-s$Infections)/
    (s$Infections[which(s$Scenario=="Baseline")])*100
  
  s$effectCI05 = (s$Infections_CI05[which(s$Scenario=="Baseline")]-s$Infections_CI05)/
    (s$Infections_CI05[which(s$Scenario=="Baseline")])*100
  
  s$effectCI95 = (s$Infections_CI95[which(s$Scenario=="Baseline")]-s$Infections_CI95)/
    (s$Infections_CI95[which(s$Scenario=="Baseline")])*100
  

  
  s$effect_real =(s$Infections[which(s$Scenario=="Baseline")]-s$Infections)/(s$pop*10000)*100
  s$infection_decline = s$Infections[which(s$Scenario=="Baseline")]-s$Infections
  return(s)
}))
write.csv(pd,paste0(outdir,"/fig3b_data.csv"),row.names = F)

S2$citysize = "Mid"
S2$citysize[which(S1$pop<500)] = "Small"
S2$citysize[which(S1$pop>=1000)] = "Large"
al2 = S2%>%group_by(Days,Scenario,citysize,variant)%>%summarise(Mean =sum(Mean),
                                                        CI05 = sum(CI05),
                                                        CI95 = sum(CI95))
pop = S2%>%group_by(citysize,variant)%>%summarise(pop = sum(unique(pop)))
al2 = do.call(rbind,lapply(split(al2,al2$Scenario),function(k){
  f = do.call(rbind,lapply(split(k,k$citysize),function(kc){
    fv = do.call(rbind,lapply(split(kc,kc$variant),function(kcv){
      kcv$cummean = cumsum(kcv$Mean)
      kcv$cumCI05 = cumsum(kcv$CI05)
      kcv$cumCI95 = cumsum(kcv$CI95)
      return(kcv)
    }))
    return(fv)
  }))
  return(f)
}))
al2 = al2[which(al2$Scenario!="withoutCC"),]
pd2 = do.call(rbind,lapply(split(al2,al2$citysize),function(c){
  c = do.call(rbind,lapply(split(c,c$variant),function(s){
    s = do.call(rbind,lapply(split(s,s$Scenario),function(v){
      e = data.frame(Infections = max(v$cummean), 
                     Infections_CI05 = max(v$cumCI05),
                     Infections_CI95 = max(v$cumCI95),
                     Scenario = unique(v$Scenario), 
                     citysize = unique(v$citysize),
                     variant = unique(v$variant),
                     pop = pop$pop[which(pop$citysize==unique(v$citysize))]) 
      return(e)
    }))
    s$effect = (s$Infections[which(s$Scenario=="Baseline")]-s$Infections)/
      unique(s$Infections[which(s$Scenario=="Baseline")])*100
    
    s$effectCI05 =(s$Infections_CI05[which(s$Scenario=="Baseline")]-s$Infections_CI05)/
      (s$Infections_CI05[which(s$Scenario=="Baseline")])*100
    
    s$effectCI95 =(s$Infections_CI95[which(s$Scenario=="Baseline")]-s$Infections_CI95)/
      (s$Infections_CI95[which(s$Scenario=="Baseline")])*100
    return(s)
  }))
  return(c)
}))
pd$citysize = "all"

pda = do.call(rbind,list(pd,pd2))

pda$citysize = factor(pda$citysize,level = c("all","Large","Mid","Small"),
                      labels = c("all"="Overall",
                                 "Large"="Large cities",
                                 "Mid"="Medium cities",
                                 "Small"="Small cities"))
pda = subset(pda,pda$Scenario!="Baseline"&pda$Scenario!="realworld")
pda$Scenario = factor(pda$Scenario,levels = c("withoutPCR","withoutFM","withoutSD"))
pda$variant = factor(pda$variant,levels = c("original&alpha","delta","omicron"))

pd = subset(pd,pd$Scenario!="realworld")
library(ggthemes)
g2 = ggplot(data = pda,aes(x=variant, y=effect, 
                           fill=Scenario, alpha = citysize,color=citysize))+
  geom_bar(position=position_dodge(), stat="identity", # Use black outlines,
           size=.6,width=0.6)+
  geom_errorbar(aes(ymin = effectCI05, ymax= effectCI95), 
                width=0.3,position=position_dodge(.6),alpha=1)+
  scale_y_continuous(#breaks=c(0,1,5,15,30,50),
    limits=c(0,100))+
 # scale_color_manual(values=c("withoutFM" = "#60ACC0",
#                              "withoutSD" = "#257187"
#                              "withoutPCR" = "#8CBCBD"),
#                     labels=c("withoutFM" = "Relative effect of FM",
#                              "withoutSD" = "Relative effect of SD",
#                              "withoutPCR" = "Relative effect of PCR"))+
  scale_fill_manual(values=c("withoutFM" = "#60ACC0",
                             "withoutSD" = "#257187",
                             "withoutPCR" = "#8CBCBD"),
                    labels=c("withoutFM" = "Relative effect of FM",
                             "withoutSD" = "Relative effect of SD",
                             "withoutPCR" = "Relative effect of PCR"))+
  scale_alpha_manual(values=c("all" = 1,
                             "Large" = 0.8,
                             "Mid" = 0.5,
                             "Small"=0.3))+
  theme_bw()+
  labs(x="",y="Relative effect of NPI (%)")+
  theme(legend.position = "right",
        legend.background = element_rect(
          fill = "transparent", # 填充色
          colour = "transparent", # 框线色
          size = 1.5),
        strip.text = element_textbox_highlight(
          size = 10, 
          fill = "transparent", box.color = "transparent", color = "gray10",
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)),
        strip.background.x = element_rect(fill = "white", colour = "transparent"),
        legend.title=element_blank(),
        panel.grid = element_blank(),
        axis.line=element_line(color="black",size=0.5),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size =unit(9,"pt")),
        axis.title.x= element_text(color="black",hjust = 0.5,vjust=0,size =unit(9,"pt")),
        axis.text.x = element_text(color="black",size = unit(9,"pt")),
        axis.title.y= element_text(color="black",hjust = 0.5,vjust=0,size =unit(9,"pt")),
        axis.text.y = element_text(color="black",size = unit(9,"pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))

outdir = "Figures"
ggsave(paste0(outdir,"/fig3c.pdf"),g2, width=16, height=7, units="cm", scale=1)

write.csv(pda,paste0(outdir,"/fig3c_data.csv"),row.names = F)

