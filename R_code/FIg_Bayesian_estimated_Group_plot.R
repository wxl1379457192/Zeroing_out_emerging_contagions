library("rstan")
library("bayesplot")
library(dplyr)
setwd("E:/zeroCOVID_NPI/Version0504")
source('code/NPI_code.R')

indir = "dataset"
##############Group Distinction#########
outdir = "Reviewed_Bayesian_model_group"
labelname = "Label_q"
NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")
auxname = c("Temp","Fully_vaccination","pop_density")
ObNPI = c("Static.management",
          "Facial_Mask",
          "PCR.screening",
          "Contact.control") 

all = read.csv(paste0(outdir,"/testdata_IndNPI.csv"),stringsAsFactors = F)
all = all[which(all$Mean.R.>0),]
al2 = do.call(rbind,lapply(unique(all$Label_q),function(g){
  allg = split(all,all$Label_q)[[g]]
  conlist = paste0("Pooled_model_Group",g,".rds")
  f = readRDS(paste0(outdir,"/",conlist))
  out = rstan::extract(f)
  al=as.data.frame(out$alpha)
  if(g ==3){
    variant = c("original&alpha","delta")
  }else{
    variant = c("original&alpha","delta","omicron")
  }
  
  
  syn2 = do.call(rbind,lapply(seq(1,length(variant)),function(x){
    v = allg[which(allg$VG==variant[x]),]
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
    synv=mcmc_intervals_data(alv,prob = .5,prob_outer= .95,point_est="mean")
    
    synv$variant = unique(v$VG)
    synv$Group = paste("Group",g)
    return(synv)
  }))
  print(paste(g,"has been processed"))
  return(syn2)
}))


al2$Group = factor(al2$Group,levels=c("Group 4","Group 3","Group 2","Group 1"))
al2$variant = factor(al2$variant,levels = c("omicron","delta","original&alpha","Overall"))

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
f = ggplot(al2)+
  geom_linerange(aes(xmin =ll, xmax=hh, y=variant,color= Group),
                 show.legend = FALSE,alpha=0.3,size=0.3,
                 position=position_dodge(width = 0.6))+
  geom_hline(show.legend = TRUE,aes(yintercept = -3, color= Group),alpha=0.3)+
  geom_linerange(aes(xmin =l, xmax=h, y=variant,color= Group),show.legend = FALSE,
                 position=position_dodge(width = 0.6),size=0.8)+
  geom_pointrange(aes(x = m, y = variant, xmin=m,xmax=m,
                      color =  Group), fatten = 3, show.legend = TRUE,alpha=0.5,
                  position= position_dodge(width = 0.6))+
  facet_grid(parameter~.)+ 
  scale_color_manual(values=c("#9b9cc9",
                              "#88A9B9",
                              "#66C2a5",
                              "#ACD851"
                              #"#3773a4"#"#6b9ec7",
  ))+theme_grey()+
  scale_fill_manual(values=c("#9b9cc9",
                             "#88A9B9",
                             "#66C2a5",
                             "#ACD851"))+
  scale_x_continuous(expand=c(0.01,0.01),limits=range(-0.05,1.00),labels=percent)+
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

outdir = "Figures_Reviewed"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}
ggsave(paste0(outdir,"/figs7_NPIeffect_group.pdf"),f, width=16, height=16, units="cm", scale=1)

write.csv(al2,paste0(outdir,"/figs7_NPIeffect_group_data.csv"),row.names = F)

