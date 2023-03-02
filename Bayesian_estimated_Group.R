library("rstan")
library("bayesplot")
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
source('code/NPI_code.R')

indir = "1220dataset"
##############Group Distinction#########
outdir = "0228Version_IndividualNPI_GroupwithoutVariant"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}
labelname = "Label_q"

all = read.csv(paste0(indir,"/Rt&smooth_NPI.csv"),stringsAsFactors = F)
all = all[which(all$Mean.R.>0),]
all$Control[which(all$Control>1)] = 1
NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")
auxname = c("Temp","fully_vaccination_rate")
mcon=stan_model('stanmodel/stanmodel_V3.stan')
outfile = data.frame()
for(c1 in split(all,all[,labelname])){
  c1[,"Temp"] = (c1[,"Temp"]-min(c1[,"Temp"]))/(max(c1[,"Temp"])-min(c1[,"Temp"]))
  testdata = domcmc_group(mcon, c1, NPIname, auxname, outdir, labelname)
  outfile  = rbind(outfile,testdata)
}
write.csv(outfile,paste0(outdir,"/testdata_IndNPI.csv"),row.names = F)


library("rstan")
library("bayesplot")
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")

NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")
auxname = c("Temp","Fully_vaccination")
ObNPI = c("Static.management",
          "Facial_Mask",
          "PCR.screening",
          "Contact.control") 

conlist = list.files(outdir,pattern="*.rds")
all = read.csv(paste0(outdir,"/testdata_IndNPI.csv"),stringsAsFactors = F)
all = all[which(all$Mean.R.>0),]


al2 = do.call(rbind,lapply(seq(1,length(conlist)),function(j){
  conf = conlist[j]
  f = readRDS(paste0(outdir,"/",conf))
  group = strsplit(conf,"_")[[1]][1]
  v = all[which(all$Label_q==substring(group,6,6)),]
  
  out = rstan::extract(f)
  al=as.data.frame(out$alpha)
  colnames(al)=NPIname
  for (i in NPIname){
    al[,i] = al[,i]*mean(v[,i])
  }
  
  al$Static.management = al$Lockdown+al$Business_Premises_Closure+al$Public_Transportation_Closure+
    al$Gathering_restriction+al$Workplace_Closure+al$School_Closure
  al$PCR.screening = al$Medicine_Management+al$Mass_screening
  colnames(al)[which(colnames(al)=="Control")] = "Contact.control"
  al= al[,ObNPI]
  for (i in ObNPI){
    al[,i] = (1-exp(-al[,i]))
  }
  parv=mcmc_intervals_data(al,prob = .5,prob_outer= .95,point_est="median")
  parv$Group = paste("Group",unique(v$Label_q))
  return(parv)
}))

al2$Group = factor(al2$Group,levels=c("Group 4","Group 3","Group 2","Group 1"))
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
  geom_linerange(aes(xmin =ll, xmax=hh, y=Group,color= Group),
                 show.legend = FALSE,alpha=0.3,size=1,
                 position=position_dodge(width = 0.5))+
  geom_hline(show.legend = TRUE,aes(yintercept = -3, color= Group),alpha=0.3)+
  geom_linerange(aes(xmin =l, xmax=h, y=Group,color= Group),show.legend = FALSE,
                 position=position_dodge(width = 0.5),size=1.5)+
  geom_pointrange(aes(x = m, y = Group, xmin=m,xmax=m,
                      color =  Group), fatten = 5, show.legend = TRUE,alpha=0.5,
                  position= position_dodge(width = 0.5))+
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
  scale_x_continuous(expand=c(0.01,0.01),limits=range(-0.01,1.00),labels=percent)+
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

outdir = "FiguresV1129"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}
ggsave(paste0(outdir,"/figs7_NPIeffect_group.pdf"),f, width=16, height=16, units="cm", scale=1)

write.csv(al2,paste0(outdir,"/figs7_NPIeffect_group_data.csv"),row.names = F)

