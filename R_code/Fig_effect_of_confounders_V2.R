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

dir =  "Reviewed_IndividualNPI_stanV5_pooled_onebeta"
conlist = list.files(dir,pattern="*.rds")
all = read.csv(paste0(dir,"/testdata_IndNPI.csv"),stringsAsFactors = F)
all = all[which(all$Mean.R.>0),]

f = readRDS(paste0(dir,"/",conlist))
out = rstan::extract(f)

be = as.data.frame(out$beta)
colnames(be)=auxname 
for (i in auxname){
  be[,i] = be[,i]*mean(pd[,i])
  be[,i] = (1-exp(-be[,i]))
}
be2=mcmc_intervals_data(be,prob = .5,prob_outer= .95,point_est="mean")


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


f = ggplot(be2)+
  geom_linerange(aes(xmin =ll, xmax=hh, y=parameter),
                 color= "grey80",
                 size = 3,
                 show.legend = FALSE,alpha=0.3,
                 position=position_dodge(width = 0.7))+
  geom_hline(show.legend = TRUE,aes(yintercept = -3),
             color= "grey80",
             size = 3,alpha=0.3)+
  geom_linerange(aes(xmin =l, xmax=h, y=parameter),
                 color= "grey80",
                 size = 3,show.legend = FALSE,
                 position=position_dodge(width = 0.8))+
  geom_pointrange(aes(x = m, y = parameter, xmin=m,xmax=m),
                  color= "grey20",
                  fatten = 8, show.legend = TRUE,alpha=0.5,
                  position= position_dodge(width = 0.7))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw()+
  scale_x_continuous(expand=c(0.01,0.01),limits=range(-0.20,0.40),labels=percent)+
  labs(x =expression("Reduction in "*R[t]),y = NULL)+
  theme(legend.position = "top",legend.key=element_rect(fill='transparent'),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = "transparent"),
        panel.grid.minor = element_line(colour = "transparent"),
        legend.title=element_blank(),
        plot.margin = margin(t = 5,  # ¶¥²¿±ßÔµ¾àÀë
                             r = 5,  # ÓÒ±ß±ßÔµ¾àÀë
                             b = 5,  # µ×²¿±ßÔµ¾àÀë
                             l = 5),
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
ggsave(paste0(outdir,"/figConfounder.pdf"),f, width=16, height=6, units="cm", scale=1)

write.csv(be2,paste0(outdir,"/figConfounder.csv"),row.names = F)



