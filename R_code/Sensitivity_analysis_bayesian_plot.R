#sensitivity analysis : Bayesian#


library("rstan")
library("bayesplot")
setwd("E:/zeroCOVID_NPI/Version0504")

alldir = "Reviewed_SE_Bayesian"

data = read.csv(paste0(alldir,"/SEdataset.csv"),stringsAsFactors = F)



conlist = list.files(alldir,pattern="*.rds")
NPIname = c("Lockdown","Business_Premises_Closure",
            "Public_Transportation_Closure","Gathering_restriction",
            "Workplace_Closure","School_Closure",
            "Medicine_Management","Mass_screening",
            "Facial_Mask","Control")
ObNPI = c("Static.management",
          "Facial_Mask",
          "PCR.screening",
          "Contact.control") 

par = do.call(rbind,lapply(conlist,function(se){
  fit = readRDS(paste0(alldir,"/",se))
  out = rstan::extract(fit)
  al=as.data.frame(out$alpha)
  colnames(al)=NPIname
  for (i in NPIname){
    al[,i] = al[,i]*mean(data[,i])
  }
  al$Static.management = al$Lockdown+al$Business_Premises_Closure+al$Public_Transportation_Closure+
    al$Gathering_restriction+al$Workplace_Closure+al$School_Closure
  al$PCR.screening = al$Medicine_Management+al$Mass_screening
  colnames(al)[which(colnames(al)=="Control")] = "Contact.control"
  al= al[,ObNPI]
  for (i in ObNPI){
    al[,i] = (1-exp(-al[,i]))
  }
  all=mcmc_intervals_data(al,prob = .5,prob_outer= .95,point_est="mean")
  all$SE = substring(se,1,3)
  return(all)
}))

alldir = "Reviewed_overallNPI_stanV5_addpop"

pd = read.csv(paste0(alldir,"/testdata_IndNPI.csv"),stringsAsFactors = F)
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
for (i in ObNPI){
  al[,i] = (1-exp(-al[,i]))
}
al1=mcmc_intervals_data(al,prob = .5,prob_outer= .95,point_est="median")
al1$SE = "DBS"
SEdata = do.call(rbind,list(par,al1))
SEdata$SE = factor(SEdata$SE,levels = c("BS4","BS3","BS2","BS1","DBS"))

fig = ggplot(SEdata)+
  geom_linerange(mapping = aes_(xmin =~ll, xmax=~hh, y=~parameter,color= ~SE),show.legend = FALSE,alpha=0.4,
                 size = 0.4, position=position_dodge(width = 0.7))+
  geom_hline(show.legend = TRUE,aes(yintercept = -1, color=SE),size = 0.5,alpha=0.4)+
  geom_linerange(mapping = aes_(xmin =~l, xmax=~h, y=~parameter,color=~SE),show.legend = FALSE,
                 size =1,position=position_dodge(width = 0.7))+
  geom_point(mapping = aes_(x = ~m, y = ~parameter, color = ~SE),show.legend = TRUE,alpha=0.7,
             size =2,position= position_dodge(width = 0.7))+
  scale_color_manual(values=c("#87a0b7","#8CBCBD","#800026","#cf221f","black"))+
  theme_bw()+
  scale_x_continuous(limits=range(-0.1,1))+
  labs(x =NULL,y = NULL)+
  theme(legend.position = "top",
        panel.grid.major = element_line(colour = "transparent"),
        panel.grid.minor = element_line(colour = "transparent"),
        legend.title=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(9, "pt")),
        axis.title.x= element_blank(),
        axis.text.x = element_text(color="black",hjust=0.5,size = unit(9, "pt")),
        axis.title.y= element_text(color="black",size = unit(9, "pt")),
        axis.text.y= element_text(color="black",hjust=1,size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
outdir = "Figures_Reviewed"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}
ggsave(paste0(outdir,"/figSE_bayesian.pdf"),fig, width=14,
       height=9, units="cm", scale=1)
