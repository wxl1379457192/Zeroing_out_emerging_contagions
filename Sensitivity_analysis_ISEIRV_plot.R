#sensitivity analysis : ISEIRV#
library(ggplot2)
library(openxlsx)
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")

alldir = "SEIR_1220_SE"

parlist= list.files(alldir,pattern="*.xlsx")

k1=parlist[[1]]

se = do.call(rbind,lapply(parlist,function(k1){
  f = read.xlsx(paste0(alldir,"/",k1), sheet = 1)
  par = do.call(rbind,lapply(seq(1,8),function(i){
      data.frame(par = colnames(f)[i], ll = quantile(f[,i])[1],
                 l = quantile(f[,i])[2],
                 m = quantile(f[,i])[3],
                 h =quantile(f[,i])[4],
                 hh = quantile(f[,i])[5])
      }))
  par$SE = substring(k1,1,3)
  return(par)
}))
se$SE = factor(se$SE,levels = c("IS3","IS2","IS1","DIS"))

fig = ggplot(se)+
  geom_linerange(mapping = aes_(xmin =~ll, xmax=~hh, y=~par,color= ~SE),show.legend = FALSE,alpha=0.4,
                 size = 0.4, position=position_dodge(width = 0.75))+
  geom_hline(show.legend = TRUE,aes(yintercept = -1, color=SE),size = 0.5,alpha=0.4)+
  geom_linerange(mapping = aes_(xmin =~l, xmax=~h, y=~par,color=~SE),show.legend = FALSE,
                 size =1,position=position_dodge(width = 0.75))+
  geom_point(mapping = aes_(x = ~m, y = ~par, color = ~SE),show.legend = TRUE,alpha=0.7,
             size =2,position= position_dodge(width = 0.75))+
  scale_color_manual(values=c("#3D729A","#5E86C1","#708090","black"))+
  theme_bw()+
  #scale_x_continuous(limits=range(0,90))+
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
outdir = "FiguresV1129"
if (dir.exists(outdir)==F){dir.create(outdir)
}else{print("This path has been exists")}
ggsave(paste0(outdir,"/figSE_SEIR.pdf"),fig, width=14,
       height=16, units="cm", scale=1)