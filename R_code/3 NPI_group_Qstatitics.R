library(dplyr)
library(geodetector)
library(ggplot2)
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")

strata_data<-function (group,dir){
  #pdf(paste0(dir,"/Starta.pdf"))
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
dir = "pre-dataset"
group = read.csv(paste0(dir,"/pandemic_shape.csv"),stringsAsFactors = F)

groupnew = strata_data(group,dir)
write.csv(groupnew,paste0(dir,"/alldata_group_Qstatistics.csv"),row.names = F)