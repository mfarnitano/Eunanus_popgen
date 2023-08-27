
### Libraries
library(tidyverse)
library(cowplot)
library(zoo)

### Load d_f data from Dinvestigate

dfs_CSFO<-read.table("./Data/Dinvestigate/CSFA_Dinvestigate_31calledSNPs.txt",header=T) %>%
  mutate(midpoint=floor((windowStart+windowEnd)/2),whichwindow=floor(midpoint/10000)*10000+1) 
dfs_BCNO<-read.table("./Data/Dinvestigate/BCNA_Dinvestigate_31calledSNPs.txt",header=T) %>%
  mutate(midpoint=floor((windowStart+windowEnd)/2),whichwindow=floor(midpoint/10000)*10000+1) 



### Load pi values in 10Kb windows from pixy

moresites_pi_10Kb<-read.table("./Data/pixy/Brevipes_submission.pass.31called.10000.bypop_pi.txt",header=T)
avg_pi_CONS_each<-moresites_pi_10Kb %>% filter(str_detect(pop,"CONS")) %>% select(chromosome,window_pos_1,avg_pi)
avg_pi_SESP_each<-moresites_pi_10Kb %>% filter(str_detect(pop,"SESP")) %>% select(chromosome,window_pos_1,avg_pi)
avg_pi_NANU_each<-moresites_pi_10Kb %>% filter(str_detect(pop,"NANU")) %>% select(chromosome,window_pos_1,avg_pi)
avg_pi_JOHN_each<-moresites_pi_10Kb %>% filter(str_detect(pop,"JOHN")) %>% select(chromosome,window_pos_1,avg_pi)
avg_pi_FREM_each<-moresites_pi_10Kb %>% filter(str_detect(pop,"FREM")) %>% select(chromosome,window_pos_1,avg_pi) %>%
  group_by(chromosome,window_pos_1) %>% summarize(avg_pi=mean(avg_pi,na.rm=T))
avg_pi_BREV_each<-moresites_pi_10Kb %>% filter(str_detect(pop,"BREV")) %>% select(chromosome,window_pos_1,avg_pi) %>%
  group_by(chromosome,window_pos_1) %>% summarize(avg_pi=mean(avg_pi,na.rm=T))


### merge d_f and pi data and filter for pi>0

dfs_CSFO_unfiltered<-dfs_CSFO %>% left_join(avg_pi_CONS_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  rename("avg_pi_CONS"="avg_pi") %>% 
  left_join(avg_pi_SESP_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  rename("avg_pi_SESP"="avg_pi") %>% 
  left_join(avg_pi_FREM_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  rename("avg_pi_FREM"="avg_pi")

dfs_CSFO_pifiltered<-dfs_CSFO %>% left_join(avg_pi_CONS_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  filter(!is.na(avg_pi),avg_pi>0)  %>% rename("avg_pi_CONS"="avg_pi") %>% 
  left_join(avg_pi_SESP_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  filter(!is.na(avg_pi),avg_pi>0)  %>% rename("avg_pi_SESP"="avg_pi") %>% 
  left_join(avg_pi_FREM_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  filter(!is.na(avg_pi),avg_pi>0)  %>% rename("avg_pi_FREM"="avg_pi") 

dfs_BCNO_unfiltered<-dfs_BCNO %>% left_join(avg_pi_CONS_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  rename("avg_pi_CONS"="avg_pi") %>% 
  left_join(avg_pi_NANU_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  rename("avg_pi_NANU"="avg_pi") %>% 
  left_join(avg_pi_BREV_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  rename("avg_pi_BREV"="avg_pi")

dfs_BCNO_pifiltered<-dfs_BCNO %>% left_join(avg_pi_CONS_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  filter(!is.na(avg_pi),avg_pi>0)  %>% rename("avg_pi_CONS"="avg_pi") %>% 
  left_join(avg_pi_NANU_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  filter(!is.na(avg_pi),avg_pi>0)  %>% rename("avg_pi_NANU"="avg_pi") %>% 
  left_join(avg_pi_BREV_each,by=c("chr"="chromosome","whichwindow"="window_pos_1")) %>% 
  filter(!is.na(avg_pi),avg_pi>0)  %>% rename("avg_pi_BREV"="avg_pi") 



### Calculate cutoff values to find outliers

CSFO_cutoff_pifiltered<-abs((dfs_CSFO_pifiltered %>% arrange(desc(abs(d_f))) %>% head(100) %>% tail(1))$d_f)
BCNO_cutoff_pifiltered<-abs((dfs_BCNO_pifiltered %>% arrange(desc(abs(d_f))) %>% head(100) %>% tail(1))$d_f)



### Get genome-wide average values

dfs_CSFO_unfiltered %>% summarize(mean(D),mean(d_f),median(D),median(d_f),mean(avg_pi_CONS,na.rm=T),mean(avg_pi_SESP,na.rm=T),mean(avg_pi_FREM,na.rm=T))
dfs_CSFO_pifiltered %>% summarize(mean(D),mean(d_f),median(D),median(d_f),mean(avg_pi_CONS),mean(avg_pi_SESP),mean(avg_pi_FREM))

dfs_BCNO_unfiltered %>% summarize(mean(D),mean(d_f),median(D),median(d_f),mean(avg_pi_BREV,na.rm=T),mean(avg_pi_CONS,na.rm=T),mean(avg_pi_NANU,na.rm=T))
dfs_BCNO_pifiltered %>% summarize(mean(D),mean(d_f),median(D),median(d_f),mean(avg_pi_BREV),mean(avg_pi_CONS),mean(avg_pi_NANU))

### summaries for outlier values
dfs_CSFO_pifiltered %>% filter(d_f>CSFO_cutoff_pifiltered) %>% 
  summarize(mean(D),mean(d_f),median(D),median(d_f),mean(avg_pi_CONS,na.rm=T),mean(avg_pi_SESP,na.rm=T),mean(avg_pi_FREM,na.rm=T))
dfs_BCNO_pifiltered %>% filter(d_f>=BCNO_cutoff_pifiltered) %>% 
  summarize(mean(D),mean(d_f),median(D),median(d_f),mean(avg_pi_BREV,na.rm=T),mean(avg_pi_CONS,na.rm=T),mean(avg_pi_NANU,na.rm=T))

### mean and median of pi values

dfs_CSFO_pifiltered %>% summarize(across(avg_pi_CONS:avg_pi_FREM,list("mean"=mean,"median"=median)))
dfs_CSFO_pifiltered %>% filter(d_f>=CSFO_cutoff_pifiltered) %>% summarize(across(avg_pi_CONS:avg_pi_FREM,list("mean"=mean,"median"=median)))

dfs_BCNO_pifiltered %>% summarize(across(avg_pi_CONS:avg_pi_BREV,list("mean"=mean,"median"=median)))
dfs_BCNO_pifiltered %>% filter(d_f>=BCNO_cutoff_pifiltered) %>% summarize(across(avg_pi_CONS:avg_pi_BREV,list("mean"=mean,"median"=median)))


### Examine distribution of pi values

h1<-ggplot(avg_pi_CONS_each,aes(x=avg_pi)) + stat_bin(binwidth=0.005,closed=c("right"),geom="step",aes(color="CONS")) + 
  stat_bin(data=avg_pi_FREM_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="FREM")) + 
  stat_bin(data=avg_pi_SESP_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="SESP")) + 
  stat_bin(data=avg_pi_NANU_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="NANU")) + 
  stat_bin(data=avg_pi_BREV_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="BREV")) + 
  stat_bin(data=avg_pi_JOHN_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="JOHN")) + 
  theme_bw() + scale_color_manual(values=c("#CE80C6","#321540","#7B0D1E","#4561BA","#CABA27","#659A61"),breaks = c("FREM","CONS","SESP","NANU","BREV","JOHN"),name="pi",guide="none")

h2<-ggplot(avg_pi_CONS_each,aes(x=avg_pi)) + stat_bin(binwidth=0.005,closed=c("right"),geom="step",aes(color="CONS")) + 
  stat_bin(data=avg_pi_FREM_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="FREM")) + 
  stat_bin(data=avg_pi_SESP_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="SESP")) + 
  stat_bin(data=avg_pi_NANU_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="NANU")) + 
  stat_bin(data=avg_pi_BREV_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="BREV")) + 
  stat_bin(data=avg_pi_JOHN_each,binwidth=0.005,closed=c("right"),geom="step",aes(color="JOHN")) + 
  theme_bw() + xlim(-0.005,0.05) + scale_color_manual(values=c("#CE80C6","#321540","#7B0D1E","#4561BA","#CABA27","#659A61"),breaks = c("FREM","CONS","SESP","NANU","BREV","JOHN"),name="pi",guide="none")

h3<-ggplot(avg_pi_CONS_each,aes(x=avg_pi)) + stat_bin(binwidth=0.001,closed=c("right"),geom="step",aes(color="CONS")) + 
  stat_bin(data=avg_pi_FREM_each,binwidth=0.001,closed=c("right"),geom="step",aes(color="FREM")) + 
  stat_bin(data=avg_pi_SESP_each,binwidth=0.001,closed=c("right"),geom="step",aes(color="SESP")) + 
  stat_bin(data=avg_pi_NANU_each,binwidth=0.001,closed=c("right"),geom="step",aes(color="NANU")) + 
  stat_bin(data=avg_pi_BREV_each,binwidth=0.001,closed=c("right"),geom="step",aes(color="BREV")) + 
  stat_bin(data=avg_pi_JOHN_each,binwidth=0.001,closed=c("right"),geom="step",aes(color="JOHN")) + 
  theme_bw() + xlim(-0.001,0.02) + scale_color_manual(values=c("#CE80C6","#321540","#7B0D1E","#4561BA","#CABA27","#659A61"),breaks = c("FREM","CONS","SESP","NANU","BREV","JOHN"),name="pi",guide="none")

h_legend<-ggdraw(cowplot::get_legend(h3+scale_color_manual(values=c("#CE80C6","#321540","#7B0D1E","#4561BA","#CABA27","#659A61"),breaks = c("FREM","CONS","SESP","NANU","BREV","JOHN"),name="pi")))

plot_grid(h1,h2,h3,h_legend,align = "hv",axis="lrtb",nrow = 4,rel_widths = c(3,3,3,1))



### Function to plot d_f traces along the genome

Dinvestigate_plotting<-function (df_frame,name,roll=F) {
  qdata<-df_frame %>%
    mutate(chr=factor(chr,levels=c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10")))
  qhistogram<-ggplot(qdata,aes(x=d_f)) + geom_histogram(binwidth=0.05,boundary=0) + geom_vline(xintercept=0) + xlim(-1,1) + 
    theme_bw() + theme(text=element_text(size=12))
  qlastwindows=qdata %>% group_by(chr) %>% summarize(lastx=max(windowEnd))
  qchoose<-qdata %>% mutate(absdf=abs(d_f)) %>% arrange(desc(absdf)) %>% head(100) %>% mutate(isRed=(d_f>0))
  qchoose_cutoff<-min(qchoose$absdf)
  qcolored<-qdata %>% mutate(cat=ifelse(d_f>=qchoose_cutoff,"red",ifelse(d_f<=-qchoose_cutoff,"black","grey")))
  qmin<-min(qchoose$d_f)
  fpr<-qchoose %>% group_by(isRed) %>% summarize(count=n())
  plot100<-ggplot(qcolored,aes(x=windowStart,y=d_f,color=cat)) + 
    geom_hline(yintercept=qmin,color="grey",linetype="dotted") + 
    geom_hline(yintercept=-qmin,color="grey",linetype="dotted") + 
    geom_point(size=0.5) + facet_wrap(~chr,nrow=1) + 
    geom_vline(data=qlastwindows,aes(xintercept=lastx,color=NULL),color="black") + 
    scale_y_continuous(limits=c(-1,1),breaks=c(-1,0,1)) + 
    scale_x_continuous(breaks=c(0,1e7,2e7),labels=function(x){x/(1e6)}) + 
    scale_color_manual(values=c("black","red","grey"),limits=c("black","red","grey"),guide="none") + 
    theme_bw()
  qratio<-fpr$count[1]/fpr$count[2]
  if (roll==T) {
    qrolled<-qcolored %>% group_by(chr) %>% mutate(rolled=rollapply(d_f,width=21,FUN=mean,fill=NA))
    plot100<-ggplot(qrolled,aes(x=windowStart,y=d_f,color=cat)) + 
      geom_hline(yintercept=qmin,color="grey",linetype="dotted") + 
      geom_hline(yintercept=-qmin,color="grey",linetype="dotted") + 
      geom_point(size=0.5) + facet_wrap(~chr,nrow=1) + 
      geom_hline(yintercept=0,color="black",linetype="dashed") + 
      geom_line(aes(x=windowStart,y=rolled,color=NULL),color="blue") + 
      geom_vline(data=qlastwindows,aes(xintercept=lastx,color=NULL),color="black") + 
      scale_y_continuous(limits=c(-1,1),breaks=c(-1,0,1)) + 
      scale_x_continuous(breaks=c(0,1e7,2e7),labels=function(x){x/(1e6)}) + 
      scale_color_manual(values=c("black","red","grey"),limits=c("black","red","grey"),guide="none") + 
      theme_bw()
  }
  
  return(list(qhistogram=qhistogram,plot100=plot100,qratio=qratio))
}


### Plot traces for CSFO

CSFO_plotting<-Dinvestigate_plotting(dfs_CSFO_pifiltered,"CSFO_pifiltered")
CSFO_roller<-Dinvestigate_plotting(dfs_CSFO_pifiltered,"CSFO_pifiltered",roll=T)

png(paste0("./Outputs/CSFO_pifiltered_histogram.png"),width=3,height=3,res=600,units="in")
CSFO_plotting$qhistogram
dev.off()

png(paste0("./Outputs/CSFO_pifiltered_dfplot_allwindows.png"),width=6,height=2,res=600,units="in")
CSFO_plotting$plot100
dev.off()

png(paste0("./Outputs/CSFO_pifiltered_rolled_dfplot_allwindows.png"),width=6,height=2,res=600,units="in")
CSFO_roller$plot100
dev.off()


### Plot traces for BCNO

BCNO_plotting<-Dinvestigate_plotting(dfs_BCNO_pifiltered,"BCNO_pifiltered")
BCNO_roller<-Dinvestigate_plotting(dfs_BCNO_pifiltered,"BCNO_pifiltered",roll=T)
png(paste0("./Outputs/BCNO_pifiltered_rolled_dfplot_allwindows.png"),width=6,height=2,res=600,units="in")
BCNO_roller$plot100
dev.off()



### Divide into deciles of pi and get mean d_f


SESP_pi_quantiles<-quantile(dfs_CSFO_pifiltered$avg_pi_SESP,c(1:10)/10)
CONS_pi_quantiles<-quantile(dfs_CSFO_pifiltered$avg_pi_CONS,c(1:10)/10)
FREM_pi_quantiles<-quantile(dfs_CSFO_pifiltered$avg_pi_FREM,c(1:10)/10)
pick_quantile<-function(value,quantiles) {
  if (value>max(quantiles)) {return(NA)}
  return(match(min(quantiles[quantiles>=value]),quantiles))
}
dfs_CSFO_pifiltered %>% mutate(SESP_pi_quantile=sapply(avg_pi_SESP,pick_quantile,quantiles=SESP_pi_quantiles)) %>%
  group_by(SESP_pi_quantile) %>% summarize(mean(d_f))
dfs_CSFO_pifiltered %>% mutate(CONS_pi_quantile=sapply(avg_pi_CONS,pick_quantile,quantiles=CONS_pi_quantiles)) %>%
  group_by(CONS_pi_quantile) %>% summarize(mean(d_f))
dfs_CSFO_pifiltered %>% mutate(FREM_pi_quantile=sapply(avg_pi_FREM,pick_quantile,quantiles=FREM_pi_quantiles)) %>%
  group_by(FREM_pi_quantile) %>% summarize(mean(d_f))

BREV_pi_quantiles<-quantile(dfs_BCNO_pifiltered$avg_pi_BREV,c(1:10)/10)
CONS_pi_quantiles_BCNO<-quantile(dfs_BCNO_pifiltered$avg_pi_CONS,c(1:10)/10)
NANU_pi_quantiles<-quantile(dfs_BCNO_pifiltered$avg_pi_NANU,c(1:10)/10)

dfs_BCNO_pifiltered %>% mutate(BREV_pi_quantile=sapply(avg_pi_BREV,pick_quantile,quantiles=BREV_pi_quantiles)) %>%
  group_by(BREV_pi_quantile) %>% summarize(mean(d_f))
dfs_BCNO_pifiltered %>% mutate(CONS_pi_quantile=sapply(avg_pi_CONS,pick_quantile,quantiles=CONS_pi_quantiles_BCNO)) %>%
  group_by(CONS_pi_quantile) %>% summarize(mean(d_f))
dfs_BCNO_pifiltered %>% mutate(NANU_pi_quantile=sapply(avg_pi_NANU,pick_quantile,quantiles=NANU_pi_quantiles)) %>%
  group_by(NANU_pi_quantile) %>% summarize(mean(d_f))

