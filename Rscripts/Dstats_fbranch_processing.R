
###Libraries
library(tidyverse)


###load data and process
fbranch_full<-read.table("./Data/Dsuite/Brevipes_submission.pass.31calledSNPs_AURA.Fbranch",header=T)
fbranch_long<-fbranch %>% pivot_longer(AURA:JOHN,names_to = "source_pop",values_to="f_branch") %>%
  mutate(source_pop=factor(source_pop,levels=c("JOHN","BREV","FREM","SESP","CONS","NANU"))) %>% 
  mutate(branch_descendants=factor(branch_descendants,levels=c("CONS","SESP","SESP,CONS","FREM","BREV","JOHN","BREV,JOHN","FREM,BREV,JOHN"))) %>%
  filter(!is.na(source_pop),!is.na(f_branch))

### plot full fbranch
fbranch_plot<-ggplot(fbranch_long,aes(x=source_pop,y=branch_descendants,fill=f_branch,label=round(f_branch,3))) + 
  geom_tile(color="black") + geom_text() + 
  scale_fill_gradient(low="white",high="red",limits=c(0,0.10),na.value = "grey") + 
  theme_minimal() + 
  theme(legend.title=element_text(face="italic"),panel.grid=element_blank())

### save plot
png("./Outputs/full_fbranch.png",width=5,height=4,units="in",res=600)
fbranch_plot
dev.off()

###ABBA-BABA and f-4 ratio data

Dstats<-read.table("./Dsuite/Brevipes_submission.pass.31calledSNPs_AURA_tree.txt",header=T) %>%
  mutate(sig1=(p.value<0.05/20),sig2=(p.value<0.01/20),sig3=(p.value<0.001/20))


###Downsampled iterations data
ITER1<-read.table("./downsample/downsample_ITER_1.pass.11calledSNPs.Fbranch",header=T)
ITER2<-read.table("./downsample/downsample_ITER_2.pass.11calledSNPs.Fbranch",header=T)
ITER3<-read.table("./downsample/downsample_ITER_3.pass.11calledSNPs.Fbranch",header=T)
ITER4<-read.table("./downsample/downsample_ITER_4.pass.11calledSNPs.Fbranch",header=T)

###process each iteration
ITER1_long<-ITER1 %>% pivot_longer(AURA:JOHN,names_to = "source_pop",values_to="f_branch") %>%
  mutate(source_pop=factor(source_pop,
                           levels=c("JOHN","BREV","FREM","SESP","CONS","NANU"))) %>% 
  mutate(branch_descendants=factor(branch_descendants,
                                   levels=c("CONS","SESP","SESP,CONS","FREM","BREV","JOHN","BREV,JOHN","FREM,BREV,JOHN"))) %>%
  filter(!is.na(source_pop),!is.na(f_branch))
ITER2_long<-ITER2 %>% pivot_longer(AURA:JOHN,names_to = "source_pop",values_to="f_branch") %>%
  mutate(source_pop=factor(source_pop,
                           levels=c("JOHN","BREV","FREM","SESP","CONS","NANU"))) %>% 
  mutate(branch_descendants=factor(branch_descendants,
                                   levels=c("CONS","SESP","SESP,CONS","FREM","BREV","JOHN","BREV,JOHN","FREM,BREV,JOHN"))) %>%
  filter(!is.na(source_pop),!is.na(f_branch))
ITER3_long<-ITER3 %>% pivot_longer(AURA:JOHN,names_to = "source_pop",values_to="f_branch") %>%
  mutate(source_pop=factor(source_pop,
                           levels=c("JOHN","BREV","FREM","SESP","CONS","NANU"))) %>% 
  mutate(branch_descendants=factor(branch_descendants,
                                   levels=c("CONS","SESP","SESP,CONS","FREM","BREV","JOHN","BREV,JOHN","FREM,BREV,JOHN"))) %>%
  filter(!is.na(source_pop),!is.na(f_branch))
ITER4_long<-ITER4 %>% pivot_longer(AURA:JOHN,names_to = "source_pop",values_to="f_branch") %>%
  mutate(source_pop=factor(source_pop,
                           levels=c("JOHN","BREV","FREM","SESP","CONS","NANU"))) %>% 
  mutate(branch_descendants=factor(branch_descendants,
                                   levels=c("CONS","SESP","SESP,CONS","FREM","BREV","JOHN","BREV,JOHN","FREM,BREV,JOHN"))) %>%
  filter(!is.na(source_pop),!is.na(f_branch))

### combine iterations and plot
iter_fbranch_all<-rbind(ITER1_long,ITER2_long,ITER3_long,ITER4_long) %>% group_by(source_pop,branch_descendants) %>%
  summarize(median_f_branch=median(f_branch),min_f_branch=min(f_branch),max_f_branch=max(f_branch)) %>%
  mutate(isbold=min_f_branch>0)

iter_medians_fbranch<-ggplot(iter_fbranch_all,aes(x=source_pop,y=branch_descendants,fill=median_f_branch)) + 
  geom_tile(color="black") + 
  geom_text(data=subset(iter_fbranch_all,isbold==T),aes(label=round(median_f_branch,3)),fontface="bold") + 
  geom_text(data=subset(iter_fbranch_all,isbold==F),aes(label=round(median_f_branch,3)),fontface="plain") + 
  scale_fill_gradient(low="white",high="red",limits=c(0,0.10),na.value = "grey") + 
  theme_minimal() + 
  theme(legend.title=element_text(face="italic"),panel.grid=element_blank())

iter_maxes_fbranch<-ggplot(iter_fbranch_all,aes(x=source_pop,y=branch_descendants,fill=max_f_branch)) + 
  geom_tile(color="black") + 
  geom_text(data=subset(iter_fbranch_all,isbold==T),aes(label=round(max_f_branch,3)),fontface="bold") + 
  geom_text(data=subset(iter_fbranch_all,isbold==F),aes(label=round(max_f_branch,3)),fontface="plain") + 
  scale_fill_gradient(low="white",high="red",limits=c(0,0.10),na.value = "grey") + 
  theme_minimal() + 
  theme(legend.title=element_text(face="italic"),panel.grid=element_blank())

###save plots
png("Results/Output_graphs/withCONS/iterations_median_fbranch.png",width=7,height=6,units="in",res=600)
iter_medians_fbranch
dev.off()

png("Results/Output_graphs/withCONS/iterations_maxes_fbranch.png",width=7,height=6,units="in",res=600)
iter_maxes_fbranch
dev.off()
