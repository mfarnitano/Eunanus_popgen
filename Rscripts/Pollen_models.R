
### Libraries
library(tidyverse)
library(beeswarm)
library(lme4)
library(emmeans)
library(multcomp)
library(multcompView)


###load data
pollen_counts<-read.csv("./Data/Pollen_hemocytometer_counts.csv")
with_backcrosses<-filter(pollen_counts,Group %in% c("Brevipes","Fremontii","BxF","BFxF","BFxB","BxBF")) %>%
  mutate(Proportion_viable=n_viable/(n_inviable+n_viable),
         Pollen_per_square=(n_viable+n_inviable)/Squares_counted,
         Pollen_per_flower=Pollen_per_square*500,
         Group=factor(Group,levels=c("Brevipes","Fremontii","BxF","BxBF","BFxB","BFxF")),
         Population_Line=factor(Population_Line),
         Individual=factor(Individual))
by_individual<-with_backcrosses %>% group_by(Group,Population_Line,Individual) %>%
  summarize(Countflowers=n(),
            n_viable=sum(n_viable),
            n_inviable=sum(n_inviable),
            Proportion_viable=mean(Proportion_viable),
            Pollen_per_flower=mean(Pollen_per_flower),
            Pollen_per_square=mean(Pollen_per_square)) %>%
  mutate(color2=cbind(Group,"_color"))
sample_size_summary<-by_individual %>% group_by(Group,Population_Line) %>% summarize(samplesize=n())

### Pollen viability model
exclude_small<-filter(by_individual,n_viable+n_inviable>=10)
pv_model<-glmer(cbind(n_viable,n_inviable) ~ Group + (1|Population_Line),data=exclude_small,family=binomial("logit"))
pv_emmeans<-emmeans(pv_model,"Group",type="response")
pv_cld<-cld(pv_emmeans,Letters="abcd")

pv_plot<-ggplot(exclude_small,aes(x=Group,y=n_viable/(n_inviable+n_viable),color=Group,fill=Group)) +
  geom_boxplot(outlier.shape=NA,fill=NA) + geom_beeswarm(cex=1.5,shape=21) +
  geom_text(data=pv_cld,aes(label=.group),y=1.1,color="black") +
  ylab("Viability proportion of counted pollen grains") + xlab(NULL) +
  scale_y_continuous(limits=c(0,1.1),breaks=c(0:4)/4) +
  scale_color_manual(limits=c("Brevipes","Fremontii","BxF","BxBF","BFxB","BFxF"),values=c("#CABA27","#CE80C6","black","grey","grey","grey"),guide="none") +
  scale_fill_manual(limits=c("Brevipes","Fremontii","BxF","BxBF","BFxB","BFxF"),values=c("#CABA27","#CE80C6","black","#CABA27","#CABA27","#CE80C6"),guide="none") +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust = 1)) + theme(text=element_text(size=12))
pv_plot

### Pollen counts model
pc_model<-glmer(round(Pollen_per_square) ~ Group + (1|Population_Line),data=by_individual,family=poisson(link="log"))
pc_emmeans<-emmeans(pc_model,"Group")
pc_cld<-cld(pc_emmeans,Letters="abcd")

pc_plot<-ggplot(by_individual,aes(x=Group,y=round(Pollen_per_square),color=Group,fill=Group)) +
  geom_boxplot(outlier.shape=NA,fill=NA) + geom_beeswarm(cex=1.5,shape=21) +
  geom_text(data=pc_cld,aes(label=.group),y=150,color="black") + ylim(0,150) +
  ylab("Total pollen grains (viable+inviable) per square") + xlab(NULL) +
  scale_color_manual(limits=c("Brevipes","Fremontii","BxF","BxBF","BFxB","BFxF"),
                     values=c("#CABA27","#CE80C6","black","grey","grey","grey"),guide="none") +
  scale_fill_manual(limits=c("Brevipes","Fremontii","BxF","BxBF","BFxB","BFxF"),
                    values=c("#CABA27","#CE80C6","black","#CABA27","#CABA27","#CE80C6"),guide="none") +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust = 1)) + theme(text=element_text(size=12))
pc_plot

### combine plots and save
both<-align_plots(pv_plot,pc_plot,align="hv",axis="tblr")
png("./Outputs/pollen_viability.png",width=6,height=6,units="in",res=600)
plot_grid(plotlist=both,nrow=1)
dev.off()


###Parental, F1 data only, without backcrosses
no_backcrosses<-filter(pollen_counts,Group %in% c("Brevipes","Fremontii","BxF")) %>%
  mutate(Proportion_viable=n_viable/(n_inviable+n_viable),
         Pollen_per_square=(n_viable+n_inviable)/Squares_counted,
         Pollen_per_flower=Pollen_per_square*500,
         Group=factor(Group,levels=c("Brevipes","Fremontii","BxF")),
         Population_Line=factor(Population_Line),
         Individual=factor(Individual))
nobc_individual<-no_backcrosses %>% group_by(Group,Population_Line,Individual) %>%
  summarize(Countflowers=n(),
            n_viable=sum(n_viable),
            n_inviable=sum(n_inviable),
            Proportion_viable=mean(Proportion_viable),
            Pollen_per_flower=mean(Pollen_per_flower),
            Pollen_per_square=mean(Pollen_per_square)) %>%
  mutate(color2=cbind(Group,"_color"))
nobc_sample_size_summary<-nobc_individual %>% group_by(Group,Population_Line) %>% summarize(samplesize=n())

###Pollen viability model, no backcrosses
nobc_exclude_small<-filter(nobc_individual,n_viable+n_inviable>=10)
nobc_pv_model<-glmer(cbind(n_viable,n_inviable) ~ Group + (1|Population_Line),data=nobc_exclude_small,family=binomial("logit"))
nobc_pv_emmeans<-emmeans(nobc_pv_model,"Group",type="response")
nobc_pv_cld<-cld(nobc_pv_emmeans,Letters="abcd")

nobc_pv_plot<-ggplot(nobc_exclude_small,aes(x=Group,y=n_viable/(n_inviable+n_viable),color=Group,fill=Group)) +
  geom_boxplot(outlier.shape=NA,fill=NA) + geom_beeswarm(cex=1.5,shape=21) +
  geom_text(data=nobc_pv_cld,aes(label=.group),y=1.1,color="black") +
  ylab("Pollen viability proportion") + xlab(NULL) +
  scale_y_continuous(limits=c(0,1.1),breaks=c(0:4)/4) +
  scale_color_manual(limits=c("Brevipes","Fremontii","BxF"),values=c("#CABA27","#CE80C6","black"),guide="none") +
  scale_fill_manual(limits=c("Brevipes","Fremontii","BxF"),values=c("#CABA27","#CE80C6","black"),guide="none") +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust = 1)) + theme(text=element_text(size=12))
nobc_pv_plot

### Pollen counts model, no backcrosses
nobc_pc_model<-glmer(round(Pollen_per_square) ~ Group + (1|Population_Line),data=nobc_individual,family=poisson(link="log"))
nobc_pc_emmeans<-emmeans(nobc_pc_model,"Group")
nobc_pc_cld<-cld(nobc_pc_emmeans,Letters="abcd")

nobc_pc_plot<-ggplot(nobc_individual,aes(x=Group,y=round(Pollen_per_square),color=Group,fill=Group)) +
  geom_boxplot(outlier.shape=NA,fill=NA) + geom_beeswarm(cex=1.5,shape=21) +
  geom_text(data=nobc_pc_cld,aes(label=.group),y=150,color="black") + ylim(0,150) +
  ylab("Total pollen abundance (count per mm2)") + xlab(NULL) +
  scale_color_manual(limits=c("Brevipes","Fremontii","BxF"),
                     values=c("#CABA27","#CE80C6","black"),guide="none") +
  scale_fill_manual(limits=c("Brevipes","Fremontii","BxF"),
                    values=c("#CABA27","#CE80C6","black"),guide="none") +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust = 1)) + theme(text=element_text(size=12))
nobc_pc_plot

###combine both and save plot
both<-align_plots(nobc_pc_plot,nobc_pv_plot,align="hv",axis="tblr")
png("./Outputs/pollen_viability_nobc.png",width=6,height=6,units="in",res=600)
plot_grid(plotlist=both,nrow=1)
dev.off()
