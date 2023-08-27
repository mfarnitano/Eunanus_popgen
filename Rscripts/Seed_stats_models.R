
### Libraries
library(tidyverse)
library(cowplot)
library(beeswarm)
library(lme4)
library(brglm)
library(brglm2)
library(emmeans)
library(multcomp)
library(multcompView)

###Load seed data and process
seedcounts<-read.csv("./Data/Brevipes_seed_counts_combined.csv") %>%
  mutate(Good_seeds=as.numeric(Good_seeds),Shriveled_seeds=as.numeric(Shriveled_seeds),
         Total_seeds=Good_seeds+Shriveled_seeds,Proportion_viable=Good_seeds/Total_seeds) %>%
  mutate(Mat_Species = str_to_upper(str_sub(Mat_Species,1,4)),
         Pat_Species = str_to_upper(str_sub(Pat_Species,1,4)),
         Cross_type=str_replace_all(Cross_type,"[VC]","S")) %>%
  mutate(Cross_type=recode(Cross_type,"JS"="JxJ","BS"="BxB","FS"="FxF","SS"="SxS"),
         Mat_Species = factor(Mat_Species,levels=c("JOHN","BREV","FREM","SESP")),
         Pat_Species = factor(Pat_Species,levels=c("JOHN","BREV","FREM","SESP")),
         fruit_ID=as.factor(c(1:length(Batch))))

correctlevels<-c("JxJ","JxB","JxF","JxS","BxJ","BxB","BxF","BxS","FxJ","FxB","FxF","FxS","SxJ","SxB","SxF","SxS")

data_all<-filter(seedcounts,Total_seeds>0,Mat_Species %in% c("JOHN","BREV","FREM","SESP"),Pat_Species %in% c("JOHN","BREV","FREM","SESP")) %>%
  mutate(Cross_type=factor(Cross_type,levels=correctlevels))


Colors<-read.csv(".Data/Hex_key_brevgroup.csv") %>% 
  filter(Level=="Species",Abbreviation %in% data_all$Mat_Species) %>%
  arrange(match(Abbreviation,levels(data_all$Mat_Species)))

Labels<-read.csv(".Data/cross-labels.csv") %>%
  mutate(Cross_type=factor(Cross_type,levels=correctlevels))

### Seed viability: firth model

bigfirth<-brglm(cbind(Good_seeds,Shriveled_seeds) ~ Cross_type, data=data_all,family=binomial("logit"),
                control.brglm=list(br.maxit=5000,br.epsilon = 1e-08, br.trace=FALSE,br.consts = NULL))
bigfirth_emmeans<-emmeans(bigfirth,"Cross_type",type="response")
bigfirth_cld<-cld(bigfirth_emmeans,Letters="abcdefghijklm") 
bigfirth_bars<-ggplot(bigfirth_cld,aes(x=Cross_type,y=prob,ymin=asymp.LCL,ymax=asymp.UCL,fill=Cross_type,label=.group)) + 
  geom_col() + geom_errorbar() + geom_text(y=1.1) +
  geom_vline(xintercept=4.5,linetype="dashed") + geom_vline(xintercept=8.5,linetype="dashed") + geom_vline(xintercept=12.5,linetype="dashed") + 
  scale_fill_manual(values=c(rep(Colors$Hexcolor[1:4],3),Colors$Hexcolor[1:2],Colors$Hexcolor[4]),guide="none") + 
  scale_y_continuous(limits=c(0,1.1),breaks=c(0:4)/4) + 
  theme_bw()
bigfirth_boxes<-ggplot(data_all,aes(x=Cross_type,y=Proportion_viable,fill=Cross_type)) + 
  geom_boxplot(outlier.shape=NA) + geom_beeswarm(cex=0.25) + 
  geom_text(data=bigfirth_cld,aes(x=Cross_type,label=.group,y=NULL),y=1.1) + 
  geom_vline(xintercept=4.5,linetype="dashed") + geom_vline(xintercept=8.5,linetype="dashed") + geom_vline(xintercept=12.5,linetype="dashed") + 
  scale_fill_manual(values=c(rep(Colors$Hexcolor[1:4],3),Colors$Hexcolor[1:2],Colors$Hexcolor[4]),guide="none") + 
  scale_y_continuous(limits=c(0,1.1),breaks=c(0:4)/4) + 
  theme_bw()
png("./Outputs/firthmodel_seed_viability_complete_bars.png",width=6,height=3,units="in",res=600)
bigfirth_bars
dev.off()
png("./Outputs/firthmodel_seed_viability_complete_boxes.png",width=6,height=3,units="in",res=600)
bigfirth_boxes
dev.off()

###Seed viability: GLM (no firth) model
nofirth<-glm(cbind(Good_seeds,Shriveled_seeds) ~ Cross_type, data=data_all,family=binomial("logit"))
nofirth_emmeans<-emmeans(nofirth,"Cross_type",type="response")
nofirth_cld<-cld(nofirth_emmeans,Letters="abcdefghijklm") 
nofirth_bars<-ggplot(nofirth_cld,aes(x=Cross_type,y=prob,ymin=asymp.LCL,ymax=asymp.UCL,fill=Cross_type,label=.group)) + 
  geom_col() + geom_errorbar() + geom_text(y=1.0) +
  geom_vline(xintercept=4.5,linetype="dashed") + geom_vline(xintercept=8.5,linetype="dashed") + geom_vline(xintercept=12.5,linetype="dashed") + 
  scale_fill_manual(values=c(Colors$Hexcolor[1],rep("grey",4),Colors$Hexcolor[2],rep("grey",4),Colors$Hexcolor[3],rep("grey",3),Colors$Hexcolor[4])) + 
  theme_bw()
nofirth_bars

###Fruit set: firth model
data_full<-filter(seedcounts,Mat_Species %in% c("JOHN","BREV","FREM","SESP"),Pat_Species %in% c("JOHN","BREV","FREM","SESP")) %>%
  mutate(Cross_type=factor(Cross_type,levels=correctlevels),
         fruit_ID=as.factor(c(1:n())),
         fruit_produced=(Total_seeds>0))

firth_seed_production<-brglm(fruit_produced ~ Cross_type,data=data_full,family=binomial("logit"))
sp_firth_emmeans<-emmeans(firth_seed_production,"Cross_type",type="response")
sp_firth_cld<-cld(sp_firth_emmeans,Letters="abcde")

sp_firth_bars<-ggplot(sp_firth_cld,aes(x=Cross_type,y=prob,ymin=asymp.LCL,ymax=asymp.UCL,fill=Cross_type,label=.group)) + 
  geom_col() + geom_errorbar(width=0.5) + geom_text(y=1.1) +
  geom_vline(xintercept=4.5,linetype="dashed") + geom_vline(xintercept=8.5,linetype="dashed") + geom_vline(xintercept=12.5,linetype="dashed") +
  scale_fill_manual(values=rep(Colors$Hexcolor[1:4],4),guide="none") + 
  scale_y_continuous(limits=c(0,1.1),breaks=c(0:4)/4) + 
  theme_bw() 

png("./Outputs/firthmodel_seed_production_complete_bars.png",width=6,height=3,units="in",res=600)
sp_firth_bars
dev.off()


###Hybrid seed production models
seedcounts_withHybrids<-read.csv("./Data/Brevipes_seed_counts_combined.csv") %>%
  mutate(Good_seeds=as.numeric(Good_seeds),Shriveled_seeds=as.numeric(Shriveled_seeds),
         Total_seeds=Good_seeds+Shriveled_seeds,Proportion_viable=Good_seeds/Total_seeds) %>%
  filter(Mat_Species %in% c("Brevipes","Fremontii","BxF"), Pat_Species %in% c("Brevipes","Fremontii","BxF")) %>%
  mutate(Mat_Species = factor(recode(Mat_Species,"Brevipes"="BREV","Fremontii"="FREM"),levels=c("BREV","BxF","FREM")),
         Pat_Species = factor(recode(Pat_Species,"Brevipes"="BREV","Fremontii"="FREM"),levels=c("BREV","BxF","FREM")),
         Cross_type=factor(recode(Cross_type,"BS"="BxB","FS"="FxF","HS"="HxH"),
                           levels=c("BxB","BxH","BxF","HxB","HxH","HxF","FxB","FxH","FxF")),
         fruit_ID=as.factor(c(1:length(Batch))),
         fruit_produced=as.integer(Total_seeds>0))

withhybrids<-glmer(fruit_produced ~ Cross_type + (1|Mat_Family),data=seedcounts_withHybrids,family=binomial("logit"))
withhybrids_emmeans<-emmeans(withhybrids,"Cross_type",type="response")
withhybrids_cld<-cld(withhybrids_emmeans,Letters="abcde") %>% 
  mutate(Cross_type=factor(Cross_type,levels=c("BxB","BxH","BxF","HxB","HxH","HxF","FxB","FxH","FxF")))

withhybrids_bars<-ggplot(withhybrids_cld,aes(x=Cross_type,y=prob,ymin=asymp.LCL,ymax=asymp.UCL,fill=Cross_type,label=.group)) + 
  geom_col() + geom_errorbar(width=0.5) + geom_text(y=1.1) +
  geom_vline(xintercept=3.5,linetype="dashed") + geom_vline(xintercept=6.5,linetype="dashed") +
  scale_fill_manual(values=rep(c(Colors$Hexcolor[2],"lightgrey",Colors$Hexcolor[3]),3),guide="none") + 
  scale_y_continuous(limits=c(0,1.1),breaks=c(0:4)/4) + 
  theme_bw() 

png("./Outputs/seed_production_hybrids.png",width=3,height=3,units="in",res=600)
withhybrids_bars
dev.off()

### Hybrid seed counts models

counts_model<-glm.nb(Total_seeds ~ Cross_type,data=seedcounts_withHybrids) ###negative binomial model
counts_emmeans<-emmeans(counts_model,"Cross_type",type="response")
counts_cld<-cld(counts_emmeans,Letters="abcdefghijkl") %>% 
  mutate(Cross_type=factor(Cross_type,levels=c("BxB","BxH","BxF","HxB","HxH","HxF","FxB","FxH","FxF")))

countmodel_bars<-ggplot(counts_cld,aes(x=Cross_type,y=response,ymin=asymp.LCL,ymax=asymp.UCL,fill=Cross_type,label=.group)) + 
  geom_col() + geom_errorbar(width=0.5) + geom_text(y=580) +
  geom_vline(xintercept=3.5,linetype="dashed") + geom_vline(xintercept=6.5,linetype="dashed") +
  scale_fill_manual(values=rep(c(Colors$Hexcolor[2],"lightgrey",Colors$Hexcolor[3]),3),guide="none") + 
  ylim(0,600) + 
  theme_bw() 

countmodel_bars

### Hybrid seed counts, separate models for each paternal species
completelevels<-c("BxB","BxH","BxF","HxB","HxH","HxF","FxB","FxH","FxF")
matplotlist<-list()
boxplotlist<-list()
for (pat in c("BREV","BxF","FREM")) {
  subdata<-filter(seedcounts_withHybrids,Pat_Species == pat)
  submodel<-glm.nb(Total_seeds ~ Cross_type,data=subdata)
  sub_emmeans<-emmeans(submodel,"Cross_type",type="response")
  sub_cld<-cld(sub_emmeans,Letters="abcdefghijkl") %>% 
    mutate(Cross_type=factor(Cross_type,levels=completelevels[Cross_type %in% completelevels]))
  subplot<-ggplot(sub_cld,aes(x=Cross_type,y=response,ymin=asymp.LCL,ymax=asymp.UCL,fill=Cross_type,label=.group)) + 
    geom_col() + geom_errorbar(width=0.5) + geom_text(y=280) +
    scale_fill_manual(values=c(Colors$Hexcolor[2],"lightgrey",Colors$Hexcolor[3]),guide="none") + 
    ylim(0,300) + ylab(NULL) + 
    theme_bw() 
  matplotlist[[pat]]<-subplot
  subplot_data<-ggplot(subdata,aes(x=Cross_type,y=Total_seeds,fill=Cross_type)) + geom_violin(scale = "width") + 
    geom_text(data=sub_cld,aes(label=.group),y=450) + 
    scale_fill_manual(values=c(Colors$Hexcolor[2],"lightgrey",Colors$Hexcolor[3]),guide="none") + 
    ylim(0,460) + ylab(NULL) + 
    theme_bw()
  boxplotlist[[pat]]<-subplot_data
}
paternal_counts<-align_plots(plotlist=boxplotlist,axis="tb",align="h")

png("./Outputs/seed_production_hybrids_bypat.png",width=5,height=5,units="in",res=600)
plot_grid(plotlist=paternal_counts,nrow=1)
dev.off()

###seed counts for female fertility, compare only four cross types
uselevels<-c("BxB","HxB","HxF","FxF")
fertilitydata<-filter(seedcounts_withHybrids,Cross_type %in% uselevels,Total_seeds>0)
fertilitymodel<-lm(Total_seeds ~ Cross_type,data=fertilitydata)
fertility_emmeans<-emmeans(fertilitymodel,"Cross_type")
fertility_cld<-cld(fertility_emmeans,Letters="abcde") %>% 
  mutate(Cross_type=factor(Cross_type,levels=uselevels))
fertilityplot_data<-ggplot(fertilitydata,aes(x=Cross_type,y=Total_seeds,color=Cross_type)) + geom_boxplot(outlier.shape=NA) + 
  geom_beeswarm() + 
  geom_text(data=fertility_cld,aes(label=.group,color=NULL),y=450,color="black") + 
  scale_color_manual(values=c("#CABA27","black","black","#CE80C6"),guide="none") + 
  ylim(0,460) + ylab(NULL) + 
  theme_bw()
fertilityplot_data<-ggplot(fertilitydata,aes(x=Cross_type,y=Total_seeds,color=Cross_type)) + geom_violin() + 
  geom_text(data=fertility_cld,aes(label=.group,color=NULL),y=450,color="black") + 
  scale_color_manual(values=c("#CABA27","black","black","#CE80C6"),guide="none") + 
  ylim(0,460) + ylab(NULL) + 
  theme_bw()
fertilityplot_data

png("./Outputs/seed_production_hybrids_fertility.png",width=4,height=6,units="in",res=600)
fertilityplot_data
dev.off()

###seed counts, relevant comparisons for postmating prezygotic barriers
uselevels<-c("BxB","BxH","HxB","HxH","HxF","FxH","FxF","BxF","FxB")
PMPZ_hybrid_data<-filter(seedcounts_withHybrids,Cross_type %in% uselevels) %>% mutate(Cross_type=factor(Cross_type,levels=uselevels))
PMPZ_hybrid_model<-glm.nb(Total_seeds ~ Cross_type,data=PMPZ_hybrid_data)
PMPZ_hybrid_emmeans<-emmeans(PMPZ_hybrid_model,"Cross_type",type="response")
PMPZ_hybrid_cld<-cld(PMPZ_hybrid_emmeans,Letters="abcdefghij") %>% 
  mutate(Cross_type=factor(Cross_type,levels=uselevels))
PMPZ_hybridplot_data<-ggplot(PMPZ_hybrid_data,aes(x=Cross_type,y=Total_seeds,color=Cross_type)) + geom_boxplot(outlier.shape=NA) + 
  geom_point() + 
  geom_text(data=PMPZ_hybrid_cld,aes(label=.group,color=NULL),y=450,color="black") + 
  geom_vline(xintercept=7.5,linetype="dashed") +
  scale_color_manual(values=c("#CABA27","grey","black","black","grey","black","#CE80C6","grey","black"),guide="none") + 
  ylim(0,460) + ylab(NULL) + 
  theme_bw()
PMPZ_hybridplot_data

png("./Outputs/seed_production_PMPZ_hybrids.png",width=4,height=6,units="in",res=600)
PMPZ_hybridplot_data
dev.off()

###seed counts for PMPZ, make relevant comparisons, exclude zeros
uselevels<-c("BxB","BxH","HxB","HxH","HxF","FxH","FxF","BxF","FxB")
PMPZ_seed_data<-filter(seedcounts_withHybrids,Cross_type %in% uselevels,Total_seeds>0) %>% 
  mutate(Cross_type=factor(Cross_type,levels=uselevels))
PMPZ_seed_model<-lm(Total_seeds ~ Cross_type,data=PMPZ_seed_data)
PMPZ_seed_emmeans<-emmeans(PMPZ_seed_model,"Cross_type")
PMPZ_seed_cld<-cld(PMPZ_seed_emmeans,Letters="abcdefghij") %>% 
  mutate(Cross_type=factor(Cross_type,levels=uselevels))
PMPZ_seedplot_data<-ggplot(PMPZ_seed_data,aes(x=Cross_type,y=Total_seeds,color=Cross_type)) + geom_boxplot(outlier.shape=NA) + 
  geom_point() + 
  geom_text(data=PMPZ_seed_cld,aes(label=.group,color=NULL),y=450,color="black") + 
  geom_vline(xintercept=7.5,linetype="dashed") +
  scale_color_manual(values=c("#CABA27","grey","black","black","grey","black","#CE80C6","grey","black"),guide="none") + 
  ylim(0,460) + ylab(NULL) + 
  theme_bw()
PMPZ_seedplot_data

png("./Outputs/seed_production_PMPZ_hybrids.png",width=4,height=6,units="in",res=600)
PMPZ_hybridplot_data
dev.off()

###PMPZ cross success, make relevant comparisons
uselevels<-c("BxB","BxH","HxB","HxH","HxF","FxH","FxF","BxF","FxB")
PMPZ_withhybrids<-glmer(fruit_produced ~ Cross_type + (1|Mat_Family),data=seedcounts_withHybrids,family=binomial("logit"))
PMPZ_withhybrids_emmeans<-emmeans(PMPZ_withhybrids,"Cross_type",type="response")
PMPZ_withhybrids_cld<-cld(PMPZ_withhybrids_emmeans,Letters="abcde") %>% 
  mutate(Cross_type=factor(Cross_type,levels=uselevels))

PMPZ_withhybrids_bars<-ggplot(PMPZ_withhybrids_cld,aes(x=Cross_type,y=prob,ymin=asymp.LCL,ymax=asymp.UCL,color=Cross_type,label=.group)) + 
  geom_col(fill="white") + geom_errorbar(width=0.5) + geom_text(color="black",y=1.1) +
  geom_vline(xintercept=7.5,linetype="dashed") +
  scale_color_manual(values=c("#CABA27","grey","black","black","grey","black","#CE80C6","grey","black"),guide="none") + 
  scale_y_continuous(limits=c(0,1.1),breaks=c(0:4)/4) + 
  theme_bw() 

png("./Outputs/PMPZ_cross_success_hybrids.png",width=5,height=3,units="in",res=600)
PMPZ_withhybrids_bars
dev.off()

###Seed counts model -- main crosses (given fruit set>0) 
seed_model<-glm.nb(Total_seeds ~ Cross_type,data=data_all)
seed_emmeans<-emmeans(seed_model,"Cross_type",type="response")
seed_cld<-cld(seed_emmeans,Letters="abcdefghij") %>% 
  mutate(Cross_type=factor(Cross_type,levels=correctlevels))
seedplot_data<-ggplot(data_all,aes(x=Cross_type,y=Total_seeds,color=Cross_type)) + geom_boxplot(outlier.shape=NA) + 
  geom_point() + 
  geom_text(data=seed_cld,aes(label=.group,color=NULL),y=450,color="black") + 
  geom_vline(xintercept=4.5,linetype="dashed") + geom_vline(xintercept=8.5,linetype="dashed") + geom_vline(xintercept=12.5,linetype="dashed") +
  scale_color_manual(values=c(rep(Colors$Hexcolor[1:4],3),Colors$Hexcolor[1:2],Colors$Hexcolor[4]),guide="none") + ylim(0,460) + ylab(NULL) + 
  theme_bw()
seedplot_data

png("./Outputs/seed_production_model.png",width=6,height=4,units="in",res=600)
seedplot_data
dev.off()

###Load seed size data
seed_sizes<-read.csv("./Data/Seed_sizes_measured.csv") %>%
  filter(Maternal!="") %>%
  mutate(cross=paste0(Maternal,"x",Paternal), 
         fruit=paste0(cross,"_",Fruit_ID),
         lw_ratio=Length/Width,
         Mat_species = factor(Mat_Species,levels=c("JOHN","BREV","FREM","SESP")),
         Pat_species = factor(Pat_species,levels=c("JOHN","BREV","FREM","SESP"))) %>%
  mutate(Cross_type=factor(Cross_type,levels=correctlevels)) %>%
  mutate(ellipse_area=pi*Length*Width/4)

seed_sizes_fruitmeans<-seed_sizes %>% group_by(Mat_species,Pat_species,Cross_type,fruit) %>%
  summarize(mean_Length=mean(Length),mean_lwratio=mean(lw_ratio),mean_ellipse_area=mean(ellipse_area))

###seed length: full model and plot
length_model<-lmer(Length ~ Cross_type + (1|fruit),data=seed_sizes)
length_emmeans<-emmeans(length_model,"Cross_type")
length_cld<-cld(length_emmeans,Letters="abcdefghi")
length_bars<-ggplot(length_cld,aes(x=Cross_type,y=emmean,ymin=lower.CL,ymax=upper.CL,fill=Cross_type,label=.group)) + 
  geom_col() + geom_errorbar(width=0.5) + geom_text(y=1.1) +
  geom_vline(xintercept=4.5,linetype="dashed") + geom_vline(xintercept=8.5,linetype="dashed") + geom_vline(xintercept=12.5,linetype="dashed") +
  scale_fill_manual(values=rep(Colors$Hexcolor[1:4],4),guide="none") + 
  scale_y_continuous(limits=c(-0.1,1.1),breaks=c(0:4)/4) + 
  theme_bw() 

png("./Outputs/length_model_bars.png",width=6,height=3,units="in",res=600)
length_bars
dev.off()

###seed length models: pairwise comparisons
modellist_lengths<-list()
emslist_lengths<-list()
plotlist_lengths<-list()
for (pair in list(c("JOHN","BREV"),c("JOHN","FREM"),c("BREV","FREM"),c("JOHN","SESP"),c("BREV","SESP"),c("FREM","SESP"))) {
  print(paste(pair,collapse="x"))
  #set up data
  pairname=str_c(pair,collapse="x")
  subdata<-seed_sizes %>% filter(Mat_species %in% pair,Pat_species %in% pair) 
  
  #linear model
  sub_model<-lmer(Length~Cross_type + (1|fruit),data=subdata)
  #sub_aov<-aov(sub_model)
  sub_ems<-emmeans(sub_model,"Cross_type")
  sub_cld<-cld(sub_ems,Letters="abcd")
  
  #set up means and CIs
  sub_cld<-sub_cld %>% mutate(Cross_type=factor(Cross_type,levels=correctlevels[correctlevels %in% sub_cld$Cross_type]))
  
  #prep colors
  sub_colors<-Colors %>% filter(Abbreviation %in% pair)
  
  #prep labels
  labeller<-Labels %>% filter(Cross_type %in% levels(sub_cld$Cross_type))
  
  #make plots
  #linear model lsmeans and CIs
  sub_plot<-ggplot(subdata,aes(x=Cross_type,y=Length,fill=Cross_type)) + 
    geom_violin() + 
    #geom_boxplot(outlier.shape=NA,lwd=0.1) + 
    #geom_beeswarm(size=0.1) + 
    geom_text(data=sub_cld,aes(label=.group),y=1.15,size=2) +
    scale_fill_manual(values=c(sub_colors$Hexcolor[1],"lightgrey","grey",sub_colors$Hexcolor[2]),guide="none") + 
    theme_bw() + xlab(NULL) + ylab(NULL) + 
    theme(axis.text.x = element_text(angle=90,hjust=1)) + 
    #scale_x_discrete(labels=labeller$Cross_type) +
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0),limits=c(-0.1,1.2))
  
  #save data and plots
  modellist_lengths[[pairname]]<-sub_model
  #aovlist_lengths[[pairname]]<-sub_aov
  emslist_lengths[[pairname]]<-sub_ems
  plotlist_lengths[[pairname]]<-sub_plot
}
#plotlist_lengths$FREMxSESP<-plotlist_lengths$FREMxSESP + 
#  scale_x_discrete(limits=c("FxF","FxS","SxF","SxS"),labels=Labels$Cross_type[c(11,12,15,16)]) + 
#  scale_fill_manual(values=c(Colors$Hexcolor[3],"lightgrey",Colors$Hexcolor[4]),guide="none")
aligned_lengths<-align_plots(plotlist=plotlist_lengths,axis="tb",align="h")


png("./Outputs/Seed_lengths_pairwise.png",width=6,height=2,units="in",res=600)
plot_grid(plotlist=aligned_lengths,nrow=1)
dev.off()

### seed length/width ratios: full model and plot
ratio_model<-lmer(lw_ratio ~ Cross_type + (1|fruit),data=seed_sizes)
ratio_emmeans<-emmeans(ratio_model,"Cross_type")
ratio_cld<-cld(ratio_emmeans,Letters="abcdefghi")
ratio_bars<-ggplot(ratio_cld,aes(x=Cross_type,y=emmean,ymin=lower.CL,ymax=upper.CL,fill=Cross_type,label=.group)) + 
  geom_col() + geom_errorbar(width=0.5) + geom_text(y=5.5) +
  geom_vline(xintercept=4.5,linetype="dashed") + geom_vline(xintercept=8.5,linetype="dashed") + geom_vline(xintercept=12.5,linetype="dashed") +
  scale_fill_manual(values=rep(Colors$Hexcolor[1:4],4),guide="none") + 
  scale_y_continuous(limits=c(0,6)) + 
  theme_bw() 

png("./Outputs/ratio_model_bars.png",width=6,height=3,units="in",res=600)
ratio_bars
dev.off()


###Seed length/width ratios: pairwise comparisons
modellist_lw<-list()
emslist_lw<-list()
plotlist_lw<-list()
for (pair in list(c("JOHN","BREV"),c("JOHN","FREM"),c("BREV","FREM"),c("JOHN","SESP"),c("BREV","SESP"),c("FREM","SESP"))) {
  print(paste(pair,collapse="x"))
  #set up data
  pairname=str_c(pair,collapse="x")
  subdata<-seed_sizes %>% filter(Mat_species %in% pair,Pat_species %in% pair) 
  
  #linear model
  sub_model<-lmer(lw_ratio~Cross_type + (1|fruit),data=subdata)
  sub_ems<-emmeans(sub_model,"Cross_type")
  sub_cld<-cld(sub_ems,Letters="abcd")
  
  #set up means and CIs
  sub_cld<-sub_cld %>% mutate(Cross_type=factor(Cross_type,levels=correctlevels[correctlevels %in% sub_cld$Cross_type]))
  
  #prep colors
  sub_colors<-Colors %>% filter(Abbreviation %in% pair)
  
  #prep labels
  labeller<-Labels %>% filter(Cross_type %in% levels(sub_cld$Cross_type))
  
  #make plots
  #linear model lsmeans and CIs
  sub_plot<-ggplot(subdata,aes(x=Cross_type,y=lw_ratio,fill=Cross_type)) + 
    geom_violin() + 
    #geom_boxplot(outlier.shape=NA,lwd=0.1) + 
    #geom_beeswarm(size=0.1) + 
    geom_text(data=sub_cld,aes(label=.group),y=7.5,size=2) +
    scale_fill_manual(values=c(sub_colors$Hexcolor[1],"lightgrey","grey",sub_colors$Hexcolor[2]),guide="none") + 
    theme_bw() + xlab(NULL) + ylab(NULL) + 
    theme(axis.text.x = element_text(angle=90,hjust=1)) + 
    #scale_x_discrete(labels=labeller$Cross_type) +
    scale_y_continuous(,limits=c(0,8))
  
  #save data and plots
  modellist_lw[[pairname]]<-sub_model
  emslist_lw[[pairname]]<-sub_ems
  plotlist_lw[[pairname]]<-sub_plot
}
#plotlist_lengths$FREMxSESP<-plotlist_lengths$FREMxSESP + 
#  scale_x_discrete(limits=c("FxF","FxS","SxF","SxS"),labels=Labels$Cross_type[c(11,12,15,16)]) + 
#  scale_fill_manual(values=c(Colors$Hexcolor[3],"lightgrey",Colors$Hexcolor[4]),guide="none")
aligned_lw<-align_plots(plotlist=plotlist_lw,axis="tb",align="h")
png("./Outputs/Seed_lw_ratios_pairwise.png",width=6,height=2,units="in",res=600)
plot_grid(plotlist=aligned_lw,nrow=1)
dev.off()

###calculating total reproductive isolation

#fruit set
fruit_set_summary<-data_full %>% group_by(Mat_Species,Pat_Species,Cross_type) %>% 
  summarize(fruit_set=sum(fruit_produced)/n()) %>% ungroup() %>%
  left_join(sp_firth_cld,by=c("Cross_type"))
conspecifics<-filter(fruit_set_summary,Mat_Species==Pat_Species) %>% dplyr::select(Pat_Species,fruit_set,prob)
get_C<-function(species) {
  getrow<-filter(conspecifics,Pat_Species == species)
  return(getrow$fruit_set)[1]
}
get_C2<-function(species) {
  getrow<-filter(conspecifics,Pat_Species == species)
  return(getrow$prob)[1]
}
fruitset_RI<-fruit_set_summary %>% mutate(H_fs=fruit_set,C_fs=sapply(Mat_Species,get_C),RI_fs=(1-(2*(H_fs/(H_fs+C_fs))))) %>%
  mutate(H_fs_model=prob,C_fs_model=sapply(Mat_Species,get_C2),RI_fs_model=(1-(2*(H_fs_model/(H_fs_model+C_fs_model))))) %>%
  filter(Mat_Species != Pat_Species) %>% dplyr::select(Mat_Species,Pat_Species,RI_fs,RI_fs_model)


##seed viability
seed_viability_summary<-data_all %>% group_by(Mat_Species,Pat_Species,Cross_type) %>%
  summarize(viability_proportion=sum(Good_seeds)/(sum(Shriveled_seeds)+sum(Good_seeds))) %>% ungroup()
conspecifics_sv<-filter(seed_viability_summary,Mat_Species==Pat_Species) %>% mutate(value=viability_proportion) %>% 
  dplyr::select(Pat_Species,viability_proportion)
get_C_sv<-function(species) {
  getrow<-filter(conspecifics_sv,Pat_Species == species)
  return(getrow$viability_proportion)[1]
}
seedviability_RI<-seed_viability_summary %>%
  mutate(H_sv=viability_proportion,C_sv=sapply(Mat_Species,get_C_sv),RI_sv=(1-(2*(H_sv/(H_sv+C_sv))))) %>%
  filter(Mat_Species != Pat_Species) %>% dplyr::select(Mat_Species,Pat_Species,RI_sv)



##pollen viability
pollen_data<-read.csv("./Data/Pollen_hemocytometer_counts.csv") %>% filter(Group %in% c("Brevipes","BxF","Fremontii")) %>% 
  group_by(Group) %>% summarize(pollen_viability=sum(n_viable)/(sum(n_inviable)+sum(n_viable))) %>% ungroup()
1-(2*(0.227/(0.227+0.839)))

