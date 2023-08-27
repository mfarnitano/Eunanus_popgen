
###Libraries
library(tidyverse)


###Load 1Mb-scale data (4-fold degenerate sites) individual comparisons
hz_1Mb<-read.table("Results/Bioinformatic_outputs/withCONS/pixy/Brevipes_FINAL_withCONS.pass.30called.4D.ind.1Mb_pi.txt",header=T)
pw_1Mb<-read.table("Results/Bioinformatic_outputs/withCONS/pixy/Brevipes_FINAL_withCONS.pass.30called.4D.ind.1Mb_dxy.txt",header=T)
colorkey<-read.csv("Data/Datasheets/Hex_key_brevgroup.csv") %>%
  filter(Level=="Species")

###Get genome-wide means
hz_1Mb_sums<-hz_1Mb %>% group_by(pop) %>% 
  summarize(total_diffs=sum(as.numeric(count_diffs),na.rm=T),total_sites=sum(as.numeric(count_comparisons),na.rm=T)) %>% 
  mutate(total_value=total_diffs/total_sites) %>%
  separate(pop,into=c("population"),sep="_",extra="drop",remove=F) %>%
  mutate(species=recode(str_sub(population,1,4),"CLEV"="AURA"),comparetype="heterozygosity")
pw_1Mb_sums<-pw_1Mb %>% group_by(pop1,pop2) %>% 
  summarize(total_diffs=sum(as.numeric(count_diffs),na.rm=T),total_sites=sum(as.numeric(count_comparisons),na.rm=T)) %>% 
  mutate(total_value=total_diffs/total_sites) %>%
  separate(pop1,into=c("population1"),sep="_",extra="drop",remove=F) %>%
  separate(pop2,into=c("population2"),sep="_",extra="drop",remove=F) %>% 
  mutate(species1=recode(str_sub(population1,1,4),"CLEV"="AURA"),
         species2=recode(str_sub(population2,1,4),"CLEV"="AURA"),
         species=ifelse(species1<species2,
                        paste0(species1,"x",species2),
                        paste0(species2,"x",species1)),
         comparetype="pairwise")
hz_simple<-hz_1Mb_sums %>% ungroup() %>% dplyr::select(c(species,comparetype,total_value))
pw_simple<-pw_1Mb_sums %>% ungroup() %>% dplyr::select(c(species,comparetype,total_value))
combined<-rbind(hz_simple,pw_simple)
combined_ordered<-combined %>% 
  mutate(species=factor(species,levels=unique(species))) %>%
  mutate(color_id=ifelse(grepl("AURA",species),"AURA",
                         ifelse(grepl("NANU",species),"NANU",
                                ifelse(grepl("SESP",species),"SESP",
                                       ifelse(grepl("CONS",species),"CONS",
                                              ifelse(grepl("FREM",species),"FREM",
                                                     ifelse(grepl("JOHN",species),"JOHN","BREV")))))))
combined_species_means<-combined_ordered %>% group_by(color_id,species) %>% 
  summarize(meanvalue=mean(total_value)) %>% ungroup() %>% arrange(meanvalue) %>% 
  mutate(species=factor(species,levels=levels(combined_ordered$species)))

### Subset to only Eunanus samples and order properly for plotting
combined_Eunanus<-combined_ordered %>% filter(!(grepl("AURA",species))) %>%
  group_by(color_id,species) %>% 
  summarize(mean_value=mean(total_value), SE=sqrt(var(total_value)/n())) %>%
  mutate(species=factor(species,levels=c("JOHN","BREV","FREM","SESP","CONS","NANU","JOHNxJOHN","BREVxBREV","FREMxFREM","SESPxSESP","NANUxNANU",
                                         "BREVxJOHN","FREMxJOHN","BREVxFREM","CONSxSESP","JOHNxSESP","BREVxSESP","FREMxSESP",
                                         "CONSxJOHN","BREVxCONS","CONSxFREM","JOHNxNANU","BREVxNANU","FREMxNANU","NANUxSESP","CONSxNANU")))
Eunanus_indivs<-combined_ordered %>% filter(!(grepl("AURA",species)),!(grepl("CLEV",species))) %>%
  mutate(species=factor(species,levels=c("JOHN","BREV","FREM","SESP","CONS","NANU","JOHNxJOHN","BREVxBREV","FREMxFREM","SESPxSESP","NANUxNANU",
                                         "BREVxJOHN","FREMxJOHN","BREVxFREM","CONSxSESP","JOHNxSESP","BREVxSESP","FREMxSESP",
                                         "CONSxJOHN","BREVxCONS","CONSxFREM","JOHNxNANU","BREVxNANU","FREMxNANU","NANUxSESP","CONSxNANU")))
###View table
combined_species_means %>% arrange(meanvalue)

###Create barplot and save
pixy_barplot<-ggplot(combined_Eunanus,aes(x=species,y=mean_value,fill=species)) + geom_col() + 
  geom_point(data=Eunanus_indivs,aes(y=total_value),size=0.5) + 
  #scale_fill_manual(breaks=colorkey$Abbreviation,values=colorkey$Hexcolor,guide="none") + 
  scale_fill_manual(values=colorkey$Hexcolor[c(6,7,5,8,4,3,6,7,5,8,3,6,6,7,4,6,7,5,6,7,5,6,7,5,8,4)],guide="none") + 
  theme_bw() + 
  scale_y_continuous(breaks=c(0:13)/100) + 
  scale_x_discrete(labels=c("J","B","F","S","C","N",
                            "J","B","F","S","N",
                            "BxJ","FxJ","FxB","SxC","SxJ","SxB","SxF","CxJ","CxB","CxF","NxJ","NxB","NxF","NxS","NxC")) + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 12, angle = 90,vjust=0.5),
        axis.text.y = element_text(size=12))

png(".Outputs/pixy_barplot.png",width=6,height=3,units="in",res=600)
pixy_barplot
dev.off()