
### Libraries
library(tidyverse)

### function to process twisst output file

twisst_summary<-function(inputfile,directory) {
  name=str_sub(inputfile,1,4)
  weights=read.table(paste0(directory,inputfile),skip=3,header=T)
  topos=data.frame(toporaw=read_lines(paste0(directory,inputfile),n_max=3)) %>% 
    separate(toporaw,into=c("topoID","topology"),sep=' ') %>% 
    mutate(topoID=str_remove(topoID,"#"))
  rowtotal=sum(weights[1,])
  summary<-weights %>% summarise(across(topo1:topo3,sum)) %>% 
    pivot_longer(cols=topo1:topo3,names_to="topoID",values_to="rawcount") %>%
    mutate(treecount=rawcount/rowtotal,proportion=rawcount/sum(rawcount)) %>% 
    left_join(topos,by="topoID")
  bintest<-binom.test(round(arrange(summary,desc(treecount))$treecount[2:3]))
  wide<-summary %>% pivot_wider(names_from=topoID,values_from=rawcount:topology,names_sep="_") %>% 
    mutate(name=name,estimate=bintest$estimate,pvalue=bintest$p.value,lower=bintest$conf.int[1],upper=bintest$conf.int[2])
  return(wide)
}

### Run function for all inputs

directory="~/OneDrive - University of Georgia/Sweigart_Lab/Ch1_Brevipes_group/Results/Bioinformatic_outputs/withCONS/twisst/"
files=list.files(directory,pattern = "weights")

#View one example
twisst_summary("FVCA_by_genes_ML_weights.txt",directory)

#All quartets
alltwisst<-do.call(rbind, lapply(files, twisst_summary,directory=directory)) %>% 
  mutate(pcorr=pvalue*20,significance=ifelse(pcorr<0.001,"***",ifelse(pcorr<0.01,"**",ifelse(pcorr<0.05,"*","n.s."))))


### Window support across the genome

inputfile="FVCA_by_genes_ML_weights.txt"

allgenes<-read.table(paste0(directory,inputfile),skip=3,header=T) %>% mutate(relsupport=topo1-topo3,geneID=1:length(relsupport))
genomic<-ggplot(allgenes,aes(x=geneID,y=relsupport)) + geom_point()

histosupport<-ggplot(allgenes,aes(x=relsupport)) + geom_histogram()

inputfile2="BVNA_by_genes_ML_weights.txt"
allgenes2<-read.table(paste0(directory,inputfile2),skip=3,header=T) %>% mutate(relsupport=topo1-topo3,geneID=1:length(relsupport))
histosupport2<-ggplot(allgenes2,aes(x=relsupport)) + geom_histogram()
