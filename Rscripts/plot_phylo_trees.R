
###Libraries
library(tidyverse)
library(ape)
library(phangorn)
library(phytools)
library(ggtree)
library(treeio)

###making and bootstrapping nj trees
fasta<-read.dna(".Data/phylo/Brevipes_submission.pass.31calledSNPs.randbase.fa",format="fasta")
phydata<-phyDat(fasta,type="DNA")
dist<-dist.ml(phydata)
newnj<-nj(dist) %>% root(AURA)
newnj_midroot<-midpoint.root(newnj)
boots<-bootstrap.phyDat(phydata, FUN =function(x) nj(dist.ml(x)), bs=1000)
withboots <- root(plotBS(newnj,boots),AURA)

### draw ng tree
bootstrap_nj_plot<-ggtree(withboots) + geom_tiplab(size=2) + geom_nodelab(size=1,nudge_x = -0.02,geom='label') + ggexpand(1.1,side='h')

png("./Outputs/simple_nj_boots.png",width=6,height=6,units="in",res=600)
bootstrap_nj_plot
dev.off()

### read and draw raxml tree
raxml_rep4<-read.raxml(".Data/phylo/RAxML_bipartitionsBranchLabels.Brevipes_submission.pass.31calledSNPs.randbase.rep4") %>%
  root(AURA)
raxml_rep4@phylo <- midpoint.root(raxml_rep4@phylo)

raxmlplot<-ggtree(raxml_rep4) + geom_tiplab(size=2) + geom_nodelab(aes(label=bootstrap),size=1,nudge_x = -0.02,geom='label') + ggexpand(1.1,side='h') + geom_treescale()
png("./Outputs/simple_raxml_boots_rep6.png",width=6,height=6,units="in",res=600)
raxmlplot
dev.off()

### fancy raxml tree
fancy_raxml_data<-as_tibble(raxml_rep4) %>%
  mutate(species=str_sub(label,1,4)) %>%
  left_join(colorkey,by=c("species" = "Abbreviation"))
fancy_raxml_plot<-ggtree(raxml_rep4) %>% flip(44,59) %<+% fancy_raxml_data +
  geom_tippoint(aes(shape=species,fill=species),size=3) +
  geom_nodepoint(aes(subset=node %in% c(51,52,60,61,62)),shape=1) +
  scale_shape_manual(breaks=colorkey$Abbreviation,values=colorkey$Shape,guide="none") +
  scale_fill_manual(breaks=colorkey$Abbreviation,values=colorkey$Hexcolor,guide="none") +
  geom_treescale() +
  ggexpand(1.1,side='h')

### draw ASTRAL tree
input_dir<-"Results/Bioinformatic_outputs/withCONS/astral/"
ASTRAL_Qc<-read.astral("./ASTRAL/Qc_astral_by_genes.namecorrected.tre") %>% root(c("CLEV_cle","AURA_gra","AURA_ari","AURA_pun","AURA_aur"))
ASTRAL_Qc@phylo$edge.length[is.na(ASTRAL_Qc@phylo$edge.length)]<-0.1
ASTRAL_Qc@data$print<-paste(round(ASTRAL_Qc@data$q1,2),round(ASTRAL_Qc@data$q2,2),round(ASTRAL_Qc@data$q3,2),sep=",")
png("./Outputs/ASTRAL_Qc.png",width=8,height=8,units='in',res=600)
ggtree(ASTRAL_Qc) + geom_label(aes(x=branch,label=print),size=2,nudge_x = -0.5) + geom_tiplab(align = T,offset = 1) + hexpand(0.3,1)
dev.off()

ASTRAL_P<-read.tree("./ASTRAL/P_astral_by_genes.namecorrected.tre") %>% root(c("CLEV_cle","AURA_gra","AURA_ari","AURA_pun","AURA_aur"))
ASTRAL_P$edge.length[is.na(ASTRAL_P$edge.length)]<-0.1

png("./Outputs/ASTRAL_P.png",width=8,height=8,units='in',res=600)
ggtree(ASTRAL_P) + geom_nodelab(nudge_x=-0.2,nudge_y=0.5,size=2) + geom_tiplab(align = T,offset = 1) + hexpand(0.3,1)
dev.off()
