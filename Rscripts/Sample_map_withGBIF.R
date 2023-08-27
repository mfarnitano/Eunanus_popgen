
###Libraries
library(tidyverse)
library(raster)
library(cowplot)
library(ggalt)
library(ggrepel)

### Data citation for GBIF.org
# Data downloaded from GBIF: 
#   GBIF.org (12 September 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.j8dxeq

### Load data
Popguide<-read.csv("./Data/Hex_key_brevgroup.csv")
GBIF_data_raw<-read.csv("./Data/GBIF_09122022_download.csv") %>% 
  dplyr::select(gbifID,species,decimalLatitude,decimalLongitude,day,month,year,stateProvince) %>%
  mutate(month=as.factor(month))
GBIF_data <- GBIF_data_raw %>% 
  filter(stateProvince %in% c("California","Oregon","Washington",
                              "Idaho","Nevada","nevada", "Arizona", 
                              "Baja california","Baja California",
                              "Baja California (norte)","Baja california sur"),
         decimalLatitude > 25,decimalLongitude < -110,decimalLatitude > 33.5 | decimalLongitude > -118) %>%
  mutate(species=factor(species,levels=c("Diplacus nanus","Diplacus constrictus","Diplacus brevipes","Diplacus fremontii","Diplacus johnstonii"))) %>% 
  arrange(species)

Colors<-Popguide %>% filter(Level=="Species",Abbreviation %in% c("BREV","FREM","JOHN","CONS","NANU","SESP")) %>% 
  dplyr::select(Abbreviation, Species, Hexcolor, M_species, Shape) %>% 
  mutate(Species=factor(Species,levels=c("Diplacus nanus","Diplacus constrictus","Diplacus brevipes","Diplacus fremontii","Diplacus johnstonii","Diplacus sp. Sespe creek"))) %>% 
  arrange(Species)
Samples<-Popguide %>% filter(Level=="Population",!is.na(Longitude),!(Abbreviation %in% c("JOHN25","BIGE25"))) %>%
  dplyr::select(Abbreviation, Species, Hexcolor, M_species, Shape,Longitude,Latitude,Shortname) %>% 
  mutate(Species=factor(Species,levels=c("Diplacus nanus","Diplacus constrictus","Diplacus brevipes","Diplacus fremontii","Diplacus johnstonii","Diplacus sp. Sespe creek")),
         textcolor=ifelse(Shortname %in% c("S25","C01"),"light","dark"))

GBIF_small<-GBIF_data %>% filter(stateProvince %in% c("California","Baja California","Baja california","Baja California (norte)"))

### Get state maps raster
us <- raster::getData("GADM",country="USA",level=1)
us.states<-us[us$NAME_1 %in% c("California"),]
us.states.plus<-us[us$NAME_1 %in% c("California","Oregon","Idaho","Nevada","Arizona"),]
mex <- raster::getData("GADM",country="MEX",level=1)
mex.states<-mex[mex$NAME_1 %in% c("Baja California"),]

###create full map
base.plus<-ggplot(us.states.plus,aes(x=long,y=lat,group=group)) + geom_path() + 
  geom_path(data=mex.states) + theme_bw() + coord_fixed(ratio=3/2,xlim=c(-125.1,-113.9),ylim=c(29.9,43.1),expand = F) + 
  geom_rect(aes(xmin=c(-121),xmax=c(-116),ymin=c(32.5),ymax=c(35)),
            color="black",alpha=0) + 
  geom_point(data=GBIF_small,
             aes(x=decimalLongitude,y=decimalLatitude,group=NULL,fill=species,shape=species),
             size=1,alpha=0.5,color="transparent",stroke=0) + 
  scale_fill_manual(values=Colors$Hexcolor,limits=Colors$Species,guide="none") +
  scale_shape_manual(values=Colors$Shape,limits=Colors$Species,guide="none") + 
  #scale_x_continuous(limits=c(-125,-114),breaks=c(-62:-57)*2,expand=c(0,0)) + 
  #scale_y_continuous(limits=c(30,43),expand=c(0,0)) + 
  ylab(NULL) + xlab(NULL)

###create inset
inset.plus<-base.plus + 
  geom_label_repel(data=subset(Samples,Abbreviation != "NANU"),
                   aes(x=Longitude,y=Latitude,label=Shortname,group=NULL,fill=Species,color=textcolor),
                   size=3) + 
  scale_color_manual(limits=c("light","dark"),values=c("white","black"),guide="none") + 
  theme(legend.text=element_text(face="italic",size=12)) + coord_fixed(ratio=3/2,xlim=c(-121,-116),ylim=c(32.5,35))

basewidth=1.1*(125-113)
baseheight=1.1*(43-30)

insetwidth=1.1*(121-116)
insetheight=1.1*(35-32.5)
scalefactor=baseheight/insetheight

###combine and save
png("./Outputs/map_with_inset_plus.png",width=8,height=8,res=600,units="in")
ggarrange(base.plus,inset.plus,nrow=1,widths=c(basewidth,insetwidth*scalefactor))
dev.off()