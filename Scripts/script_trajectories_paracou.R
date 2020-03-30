###########################################################
### Trajectories of recovery in Paracou
###########################################################

# trajectories of the 12 plots of disturbance experiment


library(data.table)
library(ggplot2)
library(BIOMASS)
library(EcoFoG)

load("DB/AGBlostDB")
DataParacou <- AllYears 

###### get data and tidying up ######
######  
DataParacou <- Guyafor2df(UID="Geraldine.Derroire",PWD="geraldine973+", WHERE="Forest='Paracou'")
# DataParacou <- Guyafor2df(UID="Geraldine.Derroire",PWD="geraldine973+", WHERE="Forest='Paracou'AND Plot='16'")
# Data <- as.data.table(read.csv("190310_DataParacou.csv", header=TRUE))
Data <- as.data.table(DataParacou)
# calculate DBH
Data$DBH <- Data$CircCorr / pi
# keep only DBH >=0 (a few tree have a circ-corr < 0, mostly Bois CathÃ©drale, cf code correction Camilia)
Data <- Data[DBH >= 10]

# keep only plot 1 to 12
# Data <- Data[Plot %in% 1:12]
# table(Data$Plot, Data$CensusYear)

# add treatment
Data$Treat <- as.factor(NA)
Data[n_parcelle %in% c("1", "6", "11"), Treat := "Control"]
Data[n_parcelle %in% c("2", "7", "9"), Treat := "T1"]
Data[n_parcelle %in% c("3", "5", "10"), Treat := "T2"]
Data[n_parcelle %in% c("4", "8", "12"), Treat := "T3"]
Data$Treat <- as.factor(Data$Treat)
# table(Data$Plot, Data$Treat)


###### set appearance of graphs ######
######  
# colour of treatments
MyColTreat <-c("Control" = "yellowgreen" ,"T1" = "gold",
               "T2" = "orangered2", "T3" = "tomato4") 


###### trajectories of AGB ######
######  
# get WD 
Taxo=correctTaxo(genus=as.character(Data$Genre), species=as.character(Data$Espece)) # correct the species names
WDdata=getWoodDensity(genus=Taxo$genusCorrected, species=Taxo$speciesCorrected, stand=Data$n_parcelle) # get WD (mean and sd)
Data$meanWD=WDdata$meanWD
Data$sdWD=WDdata$sdWD
#### calculate AGB per plot and year (!!! given per ha) (remove dead trees) without uncertainty

Data<-subset(Data,campagne==1984)

AGB=Data[,.(AGB=sum(computeAGB(D=circonf, WD=meanWD, coord = cbind(Lon,Lat)))), 
         by=.(n_parcelle, Treat)]
#### plot AGB trajectories
PlotAGB <- ggplot(data = AGB) + 
  geom_line(mapping = aes(x = CensusYear, y = AGB, group=Plot, 
                            colour=Treat), size=1.2) +
  geom_vline(aes(xintercept=1986.5), linetype="dashed") +
  facet_wrap( ~ Treat) + 
  scale_colour_manual(values = MyColTreat) +
  labs(x="",y="AGB (Mg/ha)") +
  theme_bw() + theme(legend.position="none")



####### AGB loss per plot 
# dAGB (agbloss = agb0-agbmin in next 4 years) 
AGBlogged <- AGB[Treat %in% c("T1", "T2", "T3")]
AGBpostlogged <- AGBlogged[CensusYear %in% 1986:1991] 
setorder(AGBpostlogged, Plot, CensusYear)

AGBloss <- AGBpostlogged[, .(AGBloss=(AGB[1]-min(AGB[-1]))*100/AGB[1]),by=.(Plot, Treat)]




############A verifier ##############

#### Faire une fonction pour les plots

# stem density per ha
Stemdens <- Data[CodeAlive==1,.(StemDens=.N/6.25), by=.(Plot, CensusYear, Treat)]
# Stendens[,.(meanDens=mean(V1), sdDens=sd(V1)), by=.(CensusYear, Treat)]
#### plot stem density trajectories
PlotStemDens <- ggplot(data = Stemdens) + 
  geom_line(mapping = aes(x = CensusYear, y = StemDens, group=Plot, 
                          colour=Treat), size=1.2) +
  geom_vline(aes(xintercept=1986.5), linetype="dashed") +
  facet_wrap( ~ Treat) + 
  scale_colour_manual(values = MyColTreat) +
  labs(x="",y="Stem density (per ha)") +
  theme_bw() + theme(legend.position="none")


# Basal area per ha
BA <- Data[CodeAlive==1,.(BA=sum(pi*(DBH/200)^2/6.25)), 
           by=.(Plot, CensusYear, Treat)]
#### plot BA trajectories
PlotBA <- ggplot(data = BA) + 
  geom_line(mapping = aes(x = CensusYear, y = BA, group=Plot, 
                          colour=Treat), size=1.2) +
  geom_vline(aes(xintercept=1986.5), linetype="dashed") +
  facet_wrap( ~ Treat) + 
  scale_colour_manual(values = MyColTreat) +
  labs(x="",y="Basal area (m2/per ha)") +
  theme_bw() + theme(legend.position="none")
