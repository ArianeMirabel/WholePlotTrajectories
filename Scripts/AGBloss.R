#### get data damage for Tortue and Manaré, using Guyafor
#############################################################################

#  This calculation uses the methods of Camilia but with the allometric equation fo Rejou-Mechain 2017 for AGB
# and the tarif de cubage de l'Est for vol.ext
# The methods calculate the ABG before logging and for all inventory after logging 
# then dAGB = AGB initiale - AGB minimum i nthe 4 year after logging
# AGB damage = dAGB - AGB extracted
# uses a C content in AGB of 0.5
# and average wd bellow

### open libray and set values
library(data.table)
library(BIOMASS)


###connection odbc a Guyafor sur serveur SQL
library("RODBC")
Connex=odbcConnect(dsn="Guyafor")
####s?lection des donn?es
req1<-"SELECT 
dbo.TtGuyaforShiny.NomForet AS Forest, 
dbo.TtGuyaforShiny.n_parcelle AS Plot, 
dbo.TtGuyaforShiny.Surface AS PlotSurface, 
dbo.TtGuyaforShiny.n_carre AS SubPlot, 
dbo.TtGuyaforShiny.n_arbre AS TreeFieldNum, 
dbo.TtGuyaforShiny.i_arbre AS idTree, 
dbo.TtGuyaforShiny.X AS Xfield, 
dbo.TtGuyaforShiny.Y AS Yfield, 
dbo.TtGuyaforShiny.Xutm, 
dbo.TtGuyaforShiny.Yutm, 
dbo.TtGuyaforShiny.UTMZone, 
dbo.TtGuyaforShiny.Lat, 
dbo.TtGuyaforShiny.Lon, 
dbo.TtGuyaforShiny.Famille AS Family, 
dbo.TtGuyaforShiny.Genre AS Genus, 
dbo.TtGuyaforShiny.Espece AS Species, 
dbo.TtGuyaforShiny.SourceBota AS BotaSource, 
dbo.TtGuyaforShiny.indSurete AS BotaCertainty, 
dbo.TtGuyaforShiny.n_essence, 
dbo.TtGuyaforShiny.idPilote AS idVern, 
dbo.TtGuyaforShiny.nomPilote AS VernName, 
dbo.TtGuyaforShiny.Commerciale AS CommercialSp, 
dbo.TtGuyaforShiny.Densite AS WoodDensity, 
dbo.TtGuyaforShiny.campagne AS CensusYear, 
dbo.TtGuyaforShiny.DateMesure AS CensusDate, 
dbo.TtGuyaforShiny.code_vivant AS CodeAlive, 
dbo.TtGuyaforShiny.code_mesure AS CodeMeas, 
dbo.TtGuyaforShiny.circonf AS Circ, 
dbo.taMesure_Corr.circ_corr AS CircCorr, 
dbo.taMesure_Corr.code_corr AS CodeCorr
FROM dbo.TtGuyaforShiny LEFT OUTER JOIN dbo.taMesure_Corr ON dbo.TtGuyaforShiny.idMesure = dbo.taMesure_Corr.idMesure"
sqlQuery(Connex,req1)->DataGuyafor
####cloture de la connection odbc
odbcClose(Connex)



### house keeping
DataGuyafor=data.table(DataGuyafor) # transform as a data.table
Data=DataGuyafor[(Forest=="Montagne Tortue" & Plot %in% c("P17 exploitée 1", "P17 exploitée 2")) |
                   (Forest=="Régina St Georges" & Plot %in% c("Manaré I", "Manaré II")) |
                   (Forest=="Paracou" & Plot %in% c("2", "7", "9")) ] # keep only the plots I want
Data$Plot=as.factor(as.character(Data$Plot))
Data$DBH = Data$CircCorr/pi # calculate DBH
Data=Data[DBH>=10, ] # remove less than 10 cm DBH

wd=0.736 # average Wodd density for Vol.ext

# give info on year of ref (year before logging for which tk=0) and vol.ext
Yrref=data.table(Plot=c("P17 exploitée 1", "P17 exploitée 2", "2", "7", "9", "Manaré I" , "Manaré II"), 
                 t0=c(2005, 2005, 1986,1986,1986,2009, 2009))
Data=merge(Data,Yrref, by = "Plot", all = TRUE, sort=FALSE)
# add tk (time after logging)
Data$tk=Data[,CensusYear-t0]
# keep only year of interest (t0 and all t within 4 years after logging, so tk=1 to tk=5)
Data=Data[tk>=0 & tk<=5]


# for Manaré II some trees don't have Lat and Lon => replace them by the mean of all trees
LatManII=Data[Plot =="Manaré II", mean(Lat,na.rm = TRUE)]
LonManII=Data[Plot =="Manaré II", mean(Lon,na.rm = TRUE)]
Data[Plot=="Manaré II" & is.na(Lat), Lat:=LatManII]
Data[Plot=="Manaré II" & is.na(Lon), Lon:=LonManII]



###################################################################################################
# WITHout UNCERTAINTY
# with the value of wd retreived from BIOMASS for AGB and given for extracted wood
#  with 50% C in AGB

#Data=Data[,-"wd"]
Taxo=correctTaxo(genus=as.character(Data$Genus), species=as.character(Data$Species)) # correct the species names
WDdata=getWoodDensity(genus=Taxo$genusCorrected, species=Taxo$speciesCorrected, stand=Data$Plot) # get WD (mean and sd)
Data$meanWD=WDdata$meanWD
Data$sdWD=WDdata$sdWD

#### calculate AGB per plot and year (!!! given per ha, that's why divided by area) !!! only for trees that are alive
AGB=Data[CodeAlive==1,.(AGB=sum(computeAGB(D=DBH, WD=meanWD, coord = cbind(Lon,Lat))/PlotSurface)), by=.(Forest, Plot, PlotSurface, CensusYear, tk)]
setorder(AGB, Plot, tk)

##### calculate vol.ext using tarif cubage cf Guitet guide de sylviculture 2014 (!!! given per ha, that's why divided by area) all plot together
CubEst= function (DBH) {-0.084516 + 10.461316 * ((DBH/100)^2)}
CubCentre= function (DBH) {-0.035829 + 8.7634 * ((DBH/100)^2)}
Vol.extEst=Data[Forest %in% c("Montagne Tortue", "Régina St Georges") & CodeAlive==0 & CodeMeas==4, 
                .(vol.ext=unique(sum(CubEst(DBH))/PlotSurface)), by =.(Plot)]
Vol.extCentre=Data[Forest =="Paracou" & CodeAlive==0 & CodeMeas==4, 
                   .(vol.ext=unique(sum(CubCentre(DBH))/PlotSurface)), by =.(Plot)]
Vol.ext=rbind(Vol.extCentre, Vol.extEst)
AGB=merge(AGB,Vol.ext, by = "Plot", all = TRUE, sort=FALSE)

# dAGB (agbloss = agb0-agbmin in next 4 years) and add a colum giving the year for which AGB is min
damdata=AGB[, .(dAGB=AGB[1]-min(AGB[-1]), Yrmin=CensusYear[which.min(AGB[-1])+1]), by=.(Plot, Forest, PlotSurface, vol.ext)]
damdata$AGBext = damdata$vol.ext*wd # AGB extracted (logged)
damdata$AGBdam = damdata[,dAGB-AGBext] # AGB of damage
damdata$Cdam = damdata[,AGBdam/2] # C of damage (using C = 0.5 AGB)

write.csv(damdata, file="data/damdata_TortueManareParacou_1801_WD_specific.csv", row.names = FALSE)
###################################################################################################


###################################################################################################
# With UNCERTAINTY 
  # on allometric equation
  # on WD
  # on D (Dpropag="chave2004")
  # Set a fixed year for AGBmin, ingerited from the calculation without uncertainty 
    # => that is to make sure that the Year of min AGB is the same for all iteration (see email Camila 02/01/18)
  # two calculations : (1) with 50% C in AGB (fixed) and (2) with an uncertain C content in AGB (mean=0.4713, sd=0.0206)


# keep only data for t0 et tmin
dataYrmin=damdata[,.(Plot, Yrmin)]
Data=merge(Data,dataYrmin, by = "Plot", all = TRUE, sort=FALSE)
Data_minFix=Data[CensusYear==t0 | CensusYear==Yrmin]


######### (1) Do calculation without incertainty for the C content in AGB (taken as 0.5)
#### calculate 100 iteration of AGB  per plot and year (!!! given per ha, that's why divided by area) !!! only for trees that are alive
niter=500 # set number of iterations to calculate AGB0 and AGBtk
AGBdistrib=Data_minFix[CodeAlive==1,
                .(AGB=colSums(AGBmonteCarlo(D=DBH, WD=meanWD, coord = cbind(Lon,Lat),
                                            errWD=sdWD, n=niter,Dpropag="chave2004")$AGB_simu)/PlotSurface), 
                by=.(Forest, Plot, PlotSurface, CensusYear, tk)] # get niter iterations of AGB 
setorder(AGBdistrib, Plot, tk)
AGBdistrib$iter=seq.int(1:niter)
# to test results AGBdistrib similar to results for AGB => should be close to 1
#AGBbis=merge(AGB, dataYrmin, by = "Plot", all = TRUE, sort=FALSE)
#AGBbis=AGBbis[tk==0 | CensusYear==Yrmin]
#mean(AGBdistrib[,mean(AGB), by=.(Plot, CensusYear)]$V1 / AGBbis[,AGB, by=.(Plot, CensusYear)]$AGB)
AGBdistrib=merge(AGBdistrib,Vol.ext, by = "Plot", all = TRUE, sort=FALSE)

# dAGB (agbloss = agb0-agbmin in next 4 years)
damdataIter=AGBdistrib[, .(dAGB=AGB[1]-AGB[2]), by=.(Plot, Forest, iter, PlotSurface, vol.ext)]
damdataIter$AGBext = damdataIter$vol.ext*wd # AGB extracted (logged)
damdataIter$AGBdam = damdataIter[,dAGB-AGBext] # AGB of damage
damdataIter$Cdam = damdataIter[,AGBdam/2] # C of damage (using C = 0.5 AGB)
write.csv(damdataIter, file="01_data/damdataIter_TortueManareParacou_Cfixed.csv", row.names = FALSE)

# to test results damdataIter similar to results for damdata => should be close to 1
#mean(damdataIter[,mean(dAGB), by=.(Plot)]$V1 / damdata[,dAGB, by=.(Plot)]$dAGB)
#mean(damdataIter[,mean(AGBext), by=.(Plot)]$V1 / damdata[,AGBext, by=.(Plot)]$AGBext) # should be exactly 1 because no uncertainty
#mean(damdataIter[,mean(AGBdam), by=.(Plot)]$V1 / damdata[,AGBdam, by=.(Plot)]$AGBdam)
#mean(damdataIter[,mean(Cdam), by=.(Plot)]$V1 / damdata[,Cdam, by=.(Plot)]$Cdam) # should be the same than line above
                                                                                                                                                                                                                                                                                                                                                                                                                            

######### (2) Do calculation with incertainty for the C content in AGB (taken with mean=0.4713, sd=0.0206)
# ! these values are not the one in the BIOMASS paper but in the code of the function AGBmonteCarlo
# I don't use the function AGBmonteCarlo with Carbon=TRUE to calculate AGC to avoid having a different C content value for AGB0, AGBmin et Vol.ext
# instead I convent the dAGB et AGBext to C with a set value of C content per iteration
damdataIterC = damdataIter[,-(7:9)]
damdataIterC$Ccontent = rnorm(niter, mean=0.4713, sd=0.0206) # get a value of C content used for all value in the iteration
damdataIterC$dAGC = damdataIterC[,dAGB*Ccontent]
damdataIterC$Cext = damdataIterC[,vol.ext*wd*Ccontent]
damdataIterC$Cdam = damdataIterC[,dAGC-Cext]
write.csv(damdataIterC, file="01_data/damdataIter_TortueManareParacou_Cuncertain.csv", row.names = FALSE)

# to test results damdataIter similar to results for damdata => should be close to 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       # to test results damdataIterC similar to results for damdata => should be close to 1
#mean(damdataIterC[,mean(dAGC), by=.(Plot)]$V1 / damdata[,dAGB*0.4713, by=.(Plot)]$V1)
#mean(damdataIterC[,mean(Cdam), by=.(Plot)]$V1 / damdata[,AGBdam*0.4713, by=.(Plot)]$V1) 

