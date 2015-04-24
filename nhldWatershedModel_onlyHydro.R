##### just hydrologic model for now...

rm(list=ls())

setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/NHLD_lakeCmodel")

########## load utility functions
source("nhldWatershedModel_onlyHydroSupporting.R")
require(deSolve)
require(LakeMetabolizer)

# load watershed table
#d=read.table("vilasSheds_5-18-14.txt",sep="\t",header=TRUE)

#W2L=d$shedArea_m2/d$lakeArea_m2
#toFix=which(W2L>100 & d$lakeArea_m2<10000)
#d$shedArea_m2[toFix]=d$lakeArea_m2[toFix]*1.5

#W2L=d$shedArea_m2/d$lakeArea_m2

#RTs=d$lakeVol_m3/(dailyP/1000*d$shedAreaNoLake_m2-d$shedAreaNoLake_m2*shedETs/1000+dailyP/1000*d$lakeArea_m2)


######## load forcing and flux data
setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/VICforcings&output/PROC")

cells=gsub("FLUX_","",grep("FLUX",list.files(),value=TRUE))

forceFiles=grep("FORCE",list.files(),value=TRUE)
fluxFiles=grep("FLUX",list.files(),value=TRUE)

#### UNDERC includes four VIC cells
# 46.28125_-89.53125
# 46.28125_-89.46875
# 46.21875_-89.53125
# 46.21875_-89.46875
UNDERCcells=cells[c(169,170,182,183)]

# Lat Long Year Month Day Runoff_mmDay Precip_mmDay EvapOpenWater_mmDay
flux=read.table(paste("FLUX_",UNDERCcells[1],sep=""),header=FALSE)
# Lat Long Year Month Day Hour NetLW_Wm2 NetSW_Wm2 LW_Wm2 SW_Wm2 AirTemp_dC atmPress_kPa windspeed_mS@10m relHumid_fraction
force=read.table(paste("FORCE_",UNDERCcells[1],sep=""),header=FALSE)

curForce=force[force[,3]==2012,]
curForce=curForce[curForce[,4]%in%(5:9),]
curFlux=flux[flux[,3]==2012,]
curFlux=curFlux[curFlux[,4]%in%(5:9),]

curForce=force[force[,3]==2013,]
curForce=curForce[curForce[,4]%in%(5:9),]
curFlux=flux[flux[,3]==2013,]
curFlux=curFlux[curFlux[,4]%in%(5:9),]

curFluxDOY=as.numeric(strftime(strptime(paste(curFlux[,4],curFlux[,5],curFlux[,3],sep="-"),format="%m-%d-%Y"),format="%j"))
curForceDOY=as.numeric(strftime(strptime(paste(curForce[,4],curForce[,5],curForce[,3],sep="-"),format="%m-%d-%Y"),format="%j"))

curFluxDOY=curFluxDOY-min(curFluxDOY)+1
curForceDOY=curForceDOY-min(curForceDOY)+1

setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/NHLD_lakeCmodel")
###### forcing functions --> when running longer can just use January 1, 1950 as an epoch
dailyPrecip=approxfun(curFluxDOY,curFlux[,7],method="constant")
dailyEvap=approxfun(curFluxDOY,curFlux[,8],method="constant")
dailyRunoff=approxfun(curFluxDOY,curFlux[,6],method="constant")




####### seems somewhat reasonable for summertime run

#### to run longer term need to consider:
# - ice cover:  shut off evap?  base on cumulative air temp in some way... --> look at lit.!
# - winter time precip:  does VIC discriminate between snow and not?


UNDERCsheds=read.table("../NHLDwatershedDelineations/UNDERCsheds_4-24-15.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)

GFLOWoutput=read.table("../gflowOutput_3-24-15/GFLOWperElementDischarge_4-24-15.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)

# iterate through lakes  --> perhaps in future run all lakes at once????
#for i in 1:nrow(d){
	i=17		#Long Lake from UNDERCsheds
	i=4		#Morris from UNDERCsheds
	
	curLakeID=UNDERCsheds$Permanent_[i]
	
	# current lake and shed parameters
	curLakeArea=UNDERCsheds$NHLD_lakes[i]		#m2
	V0=exp(-0.01857+1.11563*log(curLakeArea))#322027.2					#d$lakeVol_m3[i]		#m3
	curLakePerim=UNDERCsheds$Perimeter[i]
	##### when don't have volume infer from area...
	#### for Vilas Cty:  log(volume)=-0.01857+1.11563*log(area)	m^3 and m^2 for units
	curShedArea=UNDERCsheds$Area_m2[i]	#m2
	
	
	#gwIn and gwOut
	curGFLOW=GFLOWoutput[GFLOWoutput$Permanent_==curLakeID,]
	curGFLOW=curGFLOW[!is.na(curGFLOW$Permanent_),]
	GFLOWpropPerim=curGFLOW$Linesink_Length_Output/sum(curGFLOW$Linesink_Length_Output)
	GFLOWin=(curGFLOW$Linesink_PerLengthDischarge_Output>0)*1
	gwIn=sum(GFLOWin*curGFLOW$Linesink_PerLengthDischarge_Output*GFLOWpropPerim*curLakePerim)*0.0283168	#m3 d-1
	gwOut=sum((1-GFLOWin)*curGFLOW$Linesink_PerLengthDischarge_Output*GFLOWpropPerim*curLakePerim)*0.0283168	#m3 d-1
		
	#### these seem too high; for now just divide GW by 60 to make some progress...
	gwIn=gwIn/60
	gwOut=gwOut/60
	
	stage0=V0/curLakeArea
	alpha=0.99
	stageOut=alpha*stage0

	params=c(curLakeArea=curLakeArea,curShedArea=curShedArea,stageOut=stageOut,gwIn=gwIn,gwOut=gwOut)

	initialX=c(V=V0)
	
	times=curFluxDOY
	
	out<-ode(y=initialX,times=times,func=timeStep,parms=params)

	stageSim=out[,2]/curLakeArea
	C=(2/3)^1.5*9.806^0.5	# m s-2
	L=0.5	# m
	H=stageSim-stageOut
	QoutSim=ifelse(H>0,(C*L*H^1.5)*(60*60*24),0)
	
	hydroSumm=cbind(times,dailyRunoff(times)*curShedArea/1000,dailyPrecip(times)*curLakeArea/1000,rep(gwIn,length(times)),QoutSim,dailyEvap(times)*curLakeArea/1000,rep(gwOut,length(times)))
	colnames(hydroSumm)=c("time","Qin","precip","GWin","Qout","evap","GWout")

	dev.new()
	par(mfrow=c(3,3))
	plot(hydroSumm[,1],hydroSumm[,2],xlab="time",ylab="Qin (m3 d-1)",type='l',lwd=2)
	plot(hydroSumm[,1],hydroSumm[,3],xlab="time",ylab="Precip (m3 d-1)",type='l',lwd=2)
	plot(hydroSumm[,1],hydroSumm[,4],xlab="time",ylab="GWin (m3 d-1)",type='l',lwd=2)
	plot(hydroSumm[,1],hydroSumm[,5],xlab="time",ylab="Qout (m3 d-1)",type='l',lwd=2)
	plot(hydroSumm[,1],hydroSumm[,6],xlab="time",ylab="Evap (m3 d-1)",type='l',lwd=2)
	plot(hydroSumm[,1],hydroSumm[,7],xlab="time",ylab="GWout (m3 d-1)",type='l',lwd=2)
	plot(out[,1],out[,2],xlab="time",ylab="Volume (m^3)",type='l',lwd=2)
	plot(out[,1],stageSim,xlab="time",ylab="stage (m)",type='l',lwd=2)
	
	meanBudget=colMeans(hydroSumm[,-1],na.rm=TRUE)
	
	barplot(meanBudget*c(1,1,1,-1,-1,-1),names.arg=names(meanBudget))
	
