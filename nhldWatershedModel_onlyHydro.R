##### just hydrologic model for now...

# clear workspace
rm(list=ls())

setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/NHLD_lakeCmodel")

########## load utility functions and packages
source("nhldWatershedModel_onlyHydroSupporting.R")
require(deSolve)
require(LakeMetabolizer)

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

# starting year/month/day, ending year/mnth/day, & set up force/flux
startYear=2012
startMonth=1
startDay=1

endYear=2013
endMonth=12
endDay=31

curForce=force[which(((force[,3]==startYear) & (force[,4]==startMonth) & (force[,5]==startDay) & (force[,6]==0))):which(((force[,3]==endYear) & (force[,4]==endMonth) & (force[,5]==endDay) & (force[,6]==23))),]

curFlux=flux[which(((flux[,3]==startYear) & (flux[,4]==startMonth) & (flux[,5]==startDay))):which(((flux[,3]==endYear) & (flux[,4]==endMonth) & (flux[,5]==endDay))),]

curFluxDOY=1:nrow(curFlux)
curForceDOY=rep(curFluxDOY,each=24)

setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/NHLD_lakeCmodel")

###### forcing functions
dailyPrecip=approxfun(curFluxDOY,curFlux[,7],method="constant")
dailyEvap=approxfun(curFluxDOY,curFlux[,8],method="constant")
dailyRunoff=approxfun(curFluxDOY,curFlux[,6],method="constant")

####### ice cover: using regressions from Finland ()
# -> alternative: try looking at NTL data...; fit model with degree days?

###***** this stuff will have troubles if ice on happens after january 1 ******###
aprilForce=force[force[,4]==4,]
novForce=force[force[,4]==11,]

aprilAT=tapply(aprilForce[,11],aprilForce[,3],FUN=mean)
novAT=tapply(novForce[,11],novForce[3],FUN=mean)

iceOnDOY=(novAT+28.8)/0.19+182	#equation for November from inland North -> has latest ice on; pick this one because study area much higher N

iceOffDOY=(aprilAT-49.7)/-0.16+182-365	#Coastal/South model for breakup; again earliest date, but want this bias...

iceON=numeric(nrow(curFlux))
years=sort(unique(curFlux[,3]))
for(i in 1:length(years)){
	curYear=curFlux[curFlux[,3]==years[i],]
	curDOY=as.numeric(strftime(strptime(paste(curYear[,3],curYear[,4],curYear[,5],sep="-"),format="%Y-%m-%d"),format="%j"))
	iceON[curFlux[,3]==years[i]]=((curDOY<iceOffDOY[names(iceOffDOY)==years[i]]) | (curDOY>iceOnDOY[names(iceOnDOY)==years[i]]))*1
}

####### load lake/watershed information
UNDERCsheds=read.table("../NHLDwatershedDelineations/UNDERCsheds_4-24-15.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)

GFLOWoutput=read.table("../gflowOutput_3-24-15/GFLOWperElementDischarge_4-24-15.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)

###### iterate through lakes  
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
	
	hydroSumm=cbind(times,dailyRunoff(times)*curShedArea/1000,dailyPrecip(times)*curLakeArea/1000,rep(gwIn,length(times)),QoutSim,dailyEvap(times)*curLakeArea/1000*(1-iceON),rep(gwOut,length(times)))
	colnames(hydroSumm)=c("time","Qin","precip","GWin","Qout","evap","GWout")

	dev.new()
	par(mfrow=c(3,3))
	plot(hydroSumm[,1],hydroSumm[,2],xlab="time",ylab="Qin (m3 d-1)",type='l',lwd=1)
	plot(hydroSumm[,1],hydroSumm[,3],xlab="time",ylab="Precip (m3 d-1)",type='l',lwd=1)
	plot(hydroSumm[,1],hydroSumm[,4],xlab="time",ylab="GWin (m3 d-1)",type='l',lwd=1)
	plot(hydroSumm[,1],hydroSumm[,5],xlab="time",ylab="Qout (m3 d-1)",type='l',lwd=1)
	plot(hydroSumm[,1],hydroSumm[,6],xlab="time",ylab="Evap (m3 d-1)",type='l',lwd=1)
	plot(hydroSumm[,1],hydroSumm[,7],xlab="time",ylab="GWout (m3 d-1)",type='l',lwd=1)
	plot(out[,1],out[,2],xlab="time",ylab="Volume (m^3)",type='l',lwd=1)
	plot(out[,1],stageSim,xlab="time",ylab="stage (m)",type='l',lwd=1)
	
	meanBudget=colMeans(hydroSumm[,-1],na.rm=TRUE)
	
	barplot(meanBudget*c(1,1,1,-1,-1,-1),names.arg=names(meanBudget))
	
	avgDailyIn=sum(meanBudget[1:3])
	avgDailyOut=sum(meanBudget[4:6])
	# estimate of residence time
	print(mean(out[,2],na.rm=TRUE)/(mean(c(avgDailyIn,avgDailyOut))))
	
annualMin=tapply(stageSim,curFlux[,3],FUN=min,na.rm=TRUE)
dev.new()
plot(as.numeric(names(annualMin)),annualMin-mean(annualMin,na.rm=TRUE),type='o',xlab="year",ylab="min lake level anomaly")
abline(h=0,lty=2)