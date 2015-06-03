##### carbon dynamics added to pure hydrologic model

# clear workspace
rm(list=ls())

setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/NHLD_lakeCmodel")

########## load utility functions and packages
source("nhldWatershedModel_ConeSupporting.R")
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
startYear=2005
startMonth=5
startDay=1

endYear=2013
endMonth=9
endDay=30

curForce=force[which(((force[,3]==startYear) & (force[,4]==startMonth) & (force[,5]==startDay) & (force[,6]==0))):which(((force[,3]==endYear) & (force[,4]==endMonth) & (force[,5]==endDay) & (force[,6]==23))),]

curFlux=flux[which(((flux[,3]==startYear) & (flux[,4]==startMonth) & (flux[,5]==startDay))):which(((flux[,3]==endYear) & (flux[,4]==endMonth) & (flux[,5]==endDay))),]

curFluxDOY=1:nrow(curFlux)
curForceDOY=rep(curFluxDOY,each=24)

setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/NHLD_lakeCmodel")

###### forcing functions
dailyPrecip=approxfun(curFluxDOY,curFlux[,7],method="constant")
dailyEvap=approxfun(curFluxDOY,curFlux[,8],method="constant")
dailyWind=approxfun(curFluxDOY,tapply(curForce[,13],curForceDOY,FUN=mean),method="constant")

var_dailyRunoff=approxfun(curFluxDOY,5+curFlux[,6],method="constant")

const_dailyRunoff=approxfun(curFluxDOY,rep(sum(curFlux[,6])/length(curFluxDOY),length(curFluxDOY)),method="constant")

#dailySun<-function(day){
#	todayHourlySun=sw.to.par.base(curForce[curForceDOY==day,10])		# umol m2 sec; SW to PAR based on Read...
#	return(todayHourlySun)
#}
#dailyMaxWind<-function(day){
#	todayMaxWind=max(curForce[curForceDOY==day,13])*(1/10)^0.15		# convert 10m wind to 1m for preston model tPOC deposition
#	return(todayMaxWind)
#}

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
for(j in 1:length(years)){
	curYear=curFlux[curFlux[,3]==years[j],]
	curDOY=as.numeric(strftime(strptime(paste(curYear[,3],curYear[,4],curYear[,5],sep="-"),format="%Y-%m-%d"),format="%j"))
	iceON[curFlux[,3]==years[j]]=((curDOY<iceOffDOY[names(iceOffDOY)==years[j]]) | (curDOY>iceOnDOY[names(iceOnDOY)==years[j]]))*1
}

####### load lake/watershed information
UNDERCsheds=read.table("../NHLDwatershedDelineations/NHLDsheds_5-13-15.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)

#GFLOWoutput=read.table("../gflowOutput/gflowSensitivity_5-21-15/NHLD_GFLOW_UNDERC_20150520_Simplify010_RegK1_LocK0pt01_OutputLong.txt",header=TRUE,stringsAsFactors=FALSE)
GFLOWoutput=read.table("../gflowOutput/NHLDgflow_060115/NHLD_GFLOW_ALL_OutputLong.txt",header=TRUE,stringsAsFactors=FALSE)
colnames(GFLOWoutput)[1]="Permanent_"

#i=17		#Long Lake from UNDERCsheds
#i=4		#Morris from UNDERCsheds

##****** i = 214

lakes=sort(unique(GFLOWoutput[,1]))
lakes=lakes[lakes%in%UNDERCsheds$Permanent_]

for(i in 1:length(lakes)){
	print(i)
curLakeID=lakes[i]#UNDERCsheds$Permanent_[i]
	
# current lake and shed parameters
A0=UNDERCsheds$NHLD_lakes[UNDERCsheds$Permanent_==curLakeID]		#m2
V0=10^(-0.0589+1.12963*log10(A0))		#m3

# try volume model from del Giorgio group now that we have difference in elevation data

Perim0=UNDERCsheds$Perimeter[UNDERCsheds$Permanent_==curLakeID]			
DL=Perim0/(2*sqrt(pi*A0))

r0=sqrt(A0/pi)
stage0=V0/(pi*r0^2/3)
r2h=r0/stage0

alpha=0.99
stageOut=alpha*stage0

curShedArea=UNDERCsheds$Area_m2[UNDERCsheds$Permanent_==curLakeID]
	
curGFLOW=GFLOWoutput[GFLOWoutput$Permanent_==curLakeID,]
curGFLOW=curGFLOW[!is.na(curGFLOW$Permanent_),]
GFLOWpropPerim=as.numeric(curGFLOW$Linesink_Length_Output)/sum(as.numeric(curGFLOW$Linesink_Length_Output))
GFLOWin=(curGFLOW$Linesink_PerLengthDischarge_Output>0)*1
	
gwIn0=sum(GFLOWin*as.numeric(curGFLOW$Linesink_PerLengthDischarge_Output)*GFLOWpropPerim*Perim0)*0.0283168	#m3 d-1
gwOut0=-sum((1-GFLOWin)*as.numeric(curGFLOW$Linesink_PerLengthDischarge_Output)*GFLOWpropPerim*Perim0)*0.0283168	#m3 d-1
			

# calculate concentrations of DIC, DOC, TP based on landcover
#setwd("/Files/NDongoing/ScalingLakeProcessesVilasCty/vilasWatersheds_post5-17-14/waterConstituents")
#ntlGW=read.csv("NTLlter__groundwater_chemistry.csv",header=TRUE,fill=TRUE)
#ntlGW=ntlGW[,c(1:4,9:10,14)]
#ntlGW=ntlGW[!is.na(ntlGW[,5]),]
#ntlGW=ntlGW[ntlGW$doc<10,]
#lottig2011DOC=read.table("lottig2012Fig.txt",header=FALSE,sep=",")
kH=29.41	# Henry's Law constant for CO2 in water [L atm mol-1]	# temperature sensitivity
streamDOC=exp(1.3961+3.245*(UNDERCsheds$percentWetland[UNDERCsheds$Permanent_==curLakeID]/100))/12*1000/1000	# from lottig 2012; mol m-3
streamPOC=3/12*1000/1000		# ~3 mg L-1; buffam 2011; mol m-3
streamDIC=10/12*1000/1000		#10 mg L-1; lottig 2011; mol m-3
streamP=0.05/31*1000/1000#0.084/31*1000/1000#0.04/31*1000/1000	  # Long inlet has 84 ug/L	#.025 mg L-1 TDP & 0.04 mg L-1 TP; lottig 2011; mol m-3
gwDOC=13/12*1000/1000#median(ntlGW$doc,na.rm=TRUE)/12*1000/1000	# mol m-3
gwDIC=0.7025#median(ntlGW$dic,na.rm=TRUE)/12*1000/1000	# mol m-3
gwP=0.0007742#median(ntlGW$totp,na.rm=TRUE)/31*1000/1e6		# mol m-3
precipDOC=1.1/12*1000/1000		#mol m-3; from Likens et al. 1983 @ Hubbard Brook
precipDIC=400/1e6*1/kH*1000	# mol C m-3; assumed in equilibrium with atmosphere
precipP=0.01/31*1000/1000	# mol m-3; Murphy & DOskey 1976 JGLR


params=c(curShedArea=curShedArea,stageOut=stageOut,Perim0=Perim0,gwIn0=gwIn0,gwOut0=gwOut0,DL=DL,precipDIC=precipDIC,precipDOC=precipDOC,streamDIC=streamDIC,streamDOC=streamDOC,gwDIC=gwDIC,gwDOC=gwDOC,kH=kH,r2h=r2h)

initialX=c(V=V0,DIC=0.3/12*V0,DOC=8/12*V0,tPOC=0.8/12*V0,phyto=0.02/12*V0,P=0.5/1000/31*V0,Emit=0,Sed_tPOC=0,Sed_phyto=0)

times=curFluxDOY

dailyRunoff=var_dailyRunoff

# if maximum daily surface inflow is of particular values (based on UNDERC lakes) simulate accordingly
if(20*curShedArea/1000/V0<0.12){
	out<-ode(y=initialX,times=times,func=timeStep,parms=params,method="rk4")
	write.table(out,paste("V",curLakeID,".txt",sep=""),row.names=FALSE,sep="\t")
}else if(20*curShedArea/1000/V0>0.16){
	dailyRunoff=const_dailyRunoff
	out<-ode(y=initialX,times=times,func=timeStep,parms=params,method="rk4")
	write.table(out,paste("C",curLakeID,".txt",sep=""),row.names=FALSE,sep="\t")
}else{
	possibleError<-tryCatch(ode(y=initialX,times=times,func=timeStep,parms=params,method="rk4"),error=function(e) e)
	if(inherits(possibleError,"error")){
		dailyRunoff=const_dailyRunoff
		out<-ode(y=initialX,times=times,func=timeStep,parms=params,method="rk4")
		write.table(out,paste("C",curLakeID,".txt",sep=""),row.names=FALSE,sep="\t")	
	}else{
		out<-ode(y=initialX,times=times,func=timeStep,parms=params,method="rk4")
		write.table(out,paste("V",curLakeID,".txt",sep=""),row.names=FALSE,sep="\t")
	}
}
}