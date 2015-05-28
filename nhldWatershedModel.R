##### carbon dynamics added to pure hydrologic model

# clear workspace
rm(list=ls())

setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/NHLD_lakeCmodel")

########## load utility functions and packages
source("nhldWatershedModelSupporting.R")
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
startYear=1985
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
dailySun<-function(day){
	todayHourlySun=sw.to.par.base(curForce[curForceDOY==day,10])		# umol m2 sec; SW to PAR based on Read...
	return(todayHourlySun)
}
dailyMaxWind<-function(day){
	todayMaxWind=max(curForce[curForceDOY==day,13])*(1/10)^0.15		# convert 10m wind to 1m for preston model tPOC deposition
	return(todayMaxWind)
}

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
UNDERCsheds=read.table("../NHLDwatershedDelineations/UNDERCsheds_5-13-15.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)

GFLOWoutput=read.table("../gflowOutput/gflowSensitivity_5-21-15/NHLD_GFLOW_UNDERC_20150520_Simplify010_RegK1_LocK0pt01_OutputLong.txt",header=TRUE,stringsAsFactors=FALSE)
colnames(GFLOWoutput)[1]="Permanent_"

#i=17		#Long Lake from UNDERCsheds
#i=4		#Morris from UNDERCsheds
	
curLakeID=UNDERCsheds$Permanent_[i]
	
# current lake and shed parameters
A0=UNDERCsheds$NHLD_lakes[UNDERCsheds$Permanent_==curLakeID]		#m2
V0=10^(-0.0589+1.12963*log10(A0))		#m3

# try volume model from del Giorgio group now that we have difference in elevation data

Perim0=UNDERCsheds$Perimeter[UNDERCsheds$Permanent_==curLakeID]			
DR=0.45		# going with quadratic because it is simpler and close to 0.5 the two classes are almost identical 
p=6*DR-3
zbar0=V0/A0
zmax0=zbar0/DR
	
DL=Perim0/(2*sqrt(pi*A0))
	
# quadratic function peaks and begins to fall again at u=1...
#		problem for when lake rises above zmax0
#		one fix is to set zmax to something higher than value inferred from zbar and DR...
#		have to scale both zmax and zbar to maintain DR
#		after this, I think we can solve for Amax and Vmax and use these to scale Au and Vu appropriately...
#		5-12-15:  in trying UNDERC lakes Tenderfoot "overflows"  can bump this multiplier from 1.25 to 2.5 -> this changes stage dynamics and model behavior a bit because stage and Area at a given Volume are slightly different
zbar1=zbar0*2.5
zmax1=zmax0*2.5
	
uScale=zmax0/zmax1
					
A1=A0/(p*uScale^2+(1-p)*uScale)
V1=V0/((6*uScale-3*(1-p)*uScale^2-2*p*uScale^3)/(3+p))
					
curShedArea=UNDERCsheds$Area_m2[UNDERCsheds$Permanent_==curLakeID]	#m2

u0=round(uniroot(f=findU,lower=0,upper=1,p=p,Vmax=V1,Vu=V0)$root,4)
	
curGFLOW=GFLOWoutput[GFLOWoutput$Permanent_==curLakeID,]
curGFLOW=curGFLOW[!is.na(curGFLOW$Permanent_),]
GFLOWpropPerim=curGFLOW$Linesink_Length_Output/sum(curGFLOW$Linesink_Length_Output)
GFLOWin=(curGFLOW$Linesink_PerLengthDischarge_Output>0)*1
	
gwIn0=sum(GFLOWin*curGFLOW$Linesink_PerLengthDischarge_Output*GFLOWpropPerim*Perim0)*0.0283168	#m3 d-1
gwOut0=-sum((1-GFLOWin)*curGFLOW$Linesink_PerLengthDischarge_Output*GFLOWpropPerim*Perim0)*0.0283168	#m3 d-1
			
stage0=u0*zmax1
alpha=1.05#0.8
stageOut=alpha*stage0

# calculate concentrations of DIC, DOC, TP based on landcover
#setwd("/Files/NDongoing/ScalingLakeProcessesVilasCty/vilasWatersheds_post5-17-14/waterConstituents")
ntlGW=read.csv("NTLlter__groundwater_chemistry.csv",header=TRUE,fill=TRUE)
ntlGW=ntlGW[,c(1:4,9:10,14)]
ntlGW=ntlGW[!is.na(ntlGW[,5]),]
ntlGW=ntlGW[ntlGW$doc<10,]
lottig2011DOC=read.table("lottig2012Fig.txt",header=FALSE,sep=",")
kH=29.41	# Henry's Law constant for CO2 in water [L atm mol-1]	# temperature sensitivity
streamDOC=exp(1.3961+3.245*(UNDERCsheds$percentWetland[i]/100))/12*1000/1000	# from lottig 2012; mol m-3
#streamPOC=3/12*1000/1000		# ~3 mg L-1; buffam 2011; mol m-3
streamDIC=10/12*1000/1000		#10 mg L-1; lottig 2011; mol m-3
#streamP=0.04/31*1000/1000		#.025 mg L-1 TDP & 0.04 mg L-1 TP; lottig 2011; mol m-3
gwDOC=median(ntlGW$doc,na.rm=TRUE)/12*1000/1000	# mol m-3
gwDIC=median(ntlGW$dic,na.rm=TRUE)/12*1000/1000	# mol m-3
#gwP=median(ntlGW$totp,na.rm=TRUE)/31*1000/1e6		# mol m-3
precipDOC=1.1/12*1000/1000		#mol m-3; from Likens et al. 1983 @ Hubbard Brook
precipDIC=400/1e6*1/kH*1000	# mol C m-3; assumed in equilibrium with atmosphere
#precipP=0.01/31*1000/1000	# mol m-3; Murphy & DOskey 1976 JGLR
	
#DR=((curLakeArea*10^-6)^0.5)/curMeanDepth	# dynamic ratio; unitless; from Hakanson & Bryhn 2008. Water Air Soil Pollution, 187: 119-147. impacts turbulence and sedimentation
#distShore=sqrt(curLakeArea/pi)/2		#mean distance from shore; m; from GIS
#if(DR<0.26){		#adjustment to sedimentation rate based on turbulence, [unitless]; Hakanson & Bryhn 2008. Water Air Soil Pollution, 187: 119-147

#	Ydr=DR/0.26
#}else{
#	Ydr=0.26/DR
#} 	
# need a hypovolume initially to start simulation
#hypoVolume=0


params=c(Vmax=V1,Zmax=zmax1,Amax=A1,curShedArea=curShedArea,stageOut=stageOut,Perim0=Perim0,gwIn0=gwIn0,gwOut0=gwOut0,p=p,DL=DL,precipDIC=precipDIC,precipDOC=precipDOC,streamDIC=streamDIC,streamDOC=streamDOC,gwDIC=gwDIC,gwDOC=gwDOC,kH=kH)

initialX=c(V=V0,DIC=0.3/12*V0,DOC=8/12*V0)

times=curFluxDOY
	
out<-ode(y=initialX,times=times,func=timeStep,parms=params)
out=out[-(1:5),]		# always drops a lot over first few days
out=out[-nrow(out),]		# last time step always NA

outHydro=cbind(out[,1:2],NA,NA,NA)
for(j in 1:nrow(outHydro)){
	if(!is.na(outHydro[j,2])){
		uj=round(uniroot(f=findU,lower=0,upper=1,p=p,Vmax=V1,Vu=outHydro[j,2])$root,4)
		outHydro[j,3]=A1*(p*uj^2+(1-p)*uj)
		outHydro[j,4]=zmax1*uj
		outHydro[j,5]=2*pi*sqrt(outHydro[j,3]/pi)*DL
	}
}
colnames(outHydro)[3:5]=c("A","stage","perim")

C=(2/3)^1.5*9.806^0.5	# m s-2
L=0.1	# m
H=outHydro[,4]-stageOut
QoutSim=ifelse(H>0,(C*L*H^1.5)*(60*60*24),0)
	
hydroSumm=cbind(outHydro[,1],dailyRunoff(outHydro[,1])*(curShedArea-outHydro[,3])/1000,dailyPrecip(outHydro[,1])*outHydro[,3]/1000,gwIn0*outHydro[,5]/Perim0,QoutSim,dailyEvap(outHydro[,1])*outHydro[,3]/1000*(1-iceON[1:nrow(outHydro)]),gwOut0*outHydro[,5]/Perim0)
colnames(hydroSumm)=c("time","Qin","precip","GWin","Qout","evap","GWout")

startRow=out[1,1]
endRow=out[nrow(out),1]
plotDates=strptime(paste(curFlux[startRow:endRow,3],curFlux[startRow:endRow,4],curFlux[startRow:endRow,5],sep="-"),format="%Y-%m-%d")

dev.new()
par(mfrow=c(3,3))
plot(plotDates,hydroSumm[,2],xlab="time",ylab="Qin (m3 d-1)",type='l',lwd=1)
plot(plotDates,hydroSumm[,3],xlab="time",ylab="Precip (m3 d-1)",type='l',lwd=1)
plot(plotDates,hydroSumm[,4],xlab="time",ylab="GWin (m3 d-1)",type='l',lwd=1)
plot(plotDates,hydroSumm[,5],xlab="time",ylab="Qout (m3 d-1)",type='l',lwd=1)
plot(plotDates,hydroSumm[,6],xlab="time",ylab="Evap (m3 d-1)",type='l',lwd=1)
plot(plotDates,hydroSumm[,7],xlab="time",ylab="GWout (m3 d-1)",type='l',lwd=1)
plot(plotDates,outHydro[,2],xlab="time",ylab="Volume (m^3)",type='l',lwd=1)
plot(plotDates,outHydro[,4],xlab="time",ylab="stage (m)",type='l',lwd=1)
meanBudget=colMeans(hydroSumm[,-1],na.rm=TRUE)
barplot(meanBudget*c(1,1,1,-1,-1,-1),names.arg=names(meanBudget))

avgDailyIn=sum(meanBudget[1:3])
avgDailyOut=sum(meanBudget[4:6])
# estimate of residence time
print(mean(out[,2],na.rm=TRUE)/(mean(c(avgDailyIn,avgDailyOut))))


outConc=out
outConc[,3]=outConc[,3]/outConc[,2]*12		#mg L-1
outConc[,4]=outConc[,4]/outConc[,2]*12		#mg L-1

dev.new()
par(mfrow=c(2,2))
plot(plotDates,out[,2],type='l',xlab='time',ylab='Volume (m3)')
plot(plotDates,outConc[,3],type='l',xlab='time',ylab='DIC (mg C L-1)')
plot(plotDates,outConc[,4],type='l',xlab='time',ylab='DOC (mg C L-1)')

apply(outConc,2,range)
