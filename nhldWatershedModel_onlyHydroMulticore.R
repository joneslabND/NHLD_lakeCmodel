##### just hydrologic model for now...

# clear workspace
rm(list=ls())

########## load utility functions and packages
source("nhldWatershedModel_onlyHydroSupporting.R")
require(deSolve)
require(LakeMetabolizer)
require(multicore)

######## load forcing and flux data
setwd("/afs/crc.nd.edu/user/s/sjones20/Private/NHLD_lakeCmodel")

# Lat Long Year Month Day Runoff_mmDay Precip_mmDay EvapOpenWater_mmDay
flux=read.table("FLUX_46.21875_-89.46875",header=FALSE)
# Lat Long Year Month Day Hour NetLW_Wm2 NetSW_Wm2 LW_Wm2 SW_Wm2 AirTemp_dC atmPress_kPa windspeed_mS@10m relHumid_fraction
force=read.table("FORCE_46.21875_-89.46875",header=FALSE)

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
for(j in 1:length(years)){
	curYear=curFlux[curFlux[,3]==years[j],]
	curDOY=as.numeric(strftime(strptime(paste(curYear[,3],curYear[,4],curYear[,5],sep="-"),format="%Y-%m-%d"),format="%j"))
	iceON[curFlux[,3]==years[j]]=((curDOY<iceOffDOY[names(iceOffDOY)==years[j]]) | (curDOY>iceOnDOY[names(iceOnDOY)==years[j]]))*1
}

####### load lake/watershed information
UNDERCsheds=read.table("UNDERCsheds_4-24-15.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
GFLOWoutput=read.table("NHLD_GFLOW_UNDERC_20150520_Simplify010_RegK1_LocK0pt01_OutputLong.txt",header=TRUE,stringsAsFactors=FALSE)
colnames(GFLOWoutput)[1]="Permanent_"


lakeSim<-function(curLakeID){	
	
	# current lake and shed parameters
	A0=UNDERCsheds$NHLD_lakes[UNDERCsheds$Permanent_==curLakeID]		#m2
	V0=10^(-0.0589+1.12963*log10(A0))		#m3
	Perim0=UNDERCsheds$Perimeter[UNDERCsheds$Permanent_==curLakeID]			#******* how to make this dynamic???
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
	zbar1=zbar0*1.5
	zmax1=zmax0*1.5
	
	uScale=(zmax1-zmax0)/zmax1
						
	A1=A0/(1-(p*uScale^2+(1-p)*uScale))
	V1=V0/(1-((6*uScale-3*(1-p)*uScale^2-2*p*uScale^3)/(3+p)))
					
	curShedArea=UNDERCsheds$Area_m2[UNDERCsheds$Permanent_==curLakeID]	#m2

	u0=round(uniroot(f=findU,lower=0,upper=1,p=p,Vmax=V1,Vu=(V1-V0))$root,4)
	
	curGFLOW=GFLOWoutput[GFLOWoutput$Permanent_==curLakeID,]
	curGFLOW=curGFLOW[!is.na(curGFLOW$Permanent_),]
	GFLOWpropPerim=curGFLOW$Linesink_Length_Output/sum(curGFLOW$Linesink_Length_Output)
	GFLOWin=(curGFLOW$Linesink_PerLengthDischarge_Output>0)*1
	
	gwIn0=sum(GFLOWin*curGFLOW$Linesink_PerLengthDischarge_Output*GFLOWpropPerim*Perim0)*0.0283168	#m3 d-1
	gwOut0=-sum((1-GFLOWin)*curGFLOW$Linesink_PerLengthDischarge_Output*GFLOWpropPerim*Perim0)*0.0283168	#m3 d-1
		
	stage0=zmax1-u0*zmax1
	alpha=0.99
	stageOut=alpha*stage0

	params=c(Vmax=V1,Zmax=zmax1,Amax=A1,curShedArea=curShedArea,stageOut=stageOut,Perim0=Perim0,gwIn0=gwIn0,gwOut0=gwOut0,p=p,DL=DL)

	initialX=c(V=V0)
	times=curFluxDOY
	
	out<-ode(y=initialX,times=times,func=timeStep,parms=params)

	out=cbind(out,NA,NA,NA)
	for(j in 1:nrow(out)){
		if(!is.na(out[j,2])){
			uj=round(uniroot(f=findU,lower=0,upper=1,p=p,Vmax=V1,Vu=out[j,2])$root,4)
			out[j,3]=A1*(p*uj^2+(1-p)*uj)
			out[j,4]=zmax1*uj
			out[j,5]=2*pi*sqrt(out[j,3]/pi)*DL
		}
	}
	colnames(out)[3:5]=c("A","stage","perim")

	C=(2/3)^1.5*9.806^0.5	# m s-2
	L=0.1	# m
	H=out[,4]-stageOut
	QoutSim=ifelse(H>0,(C*L*H^1.5)*(60*60*24),0)
	
	hydroSumm=cbind(out[,1],dailyRunoff(out[,1])*(curShedArea-out[,3])/1000,dailyPrecip(out[,1])*out[,3]/1000,gwIn0*out[,5]/Perim0,QoutSim,dailyEvap(out[,1])*out[,3]/1000*(1-iceON[1:nrow(out)]),gwOut0*out[,5]/Perim0)
	colnames(hydroSumm)=c("time","Qin","precip","GWin","Qout","evap","GWout")
	meanBudget=colMeans(hydroSumm[,-1],na.rm=TRUE)

	avgDailyIn=sum(meanBudget[1:3])
	avgDailyOut=sum(meanBudget[4:6])
	# estimate of residence time
	print(c(curLakeID,mean(out[,2],na.rm=TRUE)/(mean(c(avgDailyIn,avgDailyOut)))))
}

lakes=sort(unique(GFLOWoutput$Permanent_))
lakes=as.list(lakes[lakes%in%UNDERCsheds$Permanent_])

mclapply(lakes,lakeSim,mc.cores=12)