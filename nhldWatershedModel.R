#### run equilibrium seepageVdrainage model for each lake in Vilas county
rm(list=ls())

setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/NHLD_lakeCmodel")

########## load utility functions
source("nhldWatershedModelSupporting.R")
require(deSolve)
require(LakeMetabolizer)

# load watershed table
#d=read.table("vilasSheds_5-18-14.txt",sep="\t",header=TRUE)

#W2L=d$shedArea_m2/d$lakeArea_m2
#toFix=which(W2L>100 & d$lakeArea_m2<10000)
#d$shedArea_m2[toFix]=d$lakeArea_m2[toFix]*1.5

#d=d[!(log10(RTs)<2 & log10(W2L)<0.5),]
#d=d[RTs>1,]

#W2L=d$shedArea_m2/d$lakeArea_m2

#RTs=d$lakeVol_m3/(dailyP/1000*d$shedAreaNoLake_m2-d$shedAreaNoLake_m2*shedETs/1000+dailyP/1000*d$lakeArea_m2)

### cut the oddball RTs????
#plot(log10(d$shedArea_m2/d$lakeArea_m2),RTs)


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

curFluxDOY=as.numeric(strftime(strptime(paste(curFlux[,4],curFlux[,5],curFlux[,3],sep="-"),format="%m-%d-%Y"),format="%j"))
curForceDOY=as.numeric(strftime(strptime(paste(curForce[,4],curForce[,5],curForce[,3],sep="-"),format="%m-%d-%Y"),format="%j"))

curFluxDOY=curFluxDOY-min(curFluxDOY)+1
curForceDOY=curForceDOY-min(curForceDOY)+1


#### repeate same forcing data 4 times just to see what happens
curForce=rbind(curForce,curForce)
curFlux=rbind(curFlux,curFlux)
curForce=rbind(curForce,curForce)
curFlux=rbind(curFlux,curFlux)

curFluxDOY=1:nrow(curFlux)
curForceDOY=rep(curFluxDOY,each=24)


startDOY=min(curFluxDOY)



setwd("/Volumes/JonesExternal/External/activeStuff/NHLD_Cmodel/NHLD_lakeCmodel")
###### forcing functions --> when running longer can just use January 1, 1950 as an epoch
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



# calculate concentrations of DIC, DOC, TP based on landcover
#setwd("/Files/NDongoing/ScalingLakeProcessesVilasCty/vilasWatersheds_post5-17-14/waterConstituents")
ntlGW=read.csv("NTLlter__groundwater_chemistry.csv",header=TRUE,fill=TRUE)
ntlGW=ntlGW[,c(1:4,9:10,14)]
ntlGW=ntlGW[!is.na(ntlGW[,5]),]
ntlGW=ntlGW[ntlGW$doc<10,]
#boxplot(ntlGW$dic~ntlGW$wellid,ylab="DIC")
#quartz()
#boxplot(ntlGW$doc~ntlGW$wellid,ylab="DOC")
#quartz()
#boxplot(ntlGW$totp~ntlGW$wellid,ylab="TP")

lottig2011DOC=read.table("lottig2012Fig.txt",header=FALSE,sep=",")

#plot(lottig2011DOC[,1],lottig2011DOC[,2],log='y',ylim=c(2,100))
#summary(lm(log(lottig2011DOC[,2])~lottig2011DOC[,1]))

streamDOCs=exp(1.3961+3.245*0.15)/12*1000/1000 #exp(1.3961+3.245*d$wetland)/12*1000/1000	# from lottig 2012; mol m-3
streamPOCs=3/12*1000/1000		# ~3 mg L-1; buffam 2011; mol m-3
streamDICs=10/12*1000/1000		#10 mg L-1; lottig 2011; mol m-3
streamPs=0.04/31*1000/1000		#.025 mg L-1 TDP & 0.04 mg L-1 TP; lottig 2011; mol m-3
gwDOCs=median(ntlGW$doc,na.rm=TRUE)/12*1000/1000	# mol m-3
gwDICs=median(ntlGW$dic,na.rm=TRUE)/12*1000/1000	# mol m-3
gwPs=median(ntlGW$totp,na.rm=TRUE)/31*1000/1e6		# mol m-3

kH=29.41	# Henry's Law constant for CO2 in water [L atm mol-1]

precipDOC=1.1/12*1000/1000		#mol m-3; from Likens et al. 1983 @ Hubbard Brook
precipDIC=400/1e6*1/kH*1000	# mol C m-3; assumed in equilibrium with atmosphere
precipP=0.01/31*1000/1000	# mol m-3; Murphy & DOskey 1976 JGLR


# iterate through lakes  --> perhaps in future run all lakes at once????
#for i in 1:nrow(d){
	# current lake and shed parameters
	curLakeArea=	80506.8				#d$lakeArea_m2[i]	#m2
	curLakeVol=322027.2					#d$lakeVol_m3[i]		#m3
	curMeanDepth=4				#curLakeVol/curLakeArea	#m
	curShedArea=	241959.1				#d$shedAreaNoLake_m2[i]	#m2
	cur_streamDOC=2.1				#streamDOCs[i]	# mol m-3
	gwQ=-179961.5*0.0283168				# m3 d-1 from GFLOW
	##******** this corresponds to 62 mm d-1; we think this is more like 1 mm d-1
		
	#### for now just divide GW by 60 to make some progress...
	gwQ=gwQ/60
	
	DR=((curLakeArea*10^-6)^0.5)/curMeanDepth	# dynamic ratio; unitless; from Hakanson & Bryhn 2008. Water Air Soil Pollution, 187: 119-147. impacts turbulence and sedimentation
	distShore=sqrt(curLakeArea/pi)/2		#mean distance from shore; m; from GIS

	if(DR<0.26){Ydr=DR/0.26}else{Ydr=0.26/DR} #adjustment to sedimentation rate based on turbulence, [unitless]; Hakanson & Bryhn 2008. Water Air Soil Pollution, 187: 119-147
	
##### WHAT TO DO WITH LITTLE LAKES WITH HUGE WATERSHEDS
	# are these "beaver ponds" with super low residence time?
	# are these delineation errors and they should be seepage lakes with WA2LA of ~1? 

# run model with epi driven by DOC-based lit. models
	# run continuously with mean open water season numbers and sum across that range for annual; assumes in winter lake is "shut off"; go with April 15 (day 89) to Dec 1 (335); 248 days

	# need a hypovolume initially to start simulation
	hypoVolume=0

	# combine parameters into vector
	params=c(DR=DR,curLakeArea=curLakeArea,curLakeVol=curLakeVol,Ydr=Ydr,curMeanDepth=curMeanDepth,curShedArea=curShedArea)

	# initial values
	initialX=c(DOC=8/12*curLakeVol,DIC=0.3/12*curLakeVol,tPOC=0.1/12*curLakeVol,Phyto=1/12*curLakeVol,SRP=1/31/1000*curLakeVol,Buried=0)

	nsteps=nrow(curFlux)
	
	times=curFluxDOY
	times=times[-length(times)]	
	
	out<-ode(y=initialX,times=times,func=timeStep,parms=params,curForce=curForce,curForceDOY=curForceDOY)
	#out<-ode(y=initialX,times=times,func=timeStep,parms=params,curForce=curForce,curForceDOY=curForceDOY,method='iteration')

	dev.new()
	par(mfrow=c(2,3))
	plot(out[,1],out[,2]/out[,9]*12,xlab="time",ylab="DOC (mg C L-1)",type='l')
	plot(out[,1],out[,3]/out[,9]*12,xlab="time",ylab="DIC (mg C L-1)",type='l')
	plot(out[,1],out[,4]/out[,9]*12,xlab="time",ylab="tPOC (mg C L-1)",type='l')
	plot(out[,1],out[,5]/out[,9]*12,xlab="time",ylab="Phyto C (mg C L-1)",type='l')
	plot(out[,1],out[,6]/out[,9]*1000*31,xlab="time",ylab="SRP (ug P L-1)",type='l')
	
	apply(out,2,range)

# store equilibrium concentrations and estimate daily/annual fluxes

#}
# compare to observations (Cole and Tranvik, classification based estimates, Hanson 2004, Hanson 2003, Hanson 2007... others?)