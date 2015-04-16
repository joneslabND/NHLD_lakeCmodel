#### run equilibrium seepageVdrainage model for each lake in Vilas county
rm(list=ls())

setwd("/Files/NDongoing/ScalingLakeProcessesVilasCty/vilasWatersheds_post5-17-14")

########## load utility functions
source("vilasWatershedModelSupporting_5-18-14.R")
require(deSolve)

# load watershed table
d=read.table("vilasSheds_5-18-14.txt",sep="\t",header=TRUE)

W2L=d$shedArea_m2/d$lakeArea_m2
toFix=which(W2L>100 & d$lakeArea_m2<10000)
d$shedArea_m2[toFix]=d$lakeArea_m2[toFix]*1.5

#d=d[!(log10(RTs)<2 & log10(W2L)<0.5),]
#d=d[RTs>1,]

W2L=d$shedArea_m2/d$lakeArea_m2


# load met data
setwd("/Files/NDongoing/ScalingLakeProcessesVilasCty/vilasWatersheds_post5-17-14/metData")
metD=read.csv("NTLwoodruffairport_daily_metdata.csv",header=TRUE)
# remove flags and unnecessary data
metD=metD[,c(1:7,9,11,18,20,24,33:35)]
metH=read.csv("NTLwoodruffairport_hourly_metdata.csv",header=TRUE)
metH=metH[,c(1:8,11,16:19)]

# run continuously with mean open water season numbers and sum across that range for annual; assumes in winter lake is "shut off"; go with April 15 (day 89) to Dec 1 (335); 248 days
seasonLength=248	#days

# annual precip
annP=tapply(metD[,5],metD[,1],sum,na.rm=TRUE)
summary(annP)
annP[names(annP)==2000]

openwaterD=metD[(metD$daynum>88 & metD$daynum<335),]
summary(tapply(openwaterD[,5],openwaterD[,1],sum,na.rm=TRUE)/seasonLength)

#dailyP=median(tapply(openwaterD[,5],openwaterD[,1],sum,na.rm=TRUE)[-23])/seasonLength	# mm day-1
# 2011 has issues
dailyP=800/seasonLength


# evaporation
allEvaps=evapCalc(openwaterD$avg_air_temp,openwaterD$daynum,lat=46.05)
dailyEvap=mean(allEvaps,na.rm=TRUE)

# PAR
openwaterH=metH[(metH$daynum>88 & metH$daynum<335),]
hourlyPAR=tapply(openwaterH$avg_par,openwaterH$hour,mean,na.rm=TRUE)

#plot(as.numeric(names(hourlyPAR))/100,hourlyPAR)

# peak wind speed for tPOC deposition
maxWind=median(openwaterD$peak_wind_speed_5sec,na.rm=TRUE)	#NOTE this ends up being above range presented in Preston 2008???

# calculate ET from landcover --> from year 2000...
setwd("/Files/NDongoing/ScalingLakeProcessesVilasCty/vilasWatersheds_post5-17-14")
et=read.table("mckayETestimates.txt",header=TRUE,sep="\t",skip=1)
etUp=median(c(et[,2],et[,3],et[,4],et[,5],et[,6],et[,7]))+median(et[,10])	#mm day-1
etLow=median(c(et[,8],et[,9]))+median(et[,11])	#mm day-1

shedETs=etUp*d$forest+etLow*d$wetland
range(shedETs)	#mm day-1

# calculate concentrations of DIC, DOC, TP based on landcover
setwd("/Files/NDongoing/ScalingLakeProcessesVilasCty/vilasWatersheds_post5-17-14/waterConstituents")
ntlGW=read.csv("NTLlter__groundwater_chemistry.csv",header=TRUE,fill=TRUE)
ntlGW=ntlGW[,c(1:4,9:10,14)]
ntlGW=ntlGW[ntlGW$doc<10,]
#boxplot(ntlGW$dic~ntlGW$wellid,ylab="DIC")
#quartz()
#boxplot(ntlGW$doc~ntlGW$wellid,ylab="DOC")
#quartz()
#boxplot(ntlGW$totp~ntlGW$wellid,ylab="TP")

lottig2011DOC=read.table("lottig2012Fig.txt",header=FALSE,sep=",")

#plot(lottig2011DOC[,1],lottig2011DOC[,2],log='y',ylim=c(2,100))
#summary(lm(log(lottig2011DOC[,2])~lottig2011DOC[,1]))

streamDOCs=exp(1.3961+3.245*d$wetland)/12*1000/1000	# from lottig 2012; mol m-3
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

RTs=d$lakeVol_m3/(dailyP/1000*d$shedAreaNoLake_m2-d$shedAreaNoLake_m2*shedETs/1000+dailyP/1000*d$lakeArea_m2)

### cut the oddball RTs????
#plot(log10(d$shedArea_m2/d$lakeArea_m2),RTs)

##############################################
#### setting lake and simulation features ####
##############################################
nsteps=1000	 	# number of days in simulation







#i=1


for i in 1:nrow(d){
	# current lake and shed parameters
	curLakeArea=d$lakeArea_m2[i]	#m2
	curLakeVol=d$lakeVol_m3[i]		#m3
	curMeanDepth=curLakeVol/curLakeArea	#m
	curShedArea=d$shedAreaNoLake_m2[i]	#m2
	curET=shedETs[i]	#mm day-1
	cur_streamDOC=streamDOCs[i]	# mol m-3
	
	DR=((curLakeArea*10^-6)^0.5)/curMeanDepth	# dynamic ratio; unitless; from Hakanson & Bryhn 2008. Water Air Soil Pollution, 187: 119-147. impacts turbulence and sedimentation
	distShore=sqrt(curLakeArea/pi)/2		#mean distance from shore; m; from GIS

	if(DR<0.26){Ydr=DR/0.26}else{Ydr=0.26/DR} #adjustment to sedimentation rate based on turbulence, [unitless]; Hakanson & Bryhn 2008. Water Air Soil Pollution, 187: 119-147
	
	
	#hydrologic inputs
	shedHydraulicLoad=dailyP/1000*curShedArea-curShedArea*curET/1000	#m3 d-1
	
	# partitioning between GW and surface -- using pseudo horton infiltration model
		#additional GW inputs???
	GWcond=1.8		#mm day-1		KNOB
	
	if(GWcond>=(dailyP-curET)){
		fractionGW=1
		fractionSW=0
	}else{
		fractionGW=GWcond/(dailyP-curET)
		fractionSW=1-fractionGW
	}
	
	streamQ=shedHydraulicLoad*fractionSW	#m3 d-1
	gwQ=shedHydraulicLoad*fractionGW		#m3 d-1
	directPrecip=dailyP/1000*curLakeArea	#m3 d-1
	
	
	#hydrologic outputs
	curEvap=dailyEvap*curLakeArea/1000	#m3 d-1
	
	hydroOutput=shedHydraulicLoad+directPrecip-curEvap
	
	# maybe this split doesn't matter for now (would only modify POC export via stream that thus wouldn't settle)
	STout=hydroOutput*fractionSW		#m3 d-1
	GWout=hydroOutput*fractionGW		#m3 d-1

##### WHAT TO DO WITH LITTLE LAKES WITH HUGE WATERSHEDS
	# are these "beaver ponds" with super low residence time?
	# are these delineation errors and they should be seepage lakes with WA2LA of ~1? 

	# daily loads of DIC, DOC, TP
	DICin=directPrecip*precipDIC+streamQ*streamDICs+gwQ*gwDICs
	DOCin=directPrecip*precipDOC+streamQ*cur_streamDOC+gwQ*gwDOCs
	Pin=directPrecip*precipP+streamQ*streamPs+gwQ*gwPs
	
	tPOCdepTotal=integrate(depWithDist,0.01,sqrt(curLakeArea/pi),precip=dailyP,maxWind=maxWind)$value/sqrt(curLakeArea/pi)*curLakeArea #tPOCdep summed across lake surface	#mol C day-1
	POCin=tPOCdepTotal+streamQ*streamPOCs

	
# run model with epi driven by DOC-based lit. models
	# run continuously with mean open water season numbers and sum across that range for annual; assumes in winter lake is "shut off"; go with April 15 (day 89) to Dec 1 (335); 248 days

	# need a hypovolume initially to start simulation
	hypoVolume=0

	# combine parameters into vector
	params=c(curLakeVol=curLakeVol,DR=DR,DOCin=DOCin,DICin=DICin,Pin=Pin,POCin=POCin,LA=curLakeArea,STout=STout,GWout=GWout,tPOCdepTotal=tPOCdepTotal,Ydr=Ydr,curMeanDepth=curMeanDepth)

	# initial values
	initialX=c(DOC=8/12*curLakeVol,DIC=0.3/12*curLakeVol,tPOC=0.1/12*curLakeVol,Phyto=1/12*curLakeVol,SRP=10/31/1000*curLakeVol,Buried=0)

	times=1:nsteps	
	
	out<-ode(y=initialX,times=times,func=timeStep,parms=params)
	print(out[nrow(out),2]/out[nrow(out),9]*12)

# store equilibrium concentrations and estimate daily/annual fluxes

}
# compare to observations (Cole and Tranvik, classification based estimates, Hanson 2004, Hanson 2003, Hanson 2007... others?)