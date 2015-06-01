#nhldWatershedModelSupporting.R supporting functions

findU<-function(u,p,Vmax,Vu){
	((6*u-3*(1-p)*u^2-2*p*u^3)/(3+p))*Vmax-Vu
} 

#light attenuation
lightAtten<-function(z,I0,kD){
	Iz=I0*exp(-kD*z)
	return(Iz)
}
lightAttenTS<-function(I0,kD,zmix){
	avgI=numeric(length(I0))
	for(i in 1:length(avgI)){
		avgI[i]=integrate(lightAtten,0,zmix,I0=I0[i],kD=kD)$value/zmix
	}
	return(avgI)
}

# dailyGPP based on hourly light, etc.
dailyGPP<-function(day,curForce,curForceDOY,chlCur,SRP,DOC,V,kD,zmix){
	hourlyPAR=sw.to.par.base(curForce[curForceDOY==day,10])		# umol m2 sec; SW to PAR based on Read...

	#### trying different GPP formulation
	PPmax=20#3.1#2.2		# hr-1
	kSRP=0.3/1000#0.25/1000/31#0.5/1000	# mol P m-3; Halmann and Stiller 1974 L&O 19(5): 774-783 (Lake Kinnerrett)
	Ik=180 		# light limitation benchmark for GPP, [umol cm-2 s-1]; Vadeboncoeur et al. 2008
			
	avgI=lightAttenTS(hourlyPAR,kD,zmix)		# hourly average light climate in mixed layer [umol cm-2 s-1]
	GPP=sum(chlCur*PPmax*tanh(avgI/Ik)*((SRP/V)/((SRP/V)+kSRP)))/(12*1000) # mol C m-3 day-1
	
	return(GPP)
}

# atmospheric deposition
depWithDist<-function(distShore,precip,maxWind){
	lPOC=10^(0.43+0.0034*precip+0.11*maxWind-1.05*distShore/(48+distShore))/(12*1000)	#mol C m-2 day-1
	sPOC=10^(0.0082+0.0068*precip+0.12*maxWind-0.6*distShore/(49+distShore))/(12*1000) #mol C m-2 day-1
	
	tPOCdep=(lPOC+sPOC) #tPOCdep summed across lake surface	#mol C m-2 day-1
	return(tPOCdep)
}

# daily tPOC deposition based on Preston et al.
dailyTPOCdep<-function(day,curLakeArea,curForce,curForceDOY){
	maxWind=max(curForce[curForceDOY==day,13])*(1/10)^0.15		# convert 10m wind to 1m for preston model tPOC deposition

	tPOCdepTotal=integrate(depWithDist,0.01,sqrt(curLakeArea/pi),precip=dailyPrecip(day),maxWind=maxWind)$value/sqrt(curLakeArea/pi)*curLakeArea #tPOCdep summed across lake surface	#mol C day-1
	return(tPOCdepTotal)
}

timeStep<-function(t,X,params){
	with(as.list(params),{
		
		#########################
		#### State variables ####
		#########################
		V=X[1]
		DIC=X[2]
		DOC=X[3]
		tPOC=X[4]
		phyto=X[5]
		P=X[6]

		u=round(uniroot(f=findU,lower=0,upper=1,p=p,Vmax=Vmax,Vu=(Vmax-V))$root,4)
		
		stage=Zmax-u*Zmax
		A=Amax*(1-(p*u^2+(1-p)*u))
		
		perim=2*pi*sqrt(A/pi)*DL

		zmix=10^(-0.518*log10(DOC/V*12)+1.006)		#fuentetaja et al. 1999
		if(zmix>Zmax){
			zmix=Zmax
		}
		
		umix=zmix/stage
		Vu=(6*umix-3*(1-p)*umix^2-2*p*umix^3)/(3+p)
		Vepi=V*Vu

		gwIn=gwIn0*perim/Perim0	#m3 d-1
		gwOut=gwOut0*perim/Perim0	#m3 d-1
				
		streamQ=dailyRunoff(t)*(curShedArea-A)/1000			#m3 day-1	****** problem because of interpolating?
		directPrecip=dailyPrecip(t)*A/1000		#m3 day-1	****** problem because of interpolating?
		
		U10=dailyWind(t)
		
		if(iceON[floor(t)]==1){
			curEvap=0	#m3
			kCur=0		#m day-1
			DOCrespired=0.0005 # fraction of DOC that is respired, [day-1]; low for temperature effect
		}else{
			curEvap=dailyEvap(t)*A/1000			#m3 day-1	
			kCur=(2.51+1.48*U10+0.39*U10*log10(A))*24/100	#m day-1 from Vachon & Prairie 2013
			DOCrespired=0.002 # fraction of DOC that is respired, [day-1]; warmer than under ice, but still below 0.005 because of spring/fall
		}
		
		# assuming ogee crest spillway
		# Q=C*L*H^(3/2), where:
		# Q = discharge [m3 s-1]
		# C = coefficient; (2/3)^1.5*g^0.5 [m^(1/2) s-1]; g = gravitational constant; 9.806 [m s-2]
		# L = spillway length [m]; call this 2-5 m?
		# H = head (difference between spillway height and reservoir stage) [m]
		C=(2/3)^1.5*9.806^0.5	
		L=0.1#0.25	# m
		if((stage-stageOut)>0){
			H=stage-stageOut
			Q=C*L*H^1.5	#m3 s-1
			STout=Q*60*60*24		#m3 day-1
		}else{
			STout=0		#m3 day-1
		}

		# C biogeochemistry
		photoOx=0#	******* THIS IS TURNED OFF RIGHT NOW!!!!!! 44/1000/12 	# photooxidation rate constant, [mol c m-2 day-1]; Graneli et al. 1996, L&O, 41: 698-706
		floc=0.#0.005		# fraction of DOC that floculates, [day-1]; von Wachenfeldt & Tranvik 2008 via Jones et al. 2012
		#DOCrespired=0.001#0.005 # fraction of DOC that is respired, [day-1]; Houser et al. 2001 via Hanson et al. 2004

		# epilimnion vs. hypolimion 
		#zmix=10^(-0.518*log10(DOC/(V*1000)*12*1000)+1.006)	#****** NEED TO MAKE SURE THE DOC UNITS ARE RIGHT HERE...	#fuentetaja et al. 1999
		phyto_C2Chl=50		# phytoplankton carbon to chlorophyll ratio, [gC gChl-1]; sort of made up/average of observations in paper (e.g. 
		chlCur=phyto/Vepi*12*1000/phyto_C2Chl # current chlorophyll concentrations [mg m-3]
		kD=0.0213+0.0177*chlCur+0.0514*(DOC/V*12) # current light attenuation coefficient [m-1]  #OR kD=-0.217+0.0537*chloro+0.186*DOC --> in M2M lightProfile Calcs r script
		
		atmCO2=400 # [ppm]
		atmEquilCO2=1*atmCO2/1e6/kH*1000	# concentration of CO2 in water at equilibrium with atmosphere, [mol C m-3]
		
		phytoC2P=106#70		#M:M; from Patrick
	
		deposition=dailyTPOCdep(floor(t),A,curForce,curForceDOY)
		
		d=5		#particle diameter [um]
		settling=0.0188*(d/2)^2/zmix
		
		GPP=dailyGPP(day=floor(t),curForce=curForce,curForceDOY=curForceDOY,chlCur=chlCur,SRP=P,DOC=DOC,V=V,kD=kD,zmix=zmix)
		
		GPPexude=0.03	# Hanson et al. 2004
		Rauto=0.8		# Hanson et al. 2004
				
		################################
		#### Differential equations ####
		################################

		dV.dt=(streamQ+directPrecip+gwIn)-(STout+curEvap+gwOut)	#[m3]
				
		dDIC.dt=streamQ*streamDIC+directPrecip*precipDIC+gwIn*gwDIC-STout*DIC/V-gwOut*DIC/V+photoOx*A+DOCrespired*DOC+kCur/zmix*(atmEquilCO2-DIC/V)*A-GPP*Vepi+GPP*Rauto*Vepi	#[mol C]

		dDOC.dt=streamQ*streamDOC+directPrecip*precipDOC+gwIn*gwDOC+GPP*Vepi*GPPexude-STout*DOC/V-gwOut*DOC/V-photoOx*A-floc*DOC-DOCrespired*DOC	#[mol C]

		dtPOC.dt=streamQ*streamPOC+deposition+floc*DOC-STout*tPOC/V-settling*tPOC	#[mol C]
		
		#dphyto.dt=GPP*zmix*A-Rauto*phyto-settling*phyto-STout*phyto/V-GPP*zmix*A*GPPexude
		dphyto.dt=GPP*Vepi-GPP*Vepi*Rauto-settling*phyto-STout*phyto/V-GPP*Vepi*GPPexude

		
		#dP.dt=streamQ*streamP+directPrecip*precipP+gwIn*gwP-STout*P/V-gwOut*P/V-GPP*zmix*A/phytoC2P+phyto*Rauto/phytoC2P*recyclingEff
		dP.dt=streamQ*streamP+directPrecip*precipP+gwIn*gwP-STout*P/V-gwOut*P/V-GPP*Vepi*(1-Rauto)/phytoC2P
		
		list(c(dV.dt,dDIC.dt,dDOC.dt,dtPOC.dt,dphyto.dt,dP.dt))
	})

}










###### relic code from original equilibrium model


# lake process model
timeStepOLD<-function(t,X,params,curForce,curForceDOY){
	with(as.list(params),{
		
		print(t)
		
		#########################
		#### State variables ####
		#########################
		DOC=X[1]
		DIC=X[2]
		tPOC=X[3]
		Phyto=X[4]
		SRP=X[5]
		Buried=X[6]
		
		streamQ=dailyRunoff(t)*curShedArea/1000			#m3 day-1	****** problem because of interpolating?
		directPrecip=dailyPrecip(t)*curLakeArea/1000		#m3 day-1	****** problem because of interpolating?
	
		if(gwQ>0){
			gwIn=gwQ
			gwOut=0
		}else{
			gwIn=0
			gwOut=abs(gwQ)
		}
		
		curEvap=dailyEvap(t)*curLakeArea/1000			#m3 day-1	****** problem because of interpolating?
		
		STout=(directPrecip+streamQ+gwIn)-(curEvap+gwOut)	#m3 day-1;  currently holds volume constant...

		# daily loads of DIC, DOC, TP
		DICin=directPrecip*precipDIC+streamQ*streamDICs+gwQ*gwDICs
		DOCin=directPrecip*precipDOC+streamQ*cur_streamDOC+gwQ*gwDOCs
		Pin=directPrecip*precipP+streamQ*streamPs+gwQ*gwPs
	
		#maxWind=dailyMaxWind(floor(t))
		
		#hourlyPAR=dailySun(floor(t))
	
		#tPOCdepTotal=integrate(depWithDist,0.01,sqrt(curLakeArea/pi),precip=dailyPrecip(t),maxWind=maxWind)$value/sqrt(curLakeArea/pi)*curLakeArea #tPOCdep summed across lake surface	#mol C day-1
		
		tPOCdepTotal=dailyTPOCdep(floor(t),curLakeArea,curForce,curForceDOY)
		POCin=tPOCdepTotal+streamQ*streamPOCs


# In lake processes
#Ik=180 		# light limitation benchmark for GPP, [umol cm-2 s-1]; Vadeboncoeur et al. 2008
#photoOx=44/1000/12#0.0042 	# photooxidation rate constant, [mol c m-2 day-1]; Graneli et al. 1996, L&O, 41: 698-706
photoOx=0
exude=0.03		# fraction of GPP released as DOC, [day-1]; Biddanda & Benner 1997 via Hanson et al. 2004
leafLeach=0.12/14	# fraction of leaf load released as DOC, [day-1]; France et al. 1997; saw 6-18% loss over two weeks  do this on an annual leaf load * leach rate OR make this seasonal; this also has to be input to sediment/hypo OM
floc=0.005		# fraction of DOC that floculates, [day-1]; von Wachenfeldt & Tranvik 2008 via Jones et al. 2012
GPPrespired=0.8	# fraction of GPP that is respired quickly, [day-1]; Quay et al. 1996, Cole et al. 2002 via Hanson et al. 2004
DOCrespired=0.005 # fraction of DOC that is respired, [day-1]; Houser et al. 2001 via Hanson et al. 2004
phyto_C2P=106		# phytoplankton carbon to posphorus molar ratio, [molC molP-1]; redfield
#phyto_C2Chl=50		# phytoplankton carbon to chlorophyll ratio, [gC gChl-1]; sort of made up/average of observations in paper (e.g. Wang et al. 2008 Biogeoscience Discussions, 5: 3869-3903; Riemann et al. 1989. JPR, 11: 1037-1045)

refPOM=2.5 # reference POM concentration to adjust sedimentation rate, [mol C m-3]; Hakanson & Bryhn 2008. Water Air Soil Pollution, 187: 119-147
vDef=72/365 # default sedimentation rate, [m d-1]; Hakanson & Bryhn 2008. Water Air Soil Pollution, 187: 119-147

V=curLakeVol-hypoVolume

zmix=10^(-0.518*log10(DOC/(V*1000)*12*1000)+1.006)		#fuentetaja et al. 1999

#print(DOC/V*12)

if(zmix>curMeanDepth){
	zmix=curMeanDepth
}

hypoVolume<<-curLakeArea*(curMeanDepth-zmix)	# m3

#print(c(zmix,hypoVolume))

autoSEDrespired=0.8	# sediment respiration rate for autochthonous carbon [day-1]; mean from Sobek et al. 2009 L&O, 54: 2243-2254 --> determines TP recycling with C2P???
alloSEDrespired=0.3  # sediment respiration rate for allochthonous carbon [day-1]; mean from Sobek et al. 2009 L&O, 54: 2243-2254 --> determines TP recycling with C2P???

atmCO2=400 # [ppm]


#### trying different GPP formulation
PPmax=15#2.2		# day-1
kSRP=0.5/1000	# mol P m-3; Halmann and Stiller 1974 L&O 19(5): 774-783 (Lake Kinnerrett)

		################################
		#### Intermediate equations ####
		################################
		kH=29.41	# Henry's Law constant for CO2 in water [L atm mol-1]	# temperature sensitivity
		atmEquilCO2=1*atmCO2/1e6/kH*1000	# concentration of CO2 in water at equilibrium with atmosphere, [mol C m-3]
	
		#mass basis:
		Yconc=1+0.75*((tPOC/V+Phyto/V)/refPOM-1)	# adjustment to sedimentation rate for POM concentration, assuming tPOC and Phyto interact similarly when aggregating [unitless]
		sedRate=(vDef/zmix)*Ydr*Yconc # sedimentation rate [d-1]
		#chlCur=Phyto/V*12*1000/phyto_C2Chl # current chlorophyll concentrations [mg m-3]
		#kD=0.0213+0.0177*chlCur+0.0514*DOC/V # current light attenuation coefficient [m-1]  #OR kD=-0.217+0.0537*chloro+0.186*DOC --> in M2M lightProfile Calcs r script
				
		#avgI=lightAttenTS(hourlyPAR,kD,zmix)		# hourly average light climate in mixed layer [umol cm-2 s-1]
		#GPP=sum(chlCur*PPmax*tanh(avgI/Ik)*((SRP/V)/((SRP/V)+kSRP)))/(12*1000) # mol C m-3 day-1
		GPP=dailyGPP(floor(t),curForce,curForceDOY,Phyto,SRP,DOC,V,zmix)
		NPP=GPP*(1-GPPrespired-exude)			# mol C m-3 day-1
		#R=GPP*GPPrespired+DOC/V*DOCrespired		# mol C m-3 day-1
		R=GPP*GPPrespired+DOC/V*DOCrespired	# mol C m-3 day-1


		kCur=0.5 	# current piston velocity, using Jordan Read/GLEON model eventually
		# should be function of area...
		
		#print(V)

		################################
		#### Differential equations ####
		################################

		dDOC.dt=DOCin+GPP*exude*V-photoOx*curLakeArea-floc*DOC-DOCrespired*DOC-(STout+gwOut)/V*DOC	#[mol C]
		
		dDIC.dt=DICin+R*V+photoOx*curLakeArea-GPP*V-(STout+gwOut)/V*DIC+kCur*(atmEquilCO2-DIC/V)*curLakeArea	#[mol C]

		dtPOC.dt=POCin+tPOCdepTotal+floc*DOC-sedRate*tPOC-STout/V*tPOC	#[mol C]

		dPhyto.dt=NPP*V-sedRate*Phyto-STout/V*Phyto	#[mol C]

		# lacking any recycling in epi...
		dSRP.dt=Pin-NPP*V/phyto_C2P-(STout+gwOut)/V*SRP	# [mol P]

		dBuried.dt=sedRate*tPOC*(1-alloSEDrespired)+sedRate*Phyto*(1-autoSEDrespired)	#[mol C]
	
		
		list(c(dV.dt,dDOC.dt,dDIC.dt,dtPOC.dt,dPhyto.dt,dSRP.dt,dBuried.dt))
		#list(c(DOC+dDOC.dt,DIC+dDIC.dt,tPOC+dtPOC.dt,Phyto+dPhyto.dt,SRP+dSRP.dt,Buried+dBuried.dt),c(zmix=zmix,V=V))
	})

}
