#nhldWatershedModelSupporting.R supporting functions

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
		Emit=X[7]
		Sed_tPOC=X[8]
		Sed_phyto=X[9]

		stage=((3*V)/(r2h^2*pi))^(1/3)
		A=pi*(r2h*stage)^2
		
		perim=2*pi*sqrt(A/pi)*DL			#### maybe just treat this as a circle???

		if(DOC>0){
			zmix=10^(-0.518*log10(DOC/V*12)+1.006)		#fuentetaja et al. 1999
		}else{
			zmix=stage
		}
		#print(c(t,V,stage,DOC,zmix,DIC,phyto,P))

		if(zmix>stage){
			zmix=stage
		}
		
		r1=r2h*stage
		r2=r2h*(stage-zmix)
		Vepi=(1/3)*pi*(r1^2+r1*r2+r2^2)*zmix

		gwIn=gwIn0*perim/Perim0	#m3 d-1
		gwOut=gwOut0*perim/Perim0	#m3 d-1
				
		streamQ=dailyRunoff(t)*(curShedArea-A)/1000		#m3 day-1	****** problem because of interpolating?
	
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
		
		phytoC2P=95#106		#M:M; from Patrick
	
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

		#dDOC.dt=streamQ*streamDOC+directPrecip*precipDOC+gwIn*gwDOC+GPP*Vepi*GPPexude-STout*DOC/V-gwOut*DOC/V-photoOx*A-floc*DOC-DOCrespired*DOC	#[mol C]
		dDOC.dt=streamQ*streamDOC+directPrecip*precipDOC+gwIn*gwDOC-STout*DOC/V-gwOut*DOC/V-photoOx*A-floc*DOC-DOCrespired*DOC	#[mol C]

		dtPOC.dt=streamQ*streamPOC+deposition+floc*DOC-STout*tPOC/V-settling*tPOC	#[mol C]
		
		#dphyto.dt=GPP*Vepi-GPP*Vepi*Rauto-settling*phyto-STout*phyto/V-GPP*Vepi*GPPexude
		dphyto.dt=GPP*Vepi-GPP*Vepi*Rauto-settling*phyto-STout*phyto/V
		
		dP.dt=streamQ*streamP+directPrecip*precipP+gwIn*gwP-STout*P/V-gwOut*P/V-GPP*Vepi*(1-Rauto)/phytoC2P
		
		dEmit.dt=-kCur/zmix*(atmEquilCO2-DIC/V)*A
		
		dSed_tPOC.dt=settling*tPOC
		
		dSed_phyto.dt=settling*phyto
		
		list(c(dV.dt,dDIC.dt,dDOC.dt,dtPOC.dt,dphyto.dt,dP.dt,dEmit.dt,dSed_tPOC.dt,dSed_phyto.dt))
	})

}