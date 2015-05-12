# nhldWatershedModel_onlyHydro supporting functions

findU<-function(u,p,Vmax,Vu){
	((6*u-3*(1-p)*u^2-2*p*u^3)/(3+p))*Vmax-Vu
}


timeStep<-function(t,X,params){
	with(as.list(params),{
		#########################
		#### State variables ####
		#########################
		V=X[1]
		
		u=round(uniroot(f=findU,lower=0,upper=1,p=p,Vmax=Vmax,Vu=V)$root,4)
		
		stage=u*Zmax
		A=Amax*(p*u^2+(1-p)*u)
		
		perim=2*pi*sqrt(A/pi)*DL

		gwIn=gwIn0*perim/Perim0	#m3 d-1
		gwOut=gwOut0*perim/Perim0	#m3 d-1
				
		# because we are accounting for direct precip isn't this the snow on ice? 
		# I guess we could "store" that precip until spring...
		streamQ=dailyRunoff(t)*(curShedArea-A)/1000			#m3 day-1	****** problem because of interpolating?
		directPrecip=dailyPrecip(t)*A/1000		#m3 day-1	****** problem because of interpolating?
		
		if(iceON[floor(t)]==1){
			curEvap=0	#m3
		}else{
			curEvap=dailyEvap(t)*A/1000			#m3 day-1	
		}
		
		# assuming ogee crest spillway
		# Q=C*L*H^(3/2), where:
		# Q = discharge [m3 s-1]
		# C = coefficient; (2/3)^1.5*g^0.5 [m^(1/2) s-1]; g = gravitational constant; 9.806 [m s-2]
		# L = spillway length [m]; call this 2-5 m?
		# H = head (difference between spillway height and reservoir stage) [m]
		C=(2/3)^1.5*9.806^0.5	
		L=0.005#0.25	# m
		if((stage-stageOut)>0){
			H=stage-stageOut
			Q=C*L*H^1.5	#m3 s-1
			STout=Q*60*60*24		#m3 day-1
		}else{
			STout=0		#m3 day-1
		}

		################################
		#### Differential equations ####
		################################

		dV.dt=(streamQ+directPrecip+gwIn)-(STout+curEvap+gwOut)
		
		list(c(dV.dt))
	})

}
