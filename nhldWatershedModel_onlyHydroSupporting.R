# nhldWatershedModel_onlyHydro supporting functions

timeStep<-function(t,X,params){
	with(as.list(params),{
			
		#########################
		#### State variables ####
		#########################
		V=X[1]
		stage=V/curLakeArea	
				
		streamQ=dailyRunoff(t)*curShedArea/1000			#m3 day-1	****** problem because of interpolating?
		directPrecip=dailyPrecip(t)*curLakeArea/1000		#m3 day-1	****** problem because of interpolating?
			
		curEvap=dailyEvap(t)*curLakeArea/1000			#m3 day-1	****** problem because of interpolating?
		
		# assuming ogee crest spillway
		# Q=C*L*H^(3/2), where:
		# Q = discharge [m3 s-1]
		# C = coefficient; (2/3)^1.5*g^0.5 [m^(1/2) s-1]; g = gravitational constant; 9.806 [m s-2]
		# L = spillway length [m]; call this 2-5 m?
		# H = head (difference between spillway height and reservoir stage) [m]
		C=(2/3)^1.5*9.806^0.5	
		L=0.5	# m
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
