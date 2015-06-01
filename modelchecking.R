summ=matrix(NA,nrow(outConc),5)

PPmax=5		# hr-1
	kSRP=0.25/1000/31#0.5/1000	# mol P m-3; Halmann and Stiller 1974 L&O 19(5): 774-783 (Lake Kinnerrett)
	Ik=180 		# light limitation benchmark for GPP, [umol cm-2 s-1]; Vadeboncoeur et al. 2008
			
for(j in 1:nrow(summ)){
	chlCur=out[j,6]/outHydro[j,7]*12*1000/phyto_C2Chl
	kD=0.0213+0.0177*chlCur+0.0514*(out[j,4]/outHydro[j,2]*12)
	
	d=5
	settling=0.0188*(d/2)^2/outHydro[j,6]
	
	summ[j,1]=dailyGPP(day=j,curForce=curForce,curForceDOY=curForceDOY,chlCur=chlCur,SRP=out[j,7],DOC=out[j,4],V=outHydro[j,2],kD=kD,zmix=outHydro[j,6])*12*1000
	summ[j,2]=((out[j,7]/out[j,2])/((out[j,7]/out[j,2])+kSRP))
	
	hourlyPAR=sw.to.par.base(curForce[curForceDOY==j,10])		# umol m2 sec; SW to PAR based on 
	avgI=lightAttenTS(hourlyPAR,kD=kD,zmix=outHydro[j,6])		# hourly average light climate in mixed layer [umol cm-2 
	Ilim=tanh(avgI/Ik)
	summ[j,3]=mean(Ilim)
	summ[j,4]=min(Ilim)
	summ[j,5]=max(Ilim)
}

colnames(summ)=c("GPP_mgCm-3d-1","Plimitation","meanIlim","minIlim","maxIlim")