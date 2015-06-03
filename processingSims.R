apply(out,2,range)
range(out[,4]/out[,2]*12)

dev.new()
par(mfrow=c(2,1))
plot(out[,4]/out[,2]*12,type='l')
mean(out[2475:3075,4]/out[2475:3075,2]*12)


stage=((3*out[,2])/(r2h^2*pi))^(1/3)
A=pi*(r2h*stage)^2
		
perim=2*pi*sqrt(A/pi)*DL			#### maybe just treat this as a circle???

zmix=10^(-0.518*log10(out[,4]/out[,2]*12)+1.006)

r1=r2h*stage
r2=r2h*(stage-zmix)
Vepi=(1/3)*pi*(r1^2+r1*r2+r2^2)*zmix

plot(out[,6]/Vepi*12/50*1000,type='l')
mean(out[2475:3075,6]/Vepi[2475:3075]*12/50*1000)



(20*UNDERCsheds$Area_m2/1000)/(10^(-0.0589+1.12963*log10(UNDERCsheds$NHLD_lakes)))

#Lake name		i		runoff mode		proportion of volume represented by max daily stream load
#Reddington 	3		variable			0.16
#MOrris			4		variable			0.11
#Mullahy		5		constant			0.2
#Ward			8		constant			0.17
#Inkpot			20		constant			0.19
#TF				21		variable?		0.11
#NG				27		constant			15