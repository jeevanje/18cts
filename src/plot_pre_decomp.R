
# All field normalized by OLR!
terms_new = c("CTS","SX","AX","GX","sum_new")
terms_old = c("CTS","EX_above","EX_below","GX","sum_old")
names_new = c(terms_new[1:4],"sum")
names_old = c("CTS",expression(E*X[above]),expression(E*X[below]),"GX","sum")
cols     = c("blue","forestgreen","orange","red","black")
ltys     = c(rep("solid",times=4),"dashed")
vars     = c(terms_new,terms_old,"tau")

tausvals = c(20)
taus     = tausvals[1]
Ntau	 = 200
Ntaus	 = length(tausvals)
Nterms	 = length(terms_new)

for (taus in tausvals){
	tau 	  	= seq(0,taus,length=Ntau)
	Bvals    	= 0.5*(1+tau)
	#Bvals[Ntau] = 0.5*(2+taus)

	CTS	= -Bvals*exp(-tau)
	SX  = rep(0,times=Ntau)
	AX	= 0.5*(-(taus-tau)*exp(-(taus-tau)) + tau*exp(-tau) - 
				exp(-(taus-tau)) + exp(-tau))
	GX  = 0.5*(1+taus-tau)*exp(-(taus-tau))
	EX_above = -0.5*(-(1+tau)*exp(-tau)+1)
	EX_below = 0.5*(-(taus-tau+1)*exp(-(taus-tau))+1)
	sum_new = CTS+SX+AX+GX
	sum_old = CTS+EX_above+EX_below+GX
		
	for (var in vars){
		assign(paste(var,taus,sep="_"),eval(as.name(var)))
	}
}


# Plot
divlab = expression("("*partialdiff[tau]*F*")"/O*L*R)
taulab = expression(tau)
lwd	   = 2
cex	   = 1.25
plot_decomp = function(taus,decomp){
	for (var in vars){
		assign(var,eval(as.name(paste(var,taus,sep="_"))))
	}	
#	divlim = range(c(AX,SX,CTS,GX))
	divlim = c(-0.75,0.75)
	taulim = rev(range(tau))
	plot(1,type="n",xlim=divlim,ylim=taulim,
			xlab = divlab,
			ylab = taulab,
			main = paste(decomp," decomposition",sep=""),
			cex.main = cex,
			cex.lab  = cex,
			cex.axis = cex
	)
	for (n in 1:(Nterms)){
		col = cols[n]
		lty = ltys[n]
		terms = eval(as.name(paste("terms_",decomp,sep="")))
		term  = eval(as.name(terms[n]))
		points(term,tau,type="l",lty=lty,col=col,lwd=lwd)
	}
}

pdf("~/Dropbox/18cts/plots/pre_decomp.pdf",width=8,height=4.5)
par(mfrow=c(1,2),mar=c(7,6,4,2))
for (decomp in c("old","new")){
	terms = eval(as.name(paste("terms_",decomp,sep="")))
	names = eval(as.name(paste("names_",decomp,sep="")))
	plot_decomp(taus,decomp)
	legend(-1.5,23,names,lwd=lwd,cex=0.8,col=cols,lty=ltys,xpd="TRUE")
}
dev.off()