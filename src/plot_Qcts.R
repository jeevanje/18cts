library(ncdf4)
setwd("~/Dropbox/18cts/src")

source("../../Rtools/thermo_tools.R")
source("../../Rtools/rfm_tools.R")
source("../../Rtools/plot_tools.R")

D	    = 1.5
ncoarse = 10
vars    = c("Qk","Qkcts","OLRk")
cols    = c("black","blue","red")
Nvars   = length(vars)
cases   = c("h2o_only_no_cont","co2_only_simple_atm","h2o_only_simple_atm")
gases   = c(expression(H[2]*O),expression(CO[2]))
ylab    = expression(Q[k]*","~OLR[k]~~"("~W/m^2/c*m^{-1}~")")
ymaxs   = c(0.45,0.3) 
legendposvec = c("topright","topright")
legend  = c(expression(Q[k]),expression(Q[k]^{CTS}),expression(OLR[k]))
ltys	= c("solid","solid","dotted")
pars   <- c('plt','usr')

pdf("../plots/Qcts.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(5,5,5,3))
for (n_case in 1:2){
	# Data
	case   = cases[n_case]
	ncpath = paste("../../17rad_cooling2/data/",case,".nc",sep="")
	nc	   = nc_open(ncpath)
	flx2d  = ncvar_get(nc,"flx")*1e-2 # SI W/m^2/m^-1
	opt2d  = D*ncvar_get(nc,"opt")
	tabs   = ncvar_get(nc,"tabs")
	k      = ncvar_get(nc,"k")
	p      = ncvar_get(nc,"p")
	np	   = length(p)
	tabs_s = rfm_i2s(tabs)
	k_coarse = coarse_grain(k,ncoarse)
	
	cts2d  = calc_cts2d(k,p,tabs_s,opt2d)
	Qkcts  = -cts2d%*%diff(p)  # vertical integral
	OLRk   = flx2d[ ,np]
	Qk	   = OLRk - flx2d[ ,1]
	OLRk_coarse  = coarse_grain(OLRk,ncoarse)
	Qkcts_coarse = coarse_grain(Qkcts,ncoarse)
	
	# Plot
	lwd  = 2
	cex  = 1.25
	ylim = c(0,ymaxs[n_case])
	
	if (n_case <= 2){
		plot(1,type="n",
			xlim = range(1e-2*k),
			ylim = ylim,
			xlab = klab,
			ylab = ylab,
			main = gases[n_case],
			cex.axis=cex,
			cex.main=cex,
			cex.lab=cex)
		legend(legendposvec[n_case],legend,col=cols,lwd=lwd)
		assign(paste("par_",n_case,sep=""),
				c(list(mfg=c(1,n_case,1,2)), par(pars)))
	} else if (n_case==3){
		par(par_1)
	}
	for (n in 1:Nvars ){
		var	  = vars[n]
		col   = cols[n]
		field = eval(as.name(var))
		field_coarse = coarse_grain(field,ncoarse)
		points(1e-2*k_coarse,1e2*field_coarse,type="l",
				col=col,lwd=lwd,lty=ltys[n_case])
	}	
	
}
dev.off()
