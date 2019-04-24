source("~/Dropbox/Rtools/gray_model.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/plot_tools.R")
load("~/Dropbox/18cts/data/gray_pre_rce.Rdata")

#=======#
# Data  #
#=======#

# Data params
Ngamma    = length(gammalist) 
Np	  	  = length(pvals)
taus      = 20

#=======#
# Plot  #
#=======#

# plot params
cex	  	  = 2.5
cex_legend = 2.5
lwd	      = 3.5
ylim_p    = 1e-2*rev(range(pvals))
ylim_tau  = c(taus,0)
xlim_p    = c(-12,12.5)  # K/day
xlim_tau  = c(-0.8,0.8)   # pptau F/OLR
xlab_p    = "H (K/day)"
xlab_tau  = expression("("*partialdiff[tau]*F*")"/OLR)
ylab_p    = expression("p (hPa)")
ylab_tau  = expression(bold(tau))
#termlist  = c("cts","sx","ax","gx","sum")
termlist  = c("cts","sx","ax","gx","pptauf")
termnames = c("CTS","SX","AX","GX","sum")
cols      = c("blue","forestgreen","orange","red","black")
ltys      = c(rep("solid",times=4),"dashed")
lwds      = c(rep(3.5,times=4),4)
Nterms    = length(termlist)
varlist   = c(termlist,"tauvals_s","tauvals","OLR")
coo1d_fac = g/Cp*3600*24   # pppf to K/day
lty_tau1  = "dashed"
lwd_tau1  = 2

plot_points_tau = function(profile,tauvals,tauvals_s,col,lty="solid",lwd=3.5){
		points(profile[-Np]/OLR,tauvals_s[-Np],type="l",lty=lty,lwd=lwd,col=col)		
}

plot_points_p = function(profile,tauvals,tauvals_s,col,lty="solid",lwd=3.5){
		points(coo1d_fac*diff(tauvals)/diff(pvals)*profile[-Np],
			1e-2*pvals_s[-Np],type="l",lty=lty,lwd=lwd,col=col)		
}


#============#
# Plot loops #
#============#

file = "~/Dropbox/18cts/plots/gray_pre_rce_decomp.pdf"
pdf(file=file,width=14,height=10)
par(mfrow=c(2,3),oma=c(0.5,8,4,3),mar=c(5,0,3,0),tcl=-0.75,yaxs="r")

for (k_coord in 1:2){
	coord = c("tau","p")[k_coord]
	for (i in 1:Ngamma){
		# Data
		gamma	= gammalist[i]
		if (is.na(gamma)){
		    case    = paste("pre_beta",beta_pre,"_taus",taus,sep="")
			casetype = "PRE"
		} else if (!is.na(gamma)){
		    case    = paste("rce_gamma",gamma,"_taus",taus,sep="")			
			casetype = "RCE"
		}
	    output  = gray_calcs[[case]]
	    for (var in varlist){
			assign(var,output[[var]])
	    }    
	 	sum = cts + gx + sx +ax

		# plot frame
		for (var in c("xlab","ylab","xlim","ylim","plot_points")){
			assign(var,eval(as.name(paste(var,coord,sep="_"))))
		}
	    plot(1,type="n",axes=TRUE,
	    	xlim = xlim,
	    	ylim = ylim,
			xlab = "",
			ylab = "",
			yaxp = c(ylim[1],ylim[2],5),
			col.axis = "white",
			cex.axis = cex,
			cex.main = cex,
			cex.lab = cex
	    )
	    mtext(xlab,side=1,line=3.5,cex=1.75)
		main = bquote(.(casetype)*", "~gamma==.(gamma))
	    axis(1,cex.axis=cex)
	    if (is.na(gamma)){
			main = bquote(.(casetype))
	    	axis(2,at=ylim[1]*c(1,0.8,0.6,0.4,0.2,0),cex.axis=cex)
	    	mtext(ylab,2,line=4.5,cex=2)
	    }
		if (coord=="tau"){
			mtext(main,3,line=2,cex=2)
		    abline(h=1,lty=lty_tau1,col="gray",lwd=lwd_tau1)
		}
		panel = paste("",letters[3*(k_coord-1)+i],"",sep="")
		mtext(panel,3,line=-2,adj=0.02,cex=1.5,col="darkgray")

		# profiles
	    for (n in 1:Nterms){
			col	= cols[n]
			lty	= ltys[n]
			lwd	= lwds[n]
			term = termlist[n]
			profile = eval(as.name(term))
			plot_points(profile,tauvals,tauvals_s,col=col,lty=lty,lwd=lwd)
	    }


		# p_tau1
	    p_tau1 = pvals[which.min(abs(tauvals - 1))]		
	    abline(h=p_tau1*1e-2,lty=lty_tau1,col="gray",lwd=lwd_tau1)
	    #text(xlim[1]+0.5,1e-2*p_tau1-30,expression(tau==1),cex=1)		
		
		# Add legend
		if (coord=="tau" & i==1){
	         legend("left",legend=c(termnames[1:5]),col=cols,
			    lty=ltys,lwd=lwd,cex=cex_legend,bty="")
		 }

	}  # gamma
} # coord
dev.off()