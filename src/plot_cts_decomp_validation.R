source("~/Dropbox/Rtools/gray_model.R")
source("~/Dropbox/Rtools/thermo_tools.R")
load("~/Dropbox/17rad_cooling2/data/gray_calcs.Rdata")

Nbeta     = length(betalist)
Np	  	  = length(pvals)
taus      = 50
logvec    = c("","")
cex	  	  = 1.5
plim      = 1e-2*rev(range(pvals))
Hlims     = list(c(-4,8),c(-4,8))
termlist  = c("cts","sx","ax","gex","sum")
colvec    = c("blue","green","orange","red","black")
Nterms    = length(termlist)
varlist   = c(termlist,"pppf","tauvals")
Nsec	  = 3600*24   # sec/day

plot_points = function(profile,col){
		points(g/Cp*Nsec*diff(tauvals)/diff(pvals)*profile[-Np],
			1e-2*pvals[-Np],type="l",
			lty="dashed",lwd = 2  ,col=col)		
}

#============#
# Plot loop  #
#============#
file = paste("~/Dropbox/17rad_cooling2/plots/cts_decomp_validation_taus",taus,".pdf",sep="")
pdf(file=file,width=9,height=5)
par(mfrow=c(1,2),mar=c(5,5,5,3))
for (i in 1:Nbeta){
	beta	= betalist[i]
	Hlim	= Hlims[[i]]
	log		= logvec[i]
    case    = paste("beta",beta,"_taus",taus,sep="")
    output  = gray_calcs[[case]]
    for (var in varlist){
		assign(var,output[[var]])
    }    
 	sum = cts + gex + sx +ax
    plot(1,type="n",xlim=Hlim,ylim=plim,
		xlab = "H (K/day)", 
		ylab = "p (hPa)",
		main = bquote(beta==.(beta)*",  "~tau[s]==.(taus)),
		log	 = log,
		cex.axis = cex,
		cex.main = cex,
		cex.lab = cex
    )
    if (beta ==betalist[1]){
       legend("bottomright",legend=c("H",termlist),col=c("black",colvec),
		   lty=c("solid",rep("dashed",times=Nterms)),lwd=2,cex=1.2,
		   bty="n")
    }	   
    points(-g/Cp*Nsec*pppf,1e-2*pvals_s[-Np],type="l",lwd = 1,col="black")	
    for (n in 1:Nterms){
		col	= colvec[n]
		profile = eval(as.name(termlist[n]))
		plot_points(profile,col)
    }
    p_tau1 = pvals[which.min(abs(tauvals - 1))]		
    abline(h=p_tau1*1e-2,lty="dotted",col="black",lwd=1.5)
    text(Hlim[1]+0.5,1e-2*p_tau1-30,expression(tau==1),cex=1)		
}
dev.off()