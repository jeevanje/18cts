source("~/Dropbox/Rtools/gray_model.R")
source("~/Dropbox/Rtools/thermo_tools.R")

# Fixed params & profiles
Ts	       = 300
Ttp	       = 210
Bs	       = B(Ts)
Btp	       = B(Ttp)
pmin  	   = 1e2   # pa
deltap     = 250*pmin
pvals	   = rev(seq(pmin,ps,by=pmin))  # pa, i levels
Np		   = length(pvals)
pvals_s    = c(pvals[-Np] + diff(pvals)/2,pvals[Np]+diff(pvals)[Np-1]/2) 
Gamma	   = 7e-3
Bvals_adbt = B(Ts*(pvals_s/ps)^(Rd*Gamma/g))  # s levels
Bvals_strat= Btp
ktp        = min(which(Bvals_adbt<Bvals_strat))
ptp        = pvals_s[ktp]
Bvals      = Bvals_adbt/2*(tanh((pvals_s-ptp)/deltap)+1) - 
			 Bvals_strat/2*(tanh((pvals_s-ptp)/deltap)-1)
Bvals      = Bvals_adbt

# case specifications
alpha      = 4
beta_pre   = 2
gammalist  = c(NA,0.5,0.1)
Ngamma     = length(gammalist)
tauslist   = c(5,10,15,20)

# Output
outputlist = c("gamma","beta","tauvals","tauvals_s","OLR","Q",
				"pptauf","cts","sx","ax","gx")
gray_calcs = list()

#============#
# Run loops  #
#============#
for (taus in tauslist){
	for (k_gamma in Ngamma:1){ # so RCE calculated first, for OLR
	    output  = list()
		gamma   = gammalist[k_gamma]
		if (!is.na(gamma)){
			# RCE
			beta    = alpha*Rd*Gamma/g/gamma
			case	= paste("rce_gamma",gamma,"_taus",taus,sep="")
			tauvals   = taus*(pvals/ps)^beta  # i levels
			tauvals_s = c(tauvals[-Np] + diff(tauvals)/2,
				      tauvals[Np]+diff(tauvals)[Np-1]/2) 
			U	    = compute_Ugray(tauvals,Bvals,Bs)
			D	    = compute_Dgray(tauvals,Bvals)
			Fnet	= U - D
			pptauf  = diff(Fnet)/diff(c(tauvals,0))     # s, up to Np !
			OLR     = Fnet[Np+1]
			Q	    = Fnet[Np+1] - Fnet[1]
			cts     = -Bvals*exp(-tauvals_s)
			gx      = (Bs-Bvals)*exp(-(taus - tauvals_s))
			ex		= compute_ex(Bvals,tauvals)
			ax	    = -ex$ax
			sx	    = -ex$sx
		} else if (is.na(gamma)){
			# PRE
			beta    = beta_pre
			case	= paste("pre_beta",beta,"_taus",taus,sep="")
			tauvals = taus*(pvals/ps)^beta  
			tauvals_s = taus*(pvals_s/ps)^beta  
			tau     = tauvals_s  # for convenience
			OLR     = 0.5*(gray_calcs$rce_gamma0.1_taus20$OLR + 
						   gray_calcs$rce_gamma0.5_taus20$OLR )			
			cts		= -OLR*0.5*(1+tau)*exp(-tau)
			sx	    = rep(0,times=Np)
			pptauf  = rep(0,times=Np)
			ax		= OLR*0.5*(-(taus-tau)*exp(-(taus-tau)) + tau*exp(-tau) - 
						exp(-(taus-tau)) + exp(-tau))
			gx	    = OLR*0.5*(1+taus-tau)*exp(-(taus-tau))
			Q       = 0
		}			
		for (var in outputlist){
		   output[[var]] <- eval(as.name(var))
		}
		gray_calcs[[case]] <- output 
	}	# k_gamma
}	# taus
save(pvals,pvals_s,gammalist,beta_pre,tauslist,gray_calcs,file="~/Dropbox/18cts/data/gray_pre_rce.Rdata")