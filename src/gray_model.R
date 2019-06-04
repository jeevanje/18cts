#================#
# Functions for  #
# computing      #
# gray radiative #
# fluxes         #
#================#

Rtoolsdir = "~/Dropbox/Rtools"
source(paste(Rtoolsdir,"/calculus_tools.R",sep=""))
source(paste(Rtoolsdir,"/thermo_tools.R",sep=""))

#===============#
# Optical depth #
# routines      #
#===============#

# WVP as a function of k, for tau calculation! 
compute_wvp_k   = function(z,rhov,lev="s"){
	    # wvp on dam i grid, z on dam lev grid, rhov on dam s grid
		if (lev=="i"){
			dz   = diff(z)
		} else if (lev=="s"){
		    zint = zinterp(z)
		    dz   = diff(zint)  # length = nz-1			
		}
	    nz   = length(z)
	    wvp  = numeric(nz)
	    wvp[nz] = rhov[nz]*dz[nz-1]  # extrapolate dzvec
	    for (k in (nz-1):1){
                 wvp[k] = wvp[k+1] + dz[k]*rhov[k]
            } 
	    return(wvp)
}

compute_trans = function(tau){
	        # All quantities ssi!
	        nz = length(tau)
	      	trans = array(dim=c(nz,nz))
		trans = outer(exp(-tau),exp(tau),"*")
		}

#====================#
# Radiative transfer #
#====================#

compute_Ugray = function(tau,Bvals,Bs=Bvals[1]){
	    # tau dam ssi, Bvals dam sss
	    ntau  = length(tau)
	    dtau  = diff(c(tau,0))
	    Ugray = numeric(ntau+1)  #i levels + 1!
   	    Ugray[1] = Bs
            for (i in 1:(ntau)){
      	       # see rad_cooling notes 10/17
	       Ugray[i+1] = ( (1+dtau[i]/2)*Ugray[i] -dtau[i]*Bvals[i] ) /(1-dtau[i]/2)
            }
            dim(Ugray) <- ntau+1		    
            return(Ugray)  
	}

compute_Dgray = function(tau,Bvals){
	    # tau ssi, Bvals sss
	    ntau  = length(tau)
	    dtau  = diff(c(tau,0))
	    Dgray = numeric(ntau+1)  # i levels + 1!
	    Dgray[ntau+1] = 0  # neglect emission from above topmost grid cell
	    for (i in (ntau):1){
		# see rad_cooling notes 10/17
		Dgray[i] = ( (1+dtau[i]/2)*Dgray[i+1] -dtau[i]*Bvals[i] ) /(1 - dtau[i]/2)
		}
	    dim(Dgray) <- ntau+1	
            return(Dgray)  
	}

#===============#
# Heating rates #
#===============#

# compute_ex outputs profiles of the ax and sx exchange terms, given 
# source function (on model levels) and optical depth profiles (on interface
# levels)

compute_ex   = function(Bvals,tau){
	       # tau ssi, Bvals sss	       
   	       epsilon=1e-10  # fudge for weird R bug in k_2tau
	       N     = length(tau)
	       ex    = list()
	       ex[['ax']] <- numeric(N)
	       ex[['sx']] <- numeric(N)		   
	       taus  = tau[1]
	       tau_s = c(tau[-N] + diff(tau)/2,tau[N]+diff(tau)[N-1]/2)
	       dtau  = -c(diff(tau),0) # slight fudge here
	       Bvals_i = s2i(1,tau,Bvals,f0=Bvals[1]+0.5*diff(Bvals)[1])
	       for (k in 1:N){   # k is s levels
 	       	   if (tau_s[k] < taus/2){

		      # ax		
	       	      k_2tau = max(which(tau>2*tau_s[k]))  # i level
	   	      kvec   = 1:(k_2tau-1)  #s levs
	       	      ax =  -sum(dtau[kvec]*(Bvals[kvec]-Bvals[k])*exp(-(tau_s[kvec]-tau_s[k])))
		      ex[['ax']][k] <- ax
				   
	 	      # sx	
		      kvec1  = (k+1):N
	       	      sx1    = -sum(dtau[kvec1]*(Bvals[kvec1]-Bvals[k])*exp(-(tau_s[k]-tau_s[kvec1])))
		      kvec2  = k_2tau:(k-1)	
		      sx2    = -sum(dtau[kvec2]*(Bvals[kvec2]-Bvals[k])*exp(-(tau_s[kvec2]-tau_s[k])))		       	      
		      gamma  = log(Bvals_i[k+1]/Bvals_i[k])/log(tau[k+1]/tau[k])
		      sx3    = dtau[k]/2*exp(-dtau[k]/4)*(Bvals_i[k]*(1-dtau[k]/4/tau[k])^gamma -2*Bvals_i[k]*(1-dtau[k]/2/tau[k])^gamma + Bvals_i[k+1]*(1+dtau[k]/4/tau[k+1])^gamma )
 		      ex[['sx']][k]  <- sx1 + sx2 + sx3
				  	
		   } else if (tau_s[k]>=taus/2){
		      # ax	
	       	      k_2tau = min(which(tau <= (2*tau_s[k]-taus + epsilon)))  # i level
		      kvec   = k_2tau:N  #s levs
	       	      ex[['ax']][k]  <- -sum(dtau[kvec]*(Bvals[kvec]-Bvals[k])*exp(-(tau_s[k]-tau_s[kvec])))
	       	      
	       	      # sx
		      kvec1  = 1:(k-1)
	       	      sx1    = -sum(dtau[kvec1]*(Bvals[kvec1]-Bvals[k])*exp(-(tau_s[kvec1]-tau_s[k])))
		      kvec2  = (k+1):(k_2tau-1)	
		      sx2    = -sum(dtau[kvec2]*(Bvals[kvec2]-Bvals[k])*exp(-(tau_s[k]-tau_s[kvec2])))		       	      
		      ex[['sx']][k]  <- sx1 + sx2
  		   }  #tau_s[k] > taus/2 		      
		}   # k
		return(ex)
}

 
compute_pptf = function(kappa,z,U,D,tabs,rhov,lapse) {
	     U_sss = i2s(3,z,U)
	     D_sss = i2s(3,z,D)
	     pptf = kappa/lapse*rhov*(U_sss + D_sss - 2*B(tabs))
	     return(pptf)
	     }

#=======#
# other #
#=======# 

tune_kappa_olr = function(wvp_k,tabs_s,Ts,olr){
	     cost  = function(kappa){
			     nz    = length(tabs_s)
		    	 Bvals = B(tabs_s) 
	 		     tau   = kappa*wvp_k  	     
	     	     U     = compute_Ugray(tau,Bvals,B(Ts))
		     	 return(U[nz]-olr)
		 }
	     kappa = uniroot(cost,interval=c(1e-2,1e1) )[[1]]
	     return(kappa)
	     }

tune_kappa_Q = function(wvp_k,tabs_s,Ts,Q){
	     cost  = function(kappa){
			     nz    = length(tabs_s)
		    	 Bvals = B(tabs_s) 
	 		     tau   = kappa*wvp_k  	     
	     	     U     = compute_Ugray(tau,Bvals,B(Ts))
	     	     D     = compute_Dgray(tau,Bvals)
				 Fnet  = U-D
		     	 return(Fnet[nz]-Fnet[1] - Q)
		 }
	     kappa = uniroot(cost,interval=c(1e-1,1e1) )[[1]]
	     return(kappa)
	     }
