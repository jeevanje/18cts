Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/18cts/src/")

library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"gray_model.R",sep=""))
source("plot_params.R")

# data params
cases = c("h2o_only_no_cont","co2_only_simple_atm")
D     = 1.5  # diffusion parameter

# plot params
coo_lim = c(-0.015,0.002) # K/day/cm^-1
Hlim_h2o = c(-3,0)
Hlim_co2 = c(-0.7,0.2)
fields 	= c("coo","sx","ax","gex","cts","sum")
fieldnames = c(expression(H[k]),"SX","AX","GX","CTS","sum")
N_fields= length(fields) # no sum
ltyvec	= c("dashed",rep("solid",times=4),"solid")
colvec  = c("black","forestgreen","orange","red","blue","black")
lwd     = 4
lwdvec  = c(rep(lwd,times=5),lwd-2.5)
cex		= 1.75
cex_leg = 1.5
plab    = "p (hPa)"
p_lim   = c(1000,0.1)
Hklab   = expression(H[k]~~"("*K/day/cm^{-1}*")")
gasnames = c(expression(H[2]*O),expression(C*O[2]))

# PDF
file = "../plots/both_taus20_profiles.pdf"
pdf(file=file,width=10,height=6,bg="white")
par(mfrow=c(1,2),mar=c(5,5,5,4))

for (i in 2:1){
	case = cases[i]

	#=======#
	# Data  #
	#=======#
	
	# Get data
	ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
	nc      = nc_open(ncpath)
	k       = ncvar_get(nc,"k")
	dk		= k[2]-k[1]
	p       = ncvar_get(nc,"p")
	np      = length(p)
	tabs    = ncvar_get(nc,"tabs")
	nk		= length(k)
	coo2d   = ncvar_get(nc,"coo")/coo2d_fac     # convert to SI 
	opt     = D*ncvar_get(nc,"opt")
	flx2d   = ncvar_get(nc,"flx")*1e-2    # W/m^2/m-1

	#=========#
	# Process #
	#=========#
	
	p_s	 	 = rfm_lni2s(p)
	tabs_s   = rfm_i2s(tabs)
	coo2d_s  = rfm_i2s(coo2d)
	opt_s	 = rfm_lni2s(opt)
	dtaudp2d = (opt[ ,2:np]-opt[ ,1:(np-1)])/(rep(1,times=nk)%o%diff(p)) #s lev
	B_i      = outer(k,tabs,planck_k)
	B_s      = outer(k,tabs_s,planck_k)
	opt_surf = opt[ ,1]%o%(rep(1,times=(np-1)))
	B_surf   = B_i[ ,1]%o%(rep(1,times=(np-1)))
	
	cts2d    = calc_cts2d(k,p,tabs_s,opt)
	gex2d    = -pi*(B_surf-B_s)*exp(-(opt_surf - opt_s))*dtaudp2d  #s, z_s
	pppf2d   = (flx2d[ ,2:np]-flx2d[ ,1:(np-1)])/(rep(1,times=nk)%o%diff(p)) #slev
	cts1d    = apply(cts2d,2,sum)*dk
	coo1d    = -apply(pppf2d,2,sum)*dk
	
	#=======#
	# Plot  #
	#=======#
	gas     = substr(case,1,3)
	gasname = gasnames[i]

	#p1vals  = 1e2*c(50,500,800)
	p1vals  = 1e2*c(300,700)
	tausvals= c(20)
	Np1     = length(tausvals)
	kmin    = 3e4
	mmin    = min(which(k>kmin))

    p1vec   = numeric(nk)
    for (m in 1:nk){
        p1vec[m] = p_s[which.min(abs(opt_s[m,]-1))]
    }
	
	for (n in 1:Np1){
	    taus    = tausvals[n]
	    m 	    = which.min(abs(opt[(k>kmin)&(k<8e4),1]-taus)) + mmin - 1
		p1      = p1vec[m]
	    kval    = k[m] 
	    tauvals = opt[m,1:(np-1)]  #i
	    tauvals_s = opt_s[m,1:(np-1)]  #s
	    Bvals   = B_s[m,] 	       #s
	    dtaudp  = dtaudp2d[m,]     #s
	    pppf    = pppf2d[m,]
	    coo   = -pppf 	# use simple difference rather than RFM's 3-stencil
	    cts   = cts2d[m,]
	    gex   = gex2d[m,]
	    taus  = round(tauvals[1],digits=1)
	
	    # gray_model.R
	    ex    = compute_ex(pi*Bvals,tauvals)  #B s levs, tau i levs
	    ax    = (ex$ax)*dtaudp
	    sx    = (ex$sx)*dtaudp
	    sum   = cts+gex+ax+sx
	
	    plot(1,type="n",
		    ylim = p_lim,
		    xlim = coo_lim,
		    xaxt = "n",
		    xlab = "",
		    ylab = plab,
		    main = gasname,
		    cex.lab  = cex,
		    cex.axis = cex,
		    cex.main = cex
	    )
	    axis(1,at=c(0,-0.01),labels=c("0","-0.01"),
	    	 cex.axis=cex,lwd.ticks=2)
		mtext(Hklab,side=1,line=3.5,cex=cex)
	    for (i in 1:N_fields){
	    	temp = eval(as.name(fields[i]))
			lty  = ltyvec[i]
			col  = colvec[i]
			lwd  = lwdvec[i]
			points(coo2d_fac*temp,1e-2*p_s,type="l",lwd=lwd,lty=lty,col=col)
	    } # field loop
	    abline(h=1e-2*p1,col="gray",lty="dashed")
	}  # mvals loop
	if (gas=="co2"){
		legend("bottomleft",legend=fieldnames[1:N_fields],cex=cex_leg, 
		      lty=ltyvec,col=colvec,lwd=lwdvec,xpd="TRUE")
	}
} #case loop
dev.off()
