args = commandArgs(trailingOnly=TRUE)

library(fields)
library(ncdf4)
Rtoolsdir = "~/Dropbox/Rtools/"
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))
source(paste(Rtoolsdir,"gray_model.R",sep=""))

gas    = args[1]
case_h2o   = "h2o_only_no_cont"
case_h2o   = "h2o_noctm_Ts300_rh0.75_gamma7"
case_co2   = "co2_only_simple_atm"
case_co2   = "co2_noctm_Ts300_rh0.75_gamma7"
Hklim_h2o = c(-0.9e-2,0.5e-2) # K/day/cm^-1
Hklim_co2 = 0.65*Hklim_h2o
Hlim_h2o   = c(-3.8,0.5)   # K/day
Hlim_co2   = c(-1,0.25)   # K/day
xunits_h2o = 1720
xunits_co2 = 900
xaxis_h2o  = c(200,400,600,800,1000,1400)
xaxis_co2  = NULL
p1lwd_h2o  = 1.15
p1lwd_co2  = 1.5
p1col_h2o  = "lightgray"
p1col_co2  = "gray"


for (var in c("case","Hklim","Hlim","xunits","xaxis","p1lwd","p1col")){
	assign(var,eval(as.name(paste(var,gas,sep="_"))))
}

D      = 1.5  # diffusion parameter
Nsec   = 3600*24   # sec/day
fac2d  = g/Cp*Nsec*1e2  # pppf to K/day/cm^-1
fac1d  = -g/Cp*Nsec      # pppf to K/day
ncoarse = 10
fields = c("gx","ax","sx","cts","coo")
fieldnames = c("GX","AX","SX","CTS","H")

#=======#
# Data  #
#=======#

ncpath  = paste("~/Dropbox/17rad_cooling2/data_paper/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k")
dk      = k[2]-k[1]
z       = ncvar_get(nc,"z")
p       = ncvar_get(nc,"p")
tabs    = ncvar_get(nc,"tabs")
nz      = length(z)
nk      = length(k)
pvec    = (nz-1):1
coo2d   = ncvar_get(nc,"coo")/fac2d     # convert to SI 
opt     = D*ncvar_get(nc,"opt")
flx2d   = ncvar_get(nc,"flx")*1e-2      # W/m^2/m^-1
ncoarse = 1000/dk

#=========#
# Process #
#=========#

p_s      = rfm_lni2s(p)
z_s      = rfm_i2s(z)
tabs_s   = rfm_i2s(tabs)
coo2d_s  = rfm_i2s(coo2d)
opt_s    = rfm_lni2s(opt)
dtaudp2d = (opt[ ,2:nz]-opt[ ,1:(nz-1)])/(rep(1,times=nk)%o%diff(p)) #s lev
B_i      = outer(k,tabs,planck_k)
B_s      = outer(k,tabs_s,planck_k)
opt_surf = opt[ ,1]%o%(rep(1,times=(nz-1)))
B_surf   = B_i[ ,1]%o%(rep(1,times=(nz-1)))
ktp      = min(which(tabs==min(tabs)))
ptp      = p[ktp]

# 2d terms
cts2d_s  = pi*B_s*exp(-opt_s)*dtaudp2d  #s, z_s, W/m^2/Pa/m^-1, pppf units
gx2d_s  = -pi*(B_surf-B_s)*exp(-(opt_surf - opt_s))*dtaudp2d  #s, z_s
pppf2d   = (flx2d[ ,2:nz]-flx2d[ ,1:(nz-1)])/(rep(1,times=nk)%o%diff(p)) #s lev
#coo2d_s = -pppf2d   # W/m^2/Pa/m^-1
ax2d_s   = array(dim=c(nk,nz-1))  # s
sx2d_s   = array(dim=c(nk,nz-1))  # s

for (m in 1:nk){
    tauvals = opt[m,1:(nz-1)]  #i
    Bvals   = B_s[m,] 	       #s
    dtaudp  = dtaudp2d[m,]     #s    
    ex    = compute_ex(pi*Bvals,tauvals)  #B s levs, tau i levs
    ax    = (ex$ax)*dtaudp
    sx    = (ex$sx)*dtaudp
    ax2d_s[m, ] = ax
    sx2d_s[m, ] = sx
    }
sum2d_s = cts2d_s + gx2d_s + ax2d_s + sx2d_s

for (field in fields){
    field2d = eval(as.name(paste(field,"2d_s",sep="")))
    assign(paste(field,"1d",sep=""), apply(field2d,2,sum)*dk)  # W/m^2/Pa
}

# tau=1
p1  = array(dim=c(nk))
for (m in 1:nk){
    p1[m] = p_s[which.min(abs(log(opt_s[m,])))]
}
p1_coarse = coarse_grain(p1,ncoarse)

#=======#
# Plot  #
#=======#

# Plot params
lwd	   = 3
ltyvec = c(rep("solid",times=4),"dashed")
colvec = c("red","orange","green","blue","black")
Ncol       = 15
xvals      = c(Hklim[1]*seq(1,1/(Ncol-1)/2,length=(Ncol-1)/2),Hklim[2]*seq(0,1,length=(Ncol+1)/2))
colbar	   = designer.colors(col=tim.colors(Ncol),x=xvals)
cex	       = 2.5
cex_lab   = 1.75
plim	   = c(1000,75)
Hk_units   = expression(K/d*a*y/cm^{-1})
col_panel  = "black"

# PDF
file = paste("~/Dropbox/18cts/plots/realgas_decomp_",gas,".pdf",sep="")
pdf(file,width=14.9,height=10,bg="white")
par(oma=c(0,6,0,5),tcl=-0.65)
set.panel(2,3)
for (k_field in 1:5){
	field	= fields[k_field]
	fieldname = fieldnames[k_field]
    varname = paste(field,"2d_s",sep="")
    var = eval(as.name(varname))
    var_coarse = coarse_grain(var,ncoarse)
    k_coarse   = coarse_grain(k,ncoarse)


	par(mar=c(5,0,5,0))
	if (field=="sx"){
		add.legend=TRUE
		par(mar=c(5,0,5,5))
	} else {
		add.legend=FALSE
	}
    my.image.plot(1e-2*k_coarse,1e-2*p_s[pvec],-fac2d*var_coarse[ ,pvec],
		add.legend=add.legend,
		ylim = plim,
		zlim = Hklim,
		xlab = "",
		ylab = "",
		col  = colbar,
		main = fieldname,
		cex.lab  = cex,
	    cex.axis = 1e-5,
		cex.main = cex+0.25,
		cex.legend = cex,
		legend.width=1.5
		)
	# Add axes
	title("")
	axis(1,at=c(0,1500), labels=c("",""), lwd.ticks=0)
	axis(1,at=xaxis,cex.axis=cex)
	mtext(klab,1,line=3.5,cex=cex_lab)
	if ( (field=="gx" | field=="cts")){
		axis(2,cex.axis=cex)
		title("")
		mtext("Pressure (hPa)",2,line=3.5,cex=cex_lab)
	}
	

	# add units
	if (field=="sx"){
		title("")
		mtext(Hk_units,3,line=0,at=xunits,cex=1.5)
		abline(h=1e-2*ptp,lty="dotted",col="gray",lwd=2)
	}

	#Add panel label
    mtext(letters[k_field],1,line=-1.5,adj=0.02,cex=1.5,col=col_panel)

	# Add p1
	points(1e-2*k_coarse,1e-2*p1_coarse,type="l",lty="dashed",lwd=p1lwd,col=p1col)	
}

# H profiles
par(mar=c(5,0,5,5))
plot(1,type="n", xlim = Hlim, ylim=plim,
	axes	 = TRUE,
	xlab 	 = "Heating (K/day)",	 
	ylab     = plab,
	cex.lab  = cex,
	cex.axis = 1e-5,
	cex.main = cex)
#axis(2,cex.axis=1e-5)
axis(1,cex.axis=cex)
abline(v=0,lty="dotted",col="gray",lwd=2)
abline(h=1e-2*ptp,lty="dotted",col="gray",lwd=2)
for (n in 1:length(fields)){
    field = fields[n]
    lty   = ltyvec[n]
    col   = colvec[n]
    var   = eval(as.name(paste(field,"1d",sep="")))
    points(fac1d*var,1e-2*p_s,type="l",lwd=lwd,lty=lty,col=col)
}
legend("topleft",legend=fieldnames,lty=ltyvec,lwd=lwd,col=colvec,cex=cex-0.25)
mtext(letters[6],1,line=-1.5,adj=0.02,cex=1.5,col=col_panel)

dev.off()
