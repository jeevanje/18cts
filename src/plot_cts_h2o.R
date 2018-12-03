
library(ncdf4)
library(fields)
Rtoolsdir = "~/Dropbox/Rtools/"
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))

case   = "h2o_only_no_cont"
D      = 1.5  # diffusion parameter
Nsec   = 3600*24   # sec/day
fac2d  = -g/Cp*Nsec*1e2  # pppf to K/day/cm^-1
fac1d  = -g/Cp*Nsec      # pppf to K/day
fields = c("coo","cts")
names1d  = c(expression(H),expression(H[CTS]))
mains2d  = c(expression("(a)"~~~~H[k]~~"("*K/day/cm^{-1}*")"),
	   	     expression("(b)"~~~~H[k]^{CTS}~~"("*K/day/cm^{-1}*")"))
ltyvec = c("solid","dashed")
colvec = two.colors(start="blue",end="red",middle="white")
plim   = c(1000,0)	  # hPa
ylim   = plim
ylab   = plab
Hlim   = c(-3,0)     # K/day
Hklim  = 0.01*c(-1,1)  # K/day/cm^-1

#=======#
# Data  #
#=======#

ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k")
dk	    = diff(k)[1]
p       = ncvar_get(nc,"p")
p_s     = rfm_lni2s(p)
tabs    = ncvar_get(nc,"tabs")
tabs_s  = rfm_i2s(tabs)
np      = length(p)
pvec    = (np-1):1
opt     = D*ncvar_get(nc,"opt")
coo2d   = ncvar_get(nc,"coo")/abs(fac2d)     # convert to SI 
coo2d_s = rfm_i2s(coo2d)
coo1d_s = apply(coo2d_s,2,sum)*dk
cts2d_s = calc_cts2d(k,p,tabs_s,opt)
cts1d_s = calc_cts1d(k,p,tabs_s,opt)

#=======#
# Plot  #
#=======#

# Plot params
cex	    = 2
cex_leg = 1.75
cex.legend = 1.75
lwd	    = 2.5
ncoarse = 10
k_coarse = coarse_grain(k,ncoarse)

# PDF
file = "~/Dropbox/18cts/plots/cts_h2o.pdf"
pdf(file,width=15,height=5,bg="white")
par(mfrow=c(1,3),mar=c(5,7,5,7))

# 2D fields
for (n_field in 1:length(fields)){
	field = fields[n_field]
	main  = mains2d[n_field]
	var2d = eval(as.name(paste(field,"2d_s",sep="")))
	var2d_coarse = coarse_grain(var2d,ncoarse)
	my.image.plot(1e-2*k_coarse,1e-2*p_s[pvec],fac2d*var2d_coarse[ ,pvec],
#		xlim 	 = , 
		ylim	 = plim,
		zlim     = Hklim,
		col		 = colvec,
		xlab 	 = klab,
		ylab     = plab,
		main     = main,
		cex.lab  = cex,
		cex.axis = cex,
		cex.main = cex,
		cex.legend = cex.legend)
}

# 1D fields
plot(1,type="n",
	xlim 	 = Hlim, 
	ylim	 = plim,
	xlab 	 = expression(H*","~H[CTS]~~"(K/day)"),	 
	ylab     = plab,
	main     = "(c)    Integrated heating",
	cex.lab  = cex,
	cex.axis = cex,
	cex.main = cex)
for (n in 1:length(fields)){
    field = fields[n]
    lty   = ltyvec[n]
    var   = eval(as.name(paste(field,"1d_s",sep="")))
    points(fac1d*var,1e-2*p_s,type="l",lwd=lwd,lty=lty)
}
legend("bottomright",legend=names1d,lty=ltyvec,lwd=lwd,cex=cex_leg)

dev.off()

