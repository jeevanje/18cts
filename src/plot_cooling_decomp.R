args = commandArgs(trailingOnly=TRUE)

library(fields)
library(ncdf4)
Rtoolsdir = "~/Dropbox/Rtools/"
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))
source(paste(Rtoolsdir,"gray_model.R",sep=""))

case   = args[1]
#case   = "h2o_only_no_cont"
D      = 1.5  # diffusion parameter
Nsec   = 3600*24   # sec/day
fac2d  = g/Cp*Nsec*1e2  # pppf to K/day/cm^-1
fac1d  = g/Cp*Nsec      # pppf to K/day
fields = c("coo","sx","ax","gex","cts","sum")
ltyvec = c("solid",rep("dashed",time=5))
colvec = c("black","green","orange","red","blue","black")

#=======#
# Data  #
#=======#

ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k")
dk      = k[2]-k[1]
z       = ncvar_get(nc,"z")
p       = ncvar_get(nc,"p")
tabs    = ncvar_get(nc,"tabs")
nz      = length(z)
nk      = length(k)
zvec    = 2:(nz-1)
coo2d   = ncvar_get(nc,"coo")/fac2d     # convert to SI 
opt     = D*ncvar_get(nc,"opt")
flx2d   = ncvar_get(nc,"flx")*1e-2      # W/m^2/m^-1

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

# 2d terms
cts2d_s  = pi*B_s*exp(-opt_s)*dtaudp2d  #s, z_s, W/m^2/Pa/m^-1, pppf units
gex2d_s  = -pi*(B_surf-B_s)*exp(-(opt_surf - opt_s))*dtaudp2d  #s, z_s
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
sum2d_s = cts2d_s + gex2d_s + ax2d_s + sx2d_s

for (field in fields){
    field2d = eval(as.name(paste(field,"2d_s",sep="")))
    assign(paste(field,"1d",sep=""), apply(field2d,2,sum)*dk)  # W/m^2/Pa
}

#=======#
# Plot  #
#=======#

# Plot params
coo_lim    = c(-0.005,0.01) # K/day/cm^-1
z_lim      = c(0,30)  # km
coo_int_lim= c(-0.25,2.5)   # K/day
cex	   	   = 2.5
ncoarse	   = 10
k_coarse   = coarse_grain(k,ncoarse)

# PDF
file = paste("~/Dropbox/18cts/plots/cooling_decomp_",case,".pdf",sep="")
pdf(file,width=18,height=10,bg="white")
par(mfrow=c(2,3),mar=c(5,5,5,11))
for (field in fields[1:5]){
    var = coarse_grain(eval(as.name(paste(field,"2d_s",sep=""))),ncoarse)
    my.image.plot(1e-2*k_coarse,1e-3*z_s,fac2d*var,
		ylim = z_lim,
		zlim = coo_lim,
		xlab = "k (cm^-1)",
		ylab = "z (km)",
		main = paste(field," (K/day/cm^-1)",sep=""),
		cex.lab  = cex,
	    cex.axis = cex,
		cex.main = cex,
		cex.legend = 2
		)
}
plot(1,type="n", xlim = coo_int_lim, ylim=z_lim,
	xlab 	 = "Cooling (K/day)",	 
	ylab     = "z (km)",
	cex.lab  = cex,
	cex.axis = cex,
	cex.main = cex)
abline(v=0,lty="dashed",col="gray")
for (n in 1:length(fields)){
    field = fields[n]
    lty   = ltyvec[n]
    col   = colvec[n]
    var   = eval(as.name(paste(field,"1d",sep="")))
    points(fac1d*var,1e-3*z_s,type="l",lwd=2,lty=lty,col=col)
}
legend("topright",legend=fields,lty=ltyvec,lwd=2,col=colvec,cex=2)
dev.off()
