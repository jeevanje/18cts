library(fields)
library(ncdf4)

Rtoolsdir = "~/Dropbox/Rtools/"
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))

      	  
#======#
# Data #
#======#

kvals  = 1e2*(1:1500)  # m^-1
Tvals  = 200:300
Nk     = length(kvals)
NT     = length(Tvals)
planck = outer(kvals,Tvals,planck_k)
alpha  = (log(planck[ ,2:NT])-log(planck[ ,1:(NT-1)]))/(rep(1,times=Nk)%o%diff(log(Tvals)))
Tvec   = 1:(NT-1)


#======#
# Plot #
#======#

# Plot params
alpha_lim  = c(0,10)
tabs_lim   = c(300,200) # K
cex        = 1.5

pdf("~/Dropbox/18cts/plots/alpha.pdf",width=7,height=5,bg="white")
par(mar=c(5,5,5,7))
contour(1e-2*kvals,Tvals[Tvec],alpha,
	     	 ylim = tabs_lim,
	   		 #zlim = alpha_lim,
             xlab = klab,
             ylab = tabslab,
             main = bquote(alpha==frac(partialdiff*ln~B(k,T),partialdiff*ln~T)),
			 labcex = 1.25,
			 lwd    = 1.5,
             cex.lab  = cex,
             cex.axis = cex,
             #cex.legend = cex
             cex.main = cex)
dev.off()

