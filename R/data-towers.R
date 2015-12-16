library(RgoogleMaps)

library(matrixStats)

library(ncf)


geo <- read.csv("SITE_ARR_LONLAT.CSV")

map <- GetMap.bbox(lonR = c(-11, -18), latR = c(12, 17), destfile = "map.png",maptype="terrain")

PlotOnStaticMap(map,lat = geo$lat, lon = geo$lon,zoom=10,pch=19,col="red",FUN=points, add=F)

map2 <- GetMap.bbox(lonR = c(-17.5, -17.35), latR = c(14.6, 14.9), destfile = "map2.png",maptype="terrain",zoom=12)

PlotOnStaticMap(map2,lat = geo$lat, lon = geo$lon,zoom=10,cex=0.5,pch=1,col="red",FUN=points, add=F)

data <- readLines("tsVis_per_25.csv")

data <- read.csv("tsVis_per_25.csv", header=F)

mat <- data.matrix(data)

cl <- correlog(geo$lon, geo$lat, mat, latlon=T, increment=2, resamp=5)

n=1666

cl <- correlog.nc(geo$lon, geo$lat, mat, latlon=T, increment=2, resamp=5)

s <- mat[1:396,1200]/sum(mat[1:396,1200])

ls <- loess(s~c(1:length(s)), span=0.15)

pr.loess <- predict(ls)

#plot(s,type="o")
lines(pr.loess~c(1:length(s)),col="grey")

rho <- matrix(0,nrow=8760,ncol=396)

for (i in c(1:8760)) {

    ssum <- sum(mat[1:396,i])
    
    if (ssum == 0) { ssum <- 1 }
    
    s <- mat[1:396,i]/ssum

    ls <- loess(s~c(1:length(s)), span=0.15)

    pr.loess <- predict(ls)
    
    rho[i,] <- pr.loess
    
    }
    
# plot daily dynamics
# 6 << 1st sunday

day = 6+7*35+4
start = (day-1)*24+12
end = start+6

mx <- max(rho[start:end,])

daymean <- colMeans(rho[start:end,])
daysd <- colSds(rho[start:end,])

plot(rho[start,], type="l", col="grey", xlim=c(0,400), ylim=c(0,mx+0.2*mx))


for (i in c( (start+1):end)) {

    lines(rho[i,],col="grey")
    
    }
    
lines(daymean,col="black", lwd=1)
lines(daymean+daysd,col="red", lwd=2)
lines(daymean-daysd,col="red", lwd=2)

# daily means

dm <- matrix(0,nrow=365,ncol=396)
dsd <- matrix(0,nrow=365,ncol=396)
asd <- c()

for (day in c(1:365)) {

    start = (day-1)*24
    end = start+24
    dm[day,] <- colMeans(rho[start:end,])
    dsd[day,] <- colSds(rho[start:end,])
    asd <- c(asd,sum(dsd[day,]))
    
    }
    
# plot daily means

mx <- max(dm)

plot(dm[1,], type="l", col="grey", xlim=c(0,400), ylim=c(0,mx+0.2*mx))

for (i in c( 2:365)) {

    lines(dm[i,],col="grey")
    
    }
    
# plot daily sd

mx <- max(dsd[,1:300])

plot(dsd[6,1:300], type="l", col="grey", ylim=c(0,mx+0.1*mx))

for (i in c(seq(6,365,7),seq(5,365,7))) {

    lines(dsd[i,1:300],col="blue")
    
    }

#mx <- max(dsd)

#plot(dsd[6,200:280], type="l", col="grey", ylim=c(0,mx+0.1*mx))

for (i in setdiff(setdiff(c(1:365),seq(6,365,7)),seq(5,365,7)) ) {

    lines(dsd[i,1:300],col="red")
    
    }



    
yearmean <- colMeans(dm)

lines(yearmean,col="black", lwd=2)

# sd for each day

plot(asd,type="l", col="grey")


# dynamic plot of signal diff

for (i in c(1:1000)) {

    name = paste("./plots/p_",formatC(i, format="d", digits=3, flag="0"),".png", sep="")

    png(file=name, bg="transparent")

    PlotOnStaticMap(map2,lat = geo$lat[1:390], lon = geo$lon[1:390],zoom=10,cex=0.5+10000*(rho[(300+i),1:390]-rho[(300+i-1),1:390]),pch=1,col="red",FUN=points, add=F)

    dev.off()

}

require(graphics)

cols <- gray.colors(150)


for (i in c(1:1000)) {

    gc <- round(100*abs(rho[(300+i),1:390]-rho[(300+i-1),1:390])/max(abs(rho[(300+i),1:390]-rho[(300+i-1),1:390])))

    name = paste("./plots/p_",formatC(i, format="d", digits=3, flag="0"),".png", sep="")

    png(file=name, bg="transparent")

    PlotOnStaticMap(map2,lat = geo$lat[1:390], lon = geo$lon[1:390],zoom=10,cex=1.5,pch=16,col=cols[gc],FUN=points, add=F)

    dev.off()

}



# daily means per 6 hours, working days

# daily means

dm <- matrix(0,nrow=365,ncol=396)
dsd <- matrix(0,nrow=365,ncol=396)
asd <- c()
idm <- matrix(0,nrow=23,ncol=396)

for (i in c(0:22)) {

for (day in setdiff(setdiff(c(1:365),seq(6,365,7)),seq(5,365,7))) {

    start = (day-1)*24+i
    end = start+2
    dm[day,] <- colMeans(rho[start:end,])
    dsd[day,] <- colSds(rho[start:end,])
    asd <- c(asd,sum(dsd[day,]))
    
    }
    
idx <- setdiff(setdiff(c(1:365),seq(6,365,7)),seq(5,365,7))

idmean <- colMeans(dm[idx,])

idm[i+1,] <- idmean

}

for (i in c(1:23)) {

    name = paste("./plots/p_",formatC(i, format="d", digits=3, flag="0"),".png", sep="")
    
    png(file=name)

    #plot(idm[i,],type="l", ylim=c(0,0.004))
    
    plot(idm[i,],type="o", ylim=c(0,0.004),yaxt="n",frame.plot=FALSE, ylab="",pch=1, main=paste("mean working day hour: ",formatC(i, format="d", digits=1, flag="0")),cex=0.5,col="black")
    
    dev.off()

}
    
# plot daily means

mx <- max(dm)

plot(dm[1,], type="l", col="grey", xlim=c(0,400), ylim=c(0,mx+0.2*mx))

for (i in setdiff(setdiff(c(1:365),seq(6,365,7)),seq(5,365,7))) {

    lines(dm[i,],col="grey")
    
    }

lines(idmean,col="black")


# daily means per 6 hours, weekends

# daily means

pdf(file="weekend_daily_means_18-24.pdf")

dm <- matrix(0,nrow=365,ncol=396)
dsd <- matrix(0,nrow=365,ncol=396)
asd <- c()

for (day in c(seq(6,365,7),seq(5,365,7))) {

    start = (day-1)*24
    end = start+6
    dm[day,] <- colMeans(rho[start:end,])
    dsd[day,] <- colSds(rho[start:end,])
    asd <- c(asd,sum(dsd[day,]))
    
    }
    
# plot daily means

mx <- max(dm)

plot(dm[1,], type="l", col="grey", xlim=c(0,400), ylim=c(0,mx+0.2*mx))

for (i in c(seq(6,365,7),seq(5,365,7))) {

    lines(dm[i,],col="grey")
    
    }
    
idx <- c(seq(6,365,7),seq(5,365,7))

idmean <- colMeans(dm[idx,])

lines(idmean,col="black")
    
dev.off()





