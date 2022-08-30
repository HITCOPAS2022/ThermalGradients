###################Load required packages##################
library(rerddapXtracto)
library(raster)
library(lubridate)
library(grec)
#Set spatial extent for extracting data
ext=extent(-79,-68,-47,-39)
#Use depth data to crete a mask
depthInfo <- rerddap::info("etopo360")#Bathymetry database
depth <- rerddapXtracto::rxtracto_3D(depthInfo,parameter ="altitude" ,xcoord = c(-79,-68), ycoord=c(-47,-39),tcoord = NULL,zcoord = NULL)
mask=raster(depth$depth)
mask=t(flip(mask,1))
mask[mask>=0]=NA
mask[mask<0]=1
mask=setExtent(mask,ext)
plot(mask)
#########Generate sst and thermal gradient rasters #########################################
sstInfo <- rerddap::info('jplMURSST41')#Daily SST database
sst <- rerddapXtracto::rxtracto_3D(sstInfo,parameter ="analysed_sst" ,xcoord = c(-79,-68), ycoord=c(-47,-39), tcoord=c(floor_date(as.POSIXct("2022-01-01 15:00:00"),unit="day"),as.POSIXct("2022-01-31 15:00:00")))
datasst=sst$analysed_sst
idsst=seq(as.Date(sst$time[1]),as.Date(sst$time[length(sst$time)]),1)

rastersst=list()
for(i in 1:length(idsst)){
  rastersst[[i]]=raster(datasst[,,i])
  rastersst[[i]]=t(flip(rastersst[[i]],1))
}
gradientsst=list()
for(i in 1:length(idsst)){
  gradientsst[[i]]=grec::detectFronts(as.matrix(rastersst[[i]]), method = "BelkinOReilly2009")
}
Grastersst=list()
for(i in 1:length(idsst)){
  Grastersst[[i]]=raster(gradientsst[[i]])
  crs(Grastersst[[i]]) ="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  Grastersst[[i]]=setExtent(Grastersst[[i]],ext)
  crs(rastersst[[i]]) ="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  rastersst[[i]]=setExtent(rastersst[[i]],ext)
}
#Use coordinates in new gradient rasters to update mask into the same resolution and extent
coord <- xyFromCell(Grastersst[[1]],1:ncell(Grastersst[[1]]))
mask=extract(mask,coord)
mask=raster(matrix(mask,ncol =ncol(Grastersst[[1]]),nrow = nrow(Grastersst[[1]]),byrow = T ))
mask=setExtent(mask,ext)
plot(Grastersst[[1]]*mask)
gradbrick=brick(unlist(Grastersst))
gradbrick=gradbrick*mask
#Animate data
animate(gradbrick, n=1,main=idsst)
