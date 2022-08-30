###################Load required packages##################
library(rerddapXtracto)
library(raster)
library(lubridate)
library(grec)
library(rtsVis)
library("tidyverse")
library("gganimate")
library("gifski")
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
names(gradbrick)=idsst
#Animate data in a single step
animate(gradbrick, n=1,main=idsst)
#Animate data in a cooler way
# convert gridded raster data dataframe
g_df <- gradbrick %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", names(gradbrick))) %>%
  pivot_longer(cols = starts_with("X20"),
               names_to = "layer",
               values_to = "val")%>%
  mutate(layer = substr(layer, 2, 14)) %>%
  mutate(date = as.POSIXct(layer, "%Y.%m.%d", tz = "UTC")
  )

world_map <- ggplot() +
  theme_void() +
  geom_tile(
    data = g_df,
    aes(
      x = x,
      y = y,
      fill = val,
      group = date)) +
  scale_fill_viridis_c(
    option = "B"
  ) +labs(x = "Longitude", y = "Latitude") +
  theme(
    plot.title = element_text(
      family = "Prata",
      size = 20,
      hjust = 0.5),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 20),
    plot.caption = element_text(
      color = "grey50",
      size = 14,
      hjust = 0.9),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    title = "Thermal gradients °C",
    subtitle = "Northern Chilean Patagonia",
    caption = "SST Data source: ERDDAP at {current_frame}")+
  transition_manual(date)

gganimate::animate(
  world_map,
  width = 400,
  height = 600,
  renderer=gifski_renderer("map.gif")
)
