# ====================================================================================
# Raster Brick from Spatial Polygon Data
# ====================================================================================

# Prilimnaries
# ====================================================================================
rm(list=ls())
gc(verbose=T)
gcinfo(FALSE)
sys <- "C:/Parth/Personal/Data Mining/sociocaRtography/"
setwd(paste0(sys, "Delhi/Data"))
# ====================================================================================

# Libraries
# ====================================================================================
library(maps)
library(mapproj)
library(sp)
library(maptools)
library(raster)
library(leaflet)
library(shiny)
library(ggplot2)
library(shinyBS)
library(rgdal)
library(dplyr)
# ====================================================================================


# Read Files
# ====================================================================================
# Locality/Ward Boundary
# -----------------------------
# adm <- readShapeSpatial("Raw/Adm/Del_locality.shp")
adm <- readShapeSpatial("Raw/Adm/Del_wards.shp")     # Wards
dim(adm)

# Import Crime Data on District Polygons
# ---------------------------
# spdf <- readShapeSpatial("Processed/PolygonStacks/Popviirs_pol.shp")    # POP@loc
# spdf <- readShapeSpatial("Processed/PolygonStacks/Popviirs.cr_pol.shp") # POP+Crime@loc
# spdf <- readShapeSpatial("Processed/PolygonStacks/Popviirs.cr_pol.wrdar.shp")   # POP+Crime@wrd 
spdf <- readShapeSpatial("Processed/PolygonStacks/Popviirs.cr_pol.vrn.shp")   # POP+Crime@wrd 
class(spdf)
names(spdf)

# Create Area-Normalised attributes
# ---------------------------
summary(spdf@data$LSpopu)
summary(spdf@data$Area_ward)
summary(spdf@data$LSpopu/spdf@data$Area_ward)

# Retain the original data
# ---------------------------
spdf.1 <- spdf
summary(spdf@data$LSpopu)
spdf@data$LSpopu <- spdf@data$LSpopu/spdf@data$Area_ward
summary(spdf@data$LSpopu)

summary(spdf@data$WrldPp)
spdf@data$WrldPp <- spdf@data$WrldPp/spdf@data$Area_ward
summary(spdf@data$WrldPp)

summary(spdf@data$VIIRS)
spdf@data$VIIRS <- spdf@data$VIIRS/spdf@data$Area_ward
summary(spdf@data$VIIRS)

summary(spdf@data$totCrm)
spdf@data$totCrm <- spdf@data$totCrm/spdf@data$Area_ward
summary(spdf@data$totCrm)

summary(spdf@data$tWrk)
spdf@data$tWrk <- spdf@data$tWrk/spdf@data$Area_ward
summary(spdf@data$tWrk)

summary(spdf@data$tSnc)
spdf@data$tSnc <- spdf@data$tSnc/spdf@data$Area_ward
summary(spdf@data$tSnc)
# ====================================================================================


# GRID Calculations
# ====================================================================================
# Import Night Light tiff
# ---------------------------
ras <- raster("Raw/NightLights/Delhi/rad2000resam1_V1.tif")
class(ras)
extent(ras)

# Extent Computations
#---------------------------
x_min <- extent(ras)[1]
x_max <- extent(ras)[2]
y_min <- extent(ras)[3]
y_max <- extent(ras)[4]
cell_res <- 0.08333333  # 10 km GRID
cell_res <- 0.008333333  # 10 km GRID
#---------------------------

# Number of Columns & Rows : Computations
#--------------------------
x_extent <-as.integer((x_max-x_min)/cell_res)
x_extent.1 <- (((x_max-x_min)/cell_res)-as.integer((x_max-x_min)/cell_res))  #Long
y_extent <- as.integer((y_max-y_min)/cell_res)
y_extent.1 <- (((y_max-y_min)/cell_res)-as.integer((y_max-y_min)/cell_res))  #Lat
n_row <- ifelse(y_extent.1>0.5,(y_extent+2),(y_extent+1))    #lat
n_col <- ifelse(x_extent.1>0.5,(x_extent+2),(x_extent+1))    #long

# Empty Raster
#--------------------------
ras1 <- raster(nrow=n_row,ncol=n_col)
extent(ras1) <- extent(ras)

# Resampling from 1km to 10km
#--------------------------
ras2 <- resample(ras,ras1,method='bilinear')
extent(ras2) <- extent(spdf)
# ====================================================================================



# Rasterize: Raster Layers
# ====================================================================================
pth <- paste0(sys, "Delhi/Data/Processed/Rasters/Popviirs")

# Rasterize: Wrldpp
#--------------------------
rbwP <- rasterize(spdf, ras2, "WrldPp", update=TRUE, CRS("+init=epsg:4326"),fun=sum)
class(rbwP)
# Write
# writeRaster(rbwP,file.path(pth,"stWrldPp"),format="GTiff",overwrite=TRUE)
# Cropping the IMagery for NA's
e=extent(adm)
subsetwP <- crop(rbwP, e, snap="out")
frwP <- rasterize(adm, subsetwP)   
lrwP <- mask(x=subsetwP, mask=frwP) 
summary(lrwP@data@values)

# Rasterize: LSpopu
#--------------------------
rblP <- rasterize(spdf, ras2, "LSpopu", update=TRUE, CRS("+init=epsg:4326"),fun=sum)
class(rblP)
# Write
# writeRaster(rblP,file.path(pth,"stLSpopu"),format="GTiff",overwrite=TRUE)
# Cropping the IMagery for NA's
subsetlP <- crop(rblP, e, snap="out")
frlP <- rasterize(adm, subsetlP)   
lrlP <- mask(x=subsetlP, mask=frlP) 
summary(lrlP@data@values)

# Rasterize: VIIRS
#--------------------------
rbnL <- rasterize(spdf, ras2, "VIIRS", update=TRUE, CRS("+init=epsg:4326"),fun=sum)
class(rbnL)
# Write
# writeRaster(rb1,file.path(pth,"stVIIRS"),format="GTiff",overwrite=TRUE)
# Cropping the IMagery for NA's
subsetnL <- crop(rbnL, e, snap="out")
frnL <- rasterize(adm, subsetnL)   
lrnL <- mask(x=subsetnL, mask=frnL) 
summary(lrnL@data@values)

# Rasterize: Delhi Total PS Crime
#--------------------------
rbtC <- rasterize(spdf, ras2, "totCrm", update=TRUE, CRS("+init=epsg:4326"),fun=sum)
class(rbtC)
# Write
# writeRaster(rb1,file.path(pth,"stIPC13"),format="GTiff",overwrite=TRUE)
# Cropping the IMagery for NA's
subsettC <- crop(rbtC, e, snap="out")
frtC <- rasterize(adm, subsettC)   
lrtC <- mask(x=subsettC, mask=frtC) 
summary(lrtC@data@values)

# Rasterize: Delhi Total PS Crime
#--------------------------
rbtS <- rasterize(spdf, ras2, "tSnc", update=TRUE, CRS("+init=epsg:4326"),fun=sum)
class(rbtS)
# Write
# writeRaster(rb1,file.path(pth,"stIPC13"),format="GTiff",overwrite=TRUE)
# Cropping the IMagery for NA's
subsettS <- crop(rbtS, e, snap="out")
frtS <- rasterize(adm, subsettS)   
lrtS <- mask(x=subsettS, mask=frtS) 
summary(lrtS@data@values)

# Rasterize: Delhi Total PS Crime
#--------------------------
rbtW <- rasterize(spdf, ras2, "tWrk", update=TRUE, CRS("+init=epsg:4326"),fun=sum)
class(rbtW)
# Write
# writeRaster(rb1,file.path(pth,"stIPC13"),format="GTiff",overwrite=TRUE)
# Cropping the IMagery for NA's
subsettW <- crop(rbtW, e, snap="out")
frtW <- rasterize(adm, subsettW)   
lrtW <- mask(x=subsettW, mask=frtW) 
summary(lrtW@data@values)

#-------------------------- Getting 1sq grid values as it is --------------------------#

# Scan Pop
# -----------------------------
spop <- raster("Raw/LandScan/Delhi/LSpopu.tif")
class(spop)
# Resolution Calibration
sp <- resample(spop, lrtC, method = "ngb")
# Cropping the IMagery for NA's
subsetsp <- crop(sp, e, snap="out")
frsp <- rasterize(adm, subsetsp)   
lrsp <- mask(x=subsetsp, mask=frsp) 
summary(lrsp@data@values)


# VIIRS
# -----------------------------
vr <- raster("Raw/NightLights/Delhi/AllIndiaVIIRS_AvgAnn2014_15.tif")
class(vr)
# Resolution Calibration
vras <- resample(vr, lrtC, method = "ngb")
# Cropping the IMagery for NA's
subsetvr <- crop(vras, e, snap="out")
frvr <- rasterize(adm, subsetvr)   
lrvr <- mask(x=subsetvr, mask=frvr) 
summary(lrvr@data@values)

plot(vras)
plot(lrvr)
# ====================================================================================


# Brick
# ====================================================================================
# Stack and Brick
#--------------------------
stck <- stack(lrwP,lrlP,lrnL,lrtC,lrtS,lrtW)
stck <- stack(lrtC,lrlP,lrtW,lrtS,lrtnL)   # only lspop and VIIRS & Del Dst Crm
stck <- stack(lrtC,lrsp,lrtW,lrtS,lrvr)   # only lspop and VIIRS & Del Dst Crm
class(stck)
br.pc <- brick(stck)
nlayers(br.pc)
x <- br.pc
class(x)

# # Variable for looping: Export individual raster to stack them
# # --------------------------------
# library(Hmisc)
# v2 <- v1[v1 %nin% c("IPC_2001")]
# t1 <- Sys.time()
# for(i in v2 )
# {
#   print(i)
#   rar <- rasterize(ipc, ras2,i, update=TRUE, CRS("+init=epsg:4326"),fun=sum)
#   print(summary(rar))
#   flnm <- paste0("st", i)
#   print(flnm)
#   writeRaster(rar,file.path(sim,flnm),format="GTiff",overwrite=TRUE)
#   r4 <<- stack(r4,rar)
#   cat("\n")
# }
# t2 <- Sys.time()- t1
# ====================================================================================

# Final Objects
# ====================================================================================
# Data I (locs)
# -----------------------------
# dt <- readShapeSpatial("Raw/Commercial/Delhi_police_stations.shp")
dt <- readShapeSpatial("Raw/CrimeNCRB/Delhi/Police_del_all.shp")
locs <- as.data.frame(dt)
locs <- locs[,c("x","y","Pol.Stn")]
locx <- locs
names(locx)[names(locx) == "Pol.Stn"] <- "loc"
names(locx)[names(locx) == "y"] <- "lat"
names(locx)[names(locx) == "x"] <- "lon"
locx$loc <- as.character(locx$loc)
locs <- locx

# Data Primers
# -----------------------------
decades <- seq(2001, 2005, by=1)
layers <- seq(1, 5, by=1)
# layers <- c("Pop", "VIIRS")
lon <- 77.1500
lat <- 28.6139
d <- d.cru$Locs[[2]] %>% filter(Month=="Jan")
# ====================================================================================

# Saving Image
# ====================================================================================
# # PopVIIRS@locality
# save(x,locs,lat,lon,decades,layers, 
#      file = "C:/Parth/Personal/Data Mining/SociocaRtography/Delhi/Sessions/DelhiCrime/PpVr_v1.RData")

# # POPVIIRS Crime@locality
# save(x,locs,lat,lon,decades,layers, 
#      file = "C:/Parth/Personal/Data Mining/SociocaRtography/Delhi/Sessions/DelhiCrime/PpVrCr_v1.RData")

# # POPVIIRS Crime@wards
# save(x,locs,lat,lon,decades,layers, 
#      file = "C:/Parth/Personal/Data Mining/SociocaRtography/Delhi/Sessions/DelhiCrime/PpVrCrwrd_v1.RData")

# POPVIIRS Crime@wards: Ward Area Normalised
# save(x,locs,lat,lon,decades,layers, 
#      file = "C:/Parth/Personal/Data Mining/SociocaRtography/Delhi/Sessions/DelhiCrime/PpVrCrwrd_v2.RData")

# POPVIIRS Crime@vornoi: Vornoi 161 PS [27th July, 2016]
# save(x,locs,lat,lon,decades,layers, 
#      file = "C:/Parth/Personal/Data Mining/SociocaRtography/Delhi/Sessions/DelhiCrime/CrPpTwvrn_v1.RData")

# POPVIIRS Crime@vornoi: Vornoi 161 PS + Original LandScan + VIIRS [30th July, 2016]
# save(x,locs,lat,lon,decades,layers, 
#      file = "C:/Parth/Personal/Data Mining/SociocaRtography/Delhi/Sessions/DelhiCrime/CrPpTwvrn_v2.RData")
# ====================================================================================
# FIN