# ====================================================================================
# Spatial Polygon Data Creation with MetaData at Adm level (locality/ward)
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
library(RColorBrewer)
source(paste0(sys, "Delhi/Codes/Source/overlay_func.R"))
# ====================================================================================

# Plan
# Create ward/locality level spatial polygon data set [flexibility of changing ward/locality]
# Meta Data
#     - Population
#     - Night Lights
#     - Green
#     - Building
#     - Barren (growth)
#     - Road Intersections
# Identify Platform -> Merge Vector Data -> Overlay -> Aggregate -> SPDF -> RasterBrickz


# Read Data
# ====================================================================================
# Locality/Ward Boundary
# -----------------------------
# # GRID --
# adm <- readShapeSpatial("Raw/Adm/Del_grid.shp")  # GRID
# # names(adm)[names(adm) == "CELLID"] <- "Sub-geography"
# # Locality --
# adm <- readShapeSpatial("Raw/Adm/Del_locality.shp")  # Locality
# # names(adm)[names(adm) == "LOCALITY"] <- "Sub-geography"
# # Ward --
# adm <- readShapeSpatial("Raw/Adm/Del_wards.shp")     # Wards
# names(adm)[names(adm) == "D_Name"] <- "Sub-geography"
# # PS Vornoi --
adm <- readShapeSpatial("Raw/Adm/Crpol_vrn.shp")  # Vornoi
names(adm)[names(adm) == "POI_NAME"] <- "Sub-geography"
dim(adm)
adm$LocID <- paste0(as.character("LocID_"),1:nrow(adm))

# Population
# -----------------------------
pop <- raster("Raw/WorldPop/Delhi/popmap15.tif")
plot(density(pop[[1]]))
pops <- rasterToPoints(pop, spatial=TRUE)
class(pops)

# Night Lights
# -----------------------------
nl <- raster("Raw/NightLights/Delhi/rad2000resam1_V1.tif")
plot(density(nl[[1]]))
nls <- rasterToPoints(nl, spatial=TRUE)
class(nls)

# Scan Pop
# -----------------------------
spop <- raster("Raw/LandScan/Delhi/LSpopu.tif")
plot(density(spop[[1]]))
spops <- rasterToPoints(spop, spatial=TRUE)
class(spop)

# VIIRS
# -----------------------------
vr <- raster("Raw/NightLights/Delhi/AllIndiaVIIRS_AvgAnn2014_15.tif")
plot(density(vr[[1]]))
vrs <- rasterToPoints(vr, spatial=TRUE)
class(vrs)
# ====================================================================================


# Raster to SPdf and Clip by Locality
# ====================================================================================
# ---------------------------------------------------------------------------------------
# Phase I: Integrating World Pop
# ---------------------------------------------------------------------------------------
# Overlay Population Data on Locality to get LOC ID
# -----------------------------
system.time(ovrly(spdf= pops, pol = adm))
loc.1 <- spdf.o
dim(loc.1)
#rm(spdf.o)

# Aggregate Data to Locality Level
# -----------------------------
names(loc.1)[names(loc.1) == "popmap15"] <- "WrldPp"
loc.d1 <- data.table(loc.1@data)
sd.cols <- "WrldPp"
loc.d2 <- loc.d1[,lapply(.SD, sum, na.rm =T), .SDcols = sd.cols, by = .(LocID)]  # simple
dim(loc.d2)

# Intersect Population Data on Locality
# -----------------------------
loc1 <- merge(x = adm, y = loc.d2, by = "LocID")
class(loc1)
rm(loc.1,loc.d1,loc.d2)

# ---------------------------------------------------------------------------------------
# Phase I: Integrating World Pop
# ---------------------------------------------------------------------------------------
# overlay Population Data on Locality to get LOC ID
# -----------------------------
system.time(ovrly(spdf= spops, pol = adm))
loc.1 <- spdf.o
dim(loc.1)
#rm(spdf.o)

# Aggregate Data to Locality Level
# -----------------------------
names(loc.1)[names(loc.1) == "LSpopu"] <- "LSpopu"
loc.d1 <- data.table(loc.1@data)
sd.cols <- "LSpopu"
loc.d2 <- loc.d1[,lapply(.SD, sum, na.rm =T), .SDcols = sd.cols, by = .(LocID)]  # simple
dim(loc.d2)

# Intersect Population Data on Locality
# -----------------------------
loc2 <- merge(x = loc1, y = loc.d2, by = "LocID")
class(loc2)
rm(loc.1,loc.d1,loc.d2)

# ---------------------------------------------------------------------------------------
# Phase III: Integrating VIIRS
# ---------------------------------------------------------------------------------------

# overlay Night Lights Data on Locality to get LOC ID
# -----------------------------
system.time(ovrly(spdf= vrs, pol = loc1))
loc.1 <- spdf.o
dim(loc.1)
rm(spdf.o)

# Aggregate Data to Locality Level
# -----------------------------
names(loc.1)[names(loc.1) == "AllIndiaVIIRS_AvgAnn2014_15"] <- "VIIRS"
loc.d1 <- data.table(loc.1@data)
sd.cols <- "VIIRS"
loc.d2 <- loc.d1[,lapply(.SD, sum, na.rm =T), .SDcols = sd.cols, by = .(LocID)]  # simple
dim(loc.d2)

# Intersect Population Data on Locality
# -----------------------------
loc3 <- merge(x = loc2, y = loc.d2, by = "LocID")
class(loc3)
spdf <- loc3[,c("LocID","LSpopu","WrldPp","VIIRS")]
names(spdf)

# Export Data into Shapefile
# ----------------------------------
# proj4string(spdf) <- CRS("+proj=longlat +datum=WGS84")
# writeOGR(spdf, "Processed/PolygonStacks","Popviirs_pol",
#          driver="ESRI Shapefile", overwrite = T)
# ====================================================================================



# Extension with Crime Attributes from NCRB@District Level
# ====================================================================================
# Import Delhi@District Data
# -----------------------------
crpdf <- readShapeSpatial("Raw/CrimeNCRB/Delhi/Police_del_all.shp")
# Total Crimes
crpdf$totCrm <- (crpdf$Murder+crpdf$Rape+crpdf$GangRap+ 
                 crpdf$Robbery+crpdf$Theft+crpdf$Asltmds+crpdf$Sxlhrsm) 
table(is.na(crpdf$totCrm))
# Total Police Peronnels Sanctioned
crpdf$tSnc <- (crpdf$InspSnc+crpdf$sI_SI_S+crpdf$aI_SI_S+
                    crpdf$HC_snc+crpdf$Cnst_sn+crpdf$Sxlhrsm) 
table(is.na(crpdf$tSnc))
# Total Police Peronnels Working
crpdf$tWrk <- (crpdf$Insp__M+crpdf$In__W_F+crpdf$sI_SI_M+crpdf$sI_SI_F+ 
                crpdf$aI_SI_M+crpdf$aI_SI_F+crpdf$HC_wrkM+crpdf$HC_wrkF+
                 crpdf$Cnst_wM+crpdf$Cnst_wF) 
table(is.na(crpdf$tWrk))
proj4string(crpdf) <- "+proj=longlat +datum=WGS84"
names(crpdf)

# Subset GIS Polygon Data for LocID
# -----------------------------
spd <- spdf[,"LocID"]
class(spd)

# Overlay Crime Delhi District Data on Prepared GIS Polygon Data to give it LOCID
# -----------------------------
system.time(ovrly(spdf= crpdf, pol = spd))
cr.spdf <- spdf.o
dim(cr.spdf)
rm(spdf.o)

# Keep only the IPC attributes data and Aggregate 
# -----------------------------
vr <- c("Totarea","Murder","Rape","GangRap","Robbery","Theft","Asltmds",
        "Sxlhrsm","InspSnc","Insp__M","In__W_F","sI_SI_S","sI_SI_M","sI_SI_F",
        "aI_SI_S","aI_SI_M","aI_SI_F","HC_snc","HC_wrkM","HC_wrkF","Cnst_sn",
        "Cnst_wM","Cnst_wF","Tl_shrt","totCrm","tSnc","tWrk") 
cr.d1 <- data.table(cr.spdf@data)
sd.cols <- vr
cr.d2 <- cr.d1[,lapply(.SD, sum, na.rm =T), .SDcols = sd.cols, by = .(LocID)]  # simple
dim(cr.d2)

# Intersect the District Delhi crime data on GIS Polygon 
# -----------------------------
crdf <- merge(x = spdf, y = cr.d2, by = "LocID")
dim(crdf)
names(crdf)

# Assign Dataset Same Name according to next codes
# -----------------------------
spdf <- crdf
table(is.na(crdf@data$totCrm))              # 75 wards with no PS
table(is.na(spdf@data$totCrm))
table(is.na(spdf@data))
spdf@data[is.na(spdf@data)] <- 0

# Export Data into Shapefile
# ----------------------------------
proj4string(spdf) <- CRS("+proj=longlat +datum=WGS84")
writeOGR(spdf, "Processed/PolygonStacks","Popviirs.cr_pol",
         driver="ESRI Shapefile", overwrite = T)
writeOGR(spdf, "Processed/PolygonStacks","Popviirs.cr_pol.wrds",
         driver="ESRI Shapefile", overwrite = T)
writeOGR(spdf, "Processed/PolygonStacks","Popviirs.cr_pol.vrn",
         driver="ESRI Shapefile", overwrite = T)
# ====================================================================================
# FIN