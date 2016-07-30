# ====================================================================================
# Public Sourced Crime Data Preparation
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


# Standaridize Point Crime Shapefile
# ====================================================================================
# Shapefile with Incident (Serene Prepared)
# ---------------------------
cr.in <- readShapeSpatial("Raw/CrimeNCRB/Delhi/Police_del_all.shp")
names(cr.in)

# Aggregate Statistics
# ---------------------------
cr.in$totCrm <- (cr.in$Murder+cr.in$Rape+cr.in$GangRap+ 
                   cr.in$Robbery+cr.in$Theft+cr.in$Asltmds+cr.in$Sxlhrsm) 
table(is.na(cr.in$totCrm))
# Total Police Peronnels Sanctioned
cr.in$tSnc <- (cr.in$InspSnc+cr.in$sI_SI_S+cr.in$aI_SI_S+
                 cr.in$HC_snc+cr.in$Cnst_sn+cr.in$Sxlhrsm) 
table(is.na(cr.in$tSnc))
# Total Police Peronnels Working
cr.in$tWrk <- (cr.in$Insp__M+cr.in$In__W_F+cr.in$sI_SI_M+cr.in$sI_SI_F+ 
                 cr.in$aI_SI_M+cr.in$aI_SI_F+cr.in$HC_wrkM+cr.in$HC_wrkF+
                 cr.in$Cnst_wM+cr.in$Cnst_wF) 
table(is.na(cr.in$tWrk))

# Export as CSV with lat longs
# ---------------------------
dt <- cr.in@data
dt$long <- cr.in@coords[,1]
dt$lat <- cr.in@coords[,2]
write.csv(dt, file = "Raw/CrimeNCRB/Delhi/PSDel.csv")


# Csv with Personnel Deployed (Serene Prepared)
# ---------------------------
cr.p <- read.csv("Raw/CrimeNCRB/Delhi/Delhi_PSPersonnel.csv")
ras <- raster("Raw/NightLights/Delhi/rad2000resam1_V1.tif")
class(ras)
extent(ras)
# ====================================================================================

# Standardize Vornoi shapefile with location of police stations
# ====================================================================================
# Shapefile with Incident (Serene Prepared)
# ---------------------------
cr.v <- readShapeSpatial("Processed/PolygonStacks/Popviirs.cr_pol.vrn.shp")
names(cr.v)

# Subset Polygon data frame
# -----------------------------
cr.v1 <- cr.v[,"LocID"]
class(cr.v1)

# Overlay Point Crime Data on Crime Vornoi to get LOC ID
# -----------------------------
system.time(ovrly(spdf= cr.in, pol = cr.v1))
ptvrn <- spdf.o
dim(ptvrn)
length(unique(ptvrn$LocID))                       # 280/297
rm(spdf.o)

# Extract lat-long data form Point Shapefile
# -----------------------------
ptvrn$lat <- ptvrn@coords[,2]
ptvrn$long <- ptvrn@coords[,1]
ptvr <- ptvrn[,c("LocID", "lat", "long")]
names(ptvr)

# Intersect with point data to get information on lat-long for georef in cartodb
# ---------------------------
cr.vrn <- merge(x = cr.v, y = ptvr@data, by = "LocID", all.x =T)
dim(cr.vrn)
class(cr.vrn)

# Summary Checks
# ---------------------------
plot(density(cr.vrn@data$LSpopu))
length(unique(cr.vrn@data$LocID))

# Export for cartodb
# ---------------------------
writeOGR(cr.vrn, "Processed/PolygonStacks","Dl.cr_pol.vrn",
         driver="ESRI Shapefile", overwrite = T)
# ====================================================================================
# FIN

