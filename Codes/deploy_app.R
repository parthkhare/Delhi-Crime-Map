library(rsconnect)
rsconnect::setAccountInfo(name='sociocartography',
                          token='26BB882140906240DE4D5B8088D15453',
                          secret='pY0rs6m97iBaeUoc8SdlWy89bdhXgjeyKBYHjzDj')

sys <- "C:/Parth/Personal/Data Mining/sociocaRtography/"
setwd(paste0(sys, "Delhi/Sessions/DelhiCrime"))
rm(list=ls())
gc(verbose=T)
deployApp(account='sociocartography')
