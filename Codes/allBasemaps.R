# All BaseMaps 
# -----------------------------------------------------
addProviderTiles('CartoDB.PositronNoLabels',group='Blank-Canvas') %>%
addProviderTiles('OpenStreetMap.BlackAndWhite', group="OSM-BlackNWhite") %>%  # Coolest
addProviderTiles('MapQuestOpen.OSM', group='MapQuest') %>%
addProviderTiles('Stamen.TonerLite', group='Stamen-Light') %>%
addProviderTiles('Esri.WorldStreetMap',group='Esri-1') %>%
addProviderTiles('Esri.DeLorme',group='Esri-2') %>%
addProviderTiles('Esri.OceanBasemap',group='Esri-3') %>%     # NA 
addProviderTiles('Esri.NatGeoWorldMap',group='NatGeo') %>%   # Grrenish 
addProviderTiles('CartoDB.Positron',group='CartoDB-1') %>%
addProviderTiles('CartoDB.PositronNoLabels',group='CartoDB-2') %>%
addProviderTiles('Stamen.TonerHybrid',group='CartoDB-2') %>%
addProviderTiles('Stamen.TonerLines',group='CartoDB-2') %>%
addProviderTiles('CartoDB.DarkMatter',group='CartoDB-3') %>%
addProviderTiles('CartoDB.DarkMatterNoLabels',group='CartoDB-4') %>%
addProviderTiles('Acetate.basemap',group='Acetate') %>%
addProviderTiles('Stamen.TonerLabels',group='Acetate') -> eventMap