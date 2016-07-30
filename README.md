# Delhi-Crime-Map
Crime map of Delhi using R

Background
Delhi’s state of crime has translated into an unfortunate adage: Delhi is state of crime.
We learn from history and data that there is iterated pattern signals us deconstructing causes and so move from curative to preventive action.

Brief description of the submission:
Open source data from WorldPop and Night Light radiance (NASA) was conjoined with Crime Statistics by police stations.
The full map is a culmination of integrating micro satellite data and integrating it with event based data. Goe-coding police station in Delhi the features are overlayed through 5 layers. There is additional information on number of police personnel deployed and sanctioned per police station so juxtaposing it with crime incident and population density can help in quick and efficient disbursal. Delhi and India by and large follows overlapping system of administrative boundaries, therefore I have stuck to creating customized Vornoi (Voronoi diagram is a partitioning of a plane into regions based on distance to points in a specific subset of the plane) boundaries. In essence the spatial and event based (converted into spatial data) was integrated together by a uniform spatial resolution. The meta data was then rasterized and condensed into a rasterBrick which is used by the slider panel in Shiny.
