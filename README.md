# Delhi-Crime-Map
## Delhi Crime Maps (R Shiny)
see interactive map (https://sociocartography.shinyapps.io/DelhiCrime/)


###Background:
Delhi's state of crime has translated into an unfortunate adage: Delhi IS a state of crime.
Information on history (data) and location of crime incident (spatial) helps in tracing/linking various crime events with community geography, which is broadly driven by income-demographics, rent rates, migration and city planning. Learning from data can help greatly in taking from curative to preventive measures at a government, society and individual level. 

###Analysis Vintage:
2014 was an year of paradigm shift for crimes in Delhi. Barely an year after to infamous Nirbhaya incident, this period was characterized with a pattern increase in urban crimes. In this looming backdrop, 2014 was the year when Delhi Police underwent a challenging transformation, by registering all cognizable offenses reported. Thus providing a germane ground for setting 2014 as a reference year for yarsdstick.  

###Data:
Open source data from WorldPop and Night Light radiance (NASA) was conjoined with Crime Statistics by police stations.
The full map is a culmination of integrating micro satellite data and integrating it with event based data. Goe-coding police station in Delhi the features are overlayed through 5 layers. There is additional information on number of police personnel deployed and sanctioned per police station so juxtaposing it with crime incident and population density can help in quick and efficient disbursal. Delhi and India by and large follows overlapping system of administrative boundaries, therefore I have stuck to creating customized Vornoi (Voronoi diagram is a partitioning of a plane into regions based on distance to points in a specific subset of the plane) boundaries. In essence the spatial and event based (converted into spatial data) was integrated together by a uniform spatial resolution. The meta data was then rasterized and condensed into a rasterBrick which is used by the slider panel in Shiny.

![crscr1](https://cloud.githubusercontent.com/assets/6264399/17649402/48e403d8-6252-11e6-9c74-6101d0337476.jpg)


###Crime Geography: Methodology Voronoi
Spatially, Delhi is divided by series of intersecting/overlapping boundaries. There are over 7 districts, 21 towns, 400+ wards overlapped with over 200+ localities, 800+ sub-localities, 7 parliamentary and 70 assembly boundaries. Under a decentralized governing structure each boundary has its own justified reason of existence but its discordance with other boundaries gravely affects the potential/possibility of realising this data for bigger analyses. 
Due to the existent multiplicity of administrative sub-geographies, it becomes an arduous exercise to reconcile data from multiple layers. This has been addressed in the SociocaRtography crime maps by first mapping out the 181 police stations and additional 99 police outposts. The co-ordinates of the police stations/outposts were collected by geocoding each and every police station address from various API sources (Google, Bing, Tor). Then following layers of information was added sequentially. Population from world population open sources data set, reporting population@1sqkm. Night-light radiance value collected from NASA. Police Personnel deployed and reported crime incidents by category from NGO’s.
None of the existing boundaries on wards, pin codes or localities follow a homogeneous pattern of either population, income or crime distribution. Concordantly the next step was to create usable/cogent/effective crime boundaries. This was done by using the police stations as central/nodal/axis points and customise boundaries by their spatial coverage across the city. This was achieved by application of a standard/known geo-spatial techniques of Voronoi polygons and specifying the algorithm for the current scenario. The result as can be seen is that size of polygon are dependent on the distance between two police stations. 
