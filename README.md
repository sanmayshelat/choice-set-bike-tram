# README #

Elimination-by-aspects based choice set generation method from [1]. Implemented here to identify access stations considered by travellers in [2].

If you use this, please cite either: 

[1] Shelat, S., Cats, O., van Oort, N., & van Lint, H. (2019, June). Calibrating route choice sets for an urban public transport network using smart card data. In 2019 6th International Conference on Models and Technologies for Intelligent Transportation Systems (MT-ITS) (pp. 1-8). IEEE.

[2] Ton, D., Shelat, S., Nijënstein, S., Rijsman, L., van Oort, N., Hoogendoorn, S., 2020. Understanding the Role of Cycling to Urban Transit Stations through a Simultaneous Access Mode and Station Choice Model. Transp. Res. Record.

### Data Required ###

* GTFS data
* This study used on-board survey data of tram users (Rijsman, L. et al, 2019). This contained information on: 
	* origin stop, 
	* destination stop, 
	* tram line, 
	* access mode, 
	* home location.
* This survey data was adjusted to show per user the access stations they had as options (a maximum access threshold was set, see [2])


### How to: Run in order ###

1. gtfsRouteExtract: Extracting routes and route attributes from GTFS data
2. topoCsgmCombined: Generating choice set based on topological rules
3. eliminationByAspects: Apply the EBA algorithm for choice set generation


### Abstract [1] ###

Identifying the set of alternatives from which travellers choose their routes is a crucial step in estimation and application of route choice models. These models are necessary for the prediction of network flows that are vital for the planning of public transport networks. However, choice set identification is typically difficult because while selected routes are observed, those considered are not. Approaches proposed in literature are not completely satisfactory, either lacking transferability across networks (observation-driven methods) or requiring strong assumptions regarding traveller behaviour (uncalibrated choice set generation methodologies (CSGM)). Therefore, this study proposes a constrained enumeration CSGM that applies the non-compensatory decision model, elimination-by-aspects, for choice set formation. Subjective assumptions of traveller preferences are avoided by calibrating the decision model using observed route choice behaviour from smart card data, which is becoming increasingly available in public transport systems around the world. The calibration procedure also returns two key insights regarding choice set formation behaviour: (i) the ranking of different attributes by their importance, and (ii) the acceptable detours for each attribute. To demonstrate the methodology and investigate choice set formation behaviour, the tram and bus networks of The Hague, Netherlands are used as a case study.

### Abstract [2] ###

Governments worldwide are aiming for an increase in sustainable mode use to increase sustainability, livability and accessibility. Integration of bicycle and transit can increase catchment areas of transit compared to walking and thus provide better competition to non-sustainable modes. To achieve this, effective measures have to be designed that require a better understanding of the factors influencing access mode and station choice. At the national/regional level this has been thoroughly studied. However, at the urban level knowledge is missing. This study aims to investigate which factors influence the joint decision for tram access mode and tram station choice. The joint investigation can identify trade-offs between the access and transit journey. Furthermore, the effect of each factor on the bicycle catchment area is investigated. Using data from tram travelers in The Hague, Netherlands, a joint simultaneous discrete choice model is estimated. Generally, walking is preferred over cycling. Our findings suggest that access distance is one of the main factors for explaining the choice, where walking distance is weighted 2.1 times cycling distance. Frequent cyclists are more likely to also cycle to the tram station, whereas frequent tram users are less inclined to do so. Bicycle parking facilities increase the cycling catchment area by 234m. The transit journey time has the largest impact on the catchment area of cyclists. Improvements to the system, such as less stops and/or higher frequency (like LRT) result with a much higher accepted cycling distance.
