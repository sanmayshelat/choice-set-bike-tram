#%% Author and Title Information
# Title: Topological route choice generation: Convert to NumPy
# Author: Sanmay Shelat
# Institute: Dept. of Transport and Planning, Delft University of Technology

#%% Imports
import pandas as pd
import numpy as np
import myGlobalFxns as gf
import pickle
kjhkjh
#%% Import GTFS
def importGtfsRouteStops(day):
    # Read in files
    datesFile = pd.read_csv('../gtfs/calendar_dates.txt')
    tripsFile = pd.read_csv('../gtfs/trips.txt', 
                            usecols = ['route_id','service_id','trip_id','direction_id','trip_headsign'])
    stopTimesFile = pd.read_csv('../gtfs/stop_times.txt',
                                usecols = ['trip_id','stop_sequence','stop_id','arrival_time','timepoint'])
    stopsFile = pd.read_csv('../gtfs/stops.txt',usecols = ['stop_id','stop_name','parent_station','stop_lat','stop_lon'])
    routesFile = pd.read_csv('../gtfs/routes.txt')#,usecols = ['route_id'])
    
    # Select network (here only tram and bus; no historic tram)
    cond = (routesFile['agency_id']=='HTM')|(routesFile['agency_id']=='HTMBUZZ')
    selectedRoutes = routesFile.loc[cond]['route_id']
    selectedRoutes = selectedRoutes.loc[selectedRoutes!='HTM:76'].reset_index(drop=True)
#    selectedRoutes = routesFile['route_id']

    # Select trips, stop-times on a date
    mondayDate = 20180212
    selectedDate = mondayDate + day #Analysing specific date
    servicesOnDate = datesFile.loc[datesFile['date'] == selectedDate]['service_id']
    
    cond = ((tripsFile['service_id'].isin(servicesOnDate)) &
            (tripsFile['route_id'].isin(selectedRoutes))) # trips relevant to selected date and routes only
    tripsOnDate = tripsFile.loc[cond]
    
    selectedTrips = tripsOnDate['trip_id']
    stopTimesOnDate = stopTimesFile.loc[stopTimesFile['trip_id'].isin(selectedTrips)]
    stopTimesOnDate = stopTimesOnDate.merge(stopsFile,on='stop_id',how='left') #merge stop names for easy manual verification
    
    # Sort by route, direction, trip, arrival time
    # All trips of a particular trip and direction together. All trips within a
    # route-direction group ordered according to the arrival time at first stop.
    # Finally, sort stop-times within each trip by arrival time.
    stopTimesOnDate = stopTimesOnDate.merge(
                    tripsFile[['trip_id','route_id','direction_id']],on='trip_id',how='left') #merge route_id, direction_id
    sortedStopTimes = stopTimesOnDate.sort_values(['route_id','direction_id','trip_id','arrival_time'])
    
    # Convert as many to integers as possible
    sortedStopTimes.loc[:,('route_id')] = sortedStopTimes['route_id'].str[4:].astype(int)
    sortedStopTimes.loc[:,('stop_id')] = sortedStopTimes['stop_id'].str[4:].astype(int)
    sortedStopTimes.loc[:,('parent_station')] = sortedStopTimes['parent_station'].str[13:].astype(int)
    sortedStopTimes.loc[:,('trip_id')] = sortedStopTimes['trip_id'].str.split(':').str[-1].astype(int)
    return(sortedStopTimes,tripsOnDate)




#%% Extract routes from GTFS
# Station is a collection of stops as defined in GTFS 'stops'->'parent_station'
def getRouteDir(routeTopology,parentStation):
    routeDir = list(routeTopology.loc[routeTopology['parent_station']==
                                          parentStation].index.unique())
    return(routeDir)

# Get sequence positions
def getSeqPos(routeTopos,parentStation): 
    #get station position in a given route+dir;
    #if stop repeated (which is always consecutive) chooses the first position
    temp = routeTopos.loc[routeTopos['parent_station']==parentStation]
    seqPos = list(temp[~temp.index.duplicated(keep='first')]['adjustedSeq'])
#    seqPos = list(temp[~temp.index.duplicated(keep='first')]['stop_sequence'])
    return(seqPos)

def extractRoutesFromGtfs(day):
    sortedStopTimes = importGtfsRouteStops(day)[0]
    # Extract complete route details
    t_groups = sortedStopTimes.groupby(['route_id','direction_id']) #group by route,direction
    t_sampleCompleteTrips = t_groups.agg(lambda x: x.value_counts().index[0])['trip_id'] #get first trip with most frequent number of stops within that route,direction
    sampleCompleteTrips = sortedStopTimes.loc[
            sortedStopTimes['trip_id'].isin(t_sampleCompleteTrips)][
                    ['route_id','direction_id','stop_sequence',
                     'stop_id','stop_name','parent_station',
                     'stop_lat','stop_lon']].reset_index(drop=True) #get all details on complete trips
    
    routeIds = sampleCompleteTrips['route_id'].unique()
    directionIds = sampleCompleteTrips['direction_id'].unique()
    
    routeTopos = sampleCompleteTrips.drop_duplicates(
            ['route_id','direction_id','parent_station']).copy()
    routeSizes = routeTopos.groupby(['route_id','direction_id']).size() # NOTE: groupby sorts grouping variables
    adjSeq = [list(range(i)) for i in routeSizes]
    
    routeTopos = routeTopos.set_index(['route_id','direction_id'])
    routeTopos.sort_index(inplace=True) # sort to append adjustedSeq column
    routeTopos['adjustedSeq'] = np.array([j+1 for i in adjSeq for j in i])
    
    # Station Details
    allStationsInfo = sampleCompleteTrips.drop_duplicates(
            subset=['parent_station'],keep='first')[
            ['parent_station','stop_name','stop_lat','stop_lon']].set_index('parent_station')
    
    numStations = len(allStationsInfo)
    
    allStationsInfo['routeDir'] = [getRouteDir(routeTopos,x) for x in allStationsInfo.index]
    
    # Transfer stations
    transferStations = [i for i in allStationsInfo.index if len(
            set(list(zip(*allStationsInfo.loc[i]['routeDir']))[0]))>1]
    
    # Station Connections: Create Indices
    t_numRouteDir = [len(allStationsInfo.loc[x]['routeDir']) for x in allStationsInfo.index]
    t_zipped = list(zip(t_numRouteDir,list(allStationsInfo.index)))
    t_stationIndex = [x for i,x in t_zipped for _ in range(i)]
    t_routeIndex = [i_route for i_routeDir in allStationsInfo['routeDir'] for i_route,_ in i_routeDir]
    t_dirIndex = [i_dir for i_routeDir in allStationsInfo['routeDir'] for _,i_dir in i_routeDir]
    allCxns_index = pd.MultiIndex.from_arrays([t_stationIndex,t_routeIndex,t_dirIndex],
                                              names=['parent_station','route','dir'])
    
    t_seqPosLists = [getSeqPos(routeTopos,x) for x in allStationsInfo.index]
    t_seqPos = gf.flattenD2(t_seqPosLists)
    allCxns = pd.DataFrame({'seqPos':t_seqPos},allCxns_index)
    
    # L' space = next adjacent stop on each route, dir
    # P space = all subsequent stops on each route, dir
    t = []
    for srd in allCxns.index:
        t_routeDir = routeTopos.loc[srd[1:]]
        if allCxns.loc[srd]['seqPos']<routeSizes.loc[srd[1:]]:
            t += [tuple(t_routeDir.loc[t_routeDir['adjustedSeq']==
                                       allCxns.loc[srd]['seqPos']+1]['parent_station'])]
        else:
            t += [np.nan]
    allCxns['ldash'] = t
    
    t = []
    for srd in allCxns.index:
        t_routeDir = routeTopos.loc[srd[1:]]
        if allCxns.loc[srd]['seqPos']<routeSizes.loc[srd[1:]]:
            t += [tuple(t_routeDir.loc[t_routeDir['adjustedSeq']>=
                                       allCxns.loc[srd]['seqPos']+1]['parent_station'])]
        else:
            t += [np.nan]
    allCxns['pspace'] = t
    
    
#    allCxns['ldash'] = [tuple((routeTopos.loc[srd[1:]]).loc[
#                        (routeTopos.loc[srd[1:]])['stop_sequence']==
#                        allCxns.loc[srd]['seqPos']+1]['parent_station'])
#                        if (allCxns.loc[srd]['seqPos']<routeSizes.loc[srd[1:]]) 
#                        else np.nan 
#                        for srd in allCxns.index] #srd:station,route,dir
#    
#    allCxns['pspace'] = [tuple((routeTopos.loc[srd[1:]]).loc[
#                        (routeTopos.loc[srd[1:]])['stop_sequence']>=
#                        allCxns.loc[srd]['seqPos']+1]['parent_station'])
#                        if (allCxns.loc[srd]['seqPos']<routeSizes.loc[srd[1:]]) 
#                        else np.nan 
#                        for srd in allCxns.index]
    
    # Adjacency matrices
    #L' space
    t_list=[[[] for x in range(numStations)] for x in range(numStations)]
    
    for i_line in np.unique(routeTopos.index.values):
        t_current = routeTopos.loc[i_line]['parent_station'].iloc[0] #18/8/9 Added iloc because of inconsistent behaviour; see: https://stackoverflow.com/questions/51767835/inconsistent-behaviour-in-multiindex-indexing
        for t_next in routeTopos.loc[i_line]['parent_station'].iloc[1:]: #same as above
            t_i = allStationsInfo.index.get_loc(t_current)
            t_j = allStationsInfo.index.get_loc(t_next)
            t_list[t_i][t_j] += [i_line]
            t_current = t_next
    
    ldash = pd.DataFrame(t_list,index=allStationsInfo.index,
                         columns=allStationsInfo.index)
    # P space
    t_list=[[[] for x in range(numStations)] for x in range(numStations)]
    
    for i_line in np.unique(routeTopos.index.values):
        for i_seq in range(routeSizes.loc[i_line]):
            t_current = routeTopos.loc[i_line]['parent_station'].iloc[i_seq] #18/8/9 Added iloc because of inconsistent behaviour; see: https://stackoverflow.com/questions/51767835/inconsistent-behaviour-in-multiindex-indexing
            for t_next in routeTopos.loc[i_line]['parent_station'].iloc[i_seq+1:]: #same as above
                t_i = allStationsInfo.index.get_loc(t_current)
                t_j = allStationsInfo.index.get_loc(t_next)
                t_list[t_i][t_j] += [i_line]
    
    pspace = pd.DataFrame(t_list,index=allStationsInfo.index,
                         columns=allStationsInfo.index)
    
    return(allStationsInfo,numStations,routeTopos,routeSizes,pspace,transferStations)




#%% Route extraction Numpy conversion
def routeExtractionNumpyConversion(day):
    (allStationsInfo,numStations,
     routeTopos,routeSizes,
     pspace,transferStations) = extractRoutesFromGtfs(day)
    
    # Stations and lines
    #Stations
    allStations = np.array(allStationsInfo.index)
    parentStationIndexMapper = dict(zip(allStations,range(numStations)))
    # Lines
    t_routeTopos = routeTopos.loc[:,['stop_sequence','adjustedSeq','parent_station']].reset_index().values
    allLines = np.unique(t_routeTopos[:,0:2],axis=0)
    numLines = len(allLines)
    # Station-Lines
    stationLines = np.zeros((numStations,numLines),dtype=int)
    for i in range(numStations):
        for j in allStationsInfo.iloc[i]['routeDir']:
            t_1 = j==allLines[:]
            t_2 = t_1[:,0] & t_1[:,1]
            stationLines[i,t_2] = 1
            
    # Transfer stations
    transferStationsMask = np.in1d(allStations,transferStations,assume_unique=True)
    
    # Walkable transfers
    stationDists = gf.latLonDistArray(allStationsInfo['stop_lat'].values,
                                     allStationsInfo['stop_lon'].values)
    
    pWalkableTransferDist = 200 # metres
    
    t_walkable = np.argwhere(stationDists<pWalkableTransferDist)
    lspaceWalk = np.zeros((numStations,numStations))
    for i in np.unique(t_walkable):
        lspaceWalk[i,t_walkable[t_walkable[:,0]==i,1]] = 1
    lspaceWalk = lspaceWalk.astype(bool)
    walkLineId = numLines
    
    # Pspace
    numLinesPspace = np.zeros((numStations,numStations),dtype=int)
    numpyPspace = np.zeros((numStations,numStations,numLines),dtype=int)
    for i in range(numStations):
        for j in range(numStations):
            numLinesPspace[i,j] = len(pspace.iloc[i,j])
            for k in range(numLinesPspace[i,j]):
                t_1 = pspace.iloc[i,j][k]==allLines[:]
                t_2 = t_1[:,0] & t_1[:,1]
                numpyPspace[i,j,t_2] = 1
    numpyPspace = numpyPspace.astype(bool)
    
    # Traversals
    routeTopos['parentStationIndex'] = [parentStationIndexMapper[i] for i in routeTopos['parent_station']]
    pspaceTraversals = [[[] for x in range(numStations)] for x in range(numStations)]
    pspaceLines = [[[] for x in range(numStations)] for x in range(numStations)]
    for i in range(numLines):
        t_routeDirStations = np.array(routeTopos.loc[tuple(allLines[i,:])]['parentStationIndex'])
        for j in range(routeSizes.iloc[i]):
            t_or = t_routeDirStations[j]
            for k in range(routeSizes.iloc[i]):
                if k-j>0:
                    pspaceLines[t_or][t_routeDirStations[k]].append(i)
                    pspaceTraversals[t_or][t_routeDirStations[k]].append(t_routeDirStations[j+1:k])
    
    # Common lines
    #import copy
    #pspaceTraversalsA = copy.deepcopy(pspaceTraversals)
    #pspaceLinesA = copy.deepcopy(pspaceLines)
    numLinesUniquePspace = np.zeros((numStations,numStations),dtype=int)
    pspaceTraversalsNumpyList = [[[] for x in range(numStations)] for x in range(numStations)]
    
    for i in range(numStations):
        for j in range(numStations):
            if i!=j:
                t_ids,t_uniqueIds = gf.uniqueElementIds(pspaceTraversals[i][j])
                pspaceTraversals[i][j] = [pspaceTraversals[i][j][x].tolist() for x in t_uniqueIds] # because the values are the same
                pspaceLines[i][j] = [[pspaceLines[i][j][y] for y in range(len(t_ids)) if t_ids[y]==x] for x in t_uniqueIds] # because the value matters here
                numLinesUniquePspace[i,j] = len(pspaceLines[i][j])
                t_numpyTraversals = np.zeros((numLinesUniquePspace[i,j],numStations)).astype(bool)
                for k in range(numLinesUniquePspace[i,j]):
                    t_numpyTraversals[k,pspaceTraversals[i][j][k]] = True
                pspaceTraversalsNumpyList[i][j] = t_numpyTraversals
    
    # Save variables
    filename = '../vars/routeExtraction'+str(day)+'.pkl'
    with open(filename,'wb') as f:
        pickle.dump([allStations,numStations,allLines,numLines,stationLines,
                     pspaceTraversals,pspaceTraversalsNumpyList,pspaceLines,
                     numLinesPspace,numLinesUniquePspace,numpyPspace,
                     walkLineId,lspaceWalk,transferStationsMask,routeSizes,routeTopos],f)



#%% GTFS Attributes
def timeDiff(tcol): 
    #need own def to handle hours > 23
    t_hour = tcol.map(lambda x: x[:2]).astype(int)
    t_min = tcol.map(lambda x: x[3:5]).astype(int)
    t_sec = tcol.map(lambda x: x[6:8]).astype(int)
    return (t_hour.diff()*3600 + t_min.diff()*60 + t_sec.diff())

def getTime(tcol): 
    #need own def to handle hours > 23
    t_hour = tcol.map(lambda x: x[:2]).astype(int)
    t_min = tcol.map(lambda x: x[3:5]).astype(int)
    t_sec = tcol.map(lambda x: x[6:8]).astype(int)
    return (t_hour*3600 + t_min*60 + t_sec)

def gtfsAttributes(day):
    (sortedStopTimes,tripsOnDate) = importGtfsRouteStops(day)
    stationLineArrivals = sortedStopTimes.copy().reset_index(drop=True).drop_duplicates(['trip_id','route_id','direction_id','parent_station'])

    with open('../vars/routeExtraction'+str(day)+'.pkl','rb') as f:# so that extract routes doesn't have to run every time
        (allStations,numStations,allLines,numLines,stationLines,
         pspaceTraversals,pspaceTraversalsNumpyList,pspaceLines,
         numLinesPspace,numLinesUniquePspace,numpyPspace,
         walkLineId,lspaceWalk,transferStationsMask,routeSizes,routeTopos) = pickle.load(f)
    
    indexStation = pd.DataFrame({'parentStation':allStations,
                                 'id':range(numStations)}).set_index('parentStation',drop=True)
    indexLines = pd.DataFrame({'lineId':allLines[:,0],'dirId':allLines[:,1],
                               'id':range(numLines)}).set_index(['lineId','dirId'],drop=True)
    routeTopos['lineId'] = indexLines.loc[routeTopos.index]
    
    # Append
    # Line Id/Station Id
    stationLineArrivals['lineId'] = indexLines.loc[
            pd.MultiIndex.from_frame(stationLineArrivals[['route_id','direction_id']]),'id'].values
    stationLineArrivals['stationId'] = indexStation.loc[stationLineArrivals['parent_station'].values,'id'].values
    # Arrival time
    stationLineArrivals['arrivalTimeSec'] = getTime(stationLineArrivals['arrival_time'])
    stationLineArrivals['arrivalHour'] = np.array(stationLineArrivals['arrivalTimeSec']/3600).astype(int)
    
    
    ivt = stationLineArrivals.sort_values(['trip_id','arrivalTimeSec']).copy()
    ivt['ivt'] = ivt['arrivalTimeSec'].diff()
    ivt.loc[~ivt.duplicated(['trip_id']),'ivt']=np.nan
    ivt = ivt.groupby(['lineId','parent_station'])['ivt'].mean().reset_index()
    ivt = pd.merge(routeTopos,ivt,on=['lineId','parent_station']).set_index(
            ['lineId','parent_station'])['ivt']
    
    headway = stationLineArrivals[['lineId','parent_station','arrivalTimeSec','arrivalHour']].copy()
    headway = headway.sort_values(['lineId','parent_station','arrivalTimeSec'])
    headway['headway'] = headway['arrivalTimeSec'].diff()
    t_once = ~headway.duplicated(['lineId','parent_station','arrivalHour'],keep=False)
    t_first = ~headway.duplicated(['lineId','parent_station','arrivalHour'])
    headway.loc[t_once,'headway'] = 3600/2 # trip takes place only once in hour so set average waiting time as 1800s; this is not necessarily true but is a frequency-based assumption!
    headway.loc[(~t_once) & t_first,'headway'] = np.nan
    headway = headway.groupby(['lineId','arrivalHour'])['headway'].mean().reset_index()
    
    headway.to_msgpack('../vars/headway'+str(day)+'.msg')
    ivt.to_msgpack('../vars/ivt'+str(day)+'.msg')


#%% Main
if __name__ == '__main__':
    for day in range(7):
        print(day)
        routeExtractionNumpyConversion(day)
        gtfsAttributes(day)
        