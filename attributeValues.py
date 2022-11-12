### Assigning attributes


#%% Import
import pandas as pd
import numpy as np
import pickle


#%% Definitions
def getMasterChoiceSetPart(filenum):
    matPathStations = pd.read_msgpack('../vars/topoGtfs/'+str(day)+'/matPathStations'+str(filenum)+'.msg')
    matPathLines = pd.read_msgpack('../vars/topoGtfs/'+str(day)+'/matPathLines'+str(filenum)+'.msg')
    return(matPathStations,matPathLines)

def intializing(odPathLines,odPathStations,maxTransfers,numLines,walkLineId,ivtPspace):
    linesUsed = np.zeros([len(odPathLines),maxTransfers+1,numLines],dtype=int)
    transfers = np.zeros([len(odPathLines)],dtype=int)
    ivt = np.zeros([len(odPathLines)])
    stationChanges = np.zeros([len(odPathLines)])
    for x in range(len(odPathLines)):
        for y in range(len(odPathLines[x])):
            if odPathLines[x][y]==(walkLineId,):
                stationChanges[x] += 1
                continue
            ivt[x] += np.mean(ivtPspace[odPathStations[x][y],odPathStations[x][y+1],
                              odPathLines[x][y]])
            linesUsed[x,transfers[x],odPathLines[x][y]] = 1
            transfers[x] += 1
    transfers += -1
    return(linesUsed,transfers,ivt,stationChanges)

def noWalkStuff(odPathLines,odPathStations,allStations,allLines,walkLineId):
    noWalkLines = np.array([[list(allLines[j,0]) for j in i if j!=(walkLineId,)] for i in odPathLines]) # Remove all walking links [(48,)] for easy WT, num Transfer
    noWalkStations = [[(odPathStations[i][j],odPathStations[i][j+1]) 
                            for j in range(len(odPathLines[i])) if odPathLines[i][j]!=(walkLineId,)] 
                            for i in range(len(odPathLines))]
    checkInStations = np.array([[allStations[j[0]] for j in i] for i in noWalkStations]) # change to original IDs to match with AFC data
    checkOutStations = np.array([[allStations[j[1]] for j in i] for i in noWalkStations])
    return(noWalkLines,noWalkStations,checkInStations,checkOutStations)

def getActiveLinesWt(timeH,headway,linesUsed,transfers,numLines):
    linesF = np.zeros(numLines)
    activeLines = headway.loc[timeH].index.values
    activeLinesF = (3600/headway.loc[timeH].values).T[0]
    linesF[activeLines] = activeLinesF
    linesUsedF = linesUsed * linesF
    linesUsedFSum = np.sum(linesUsedF,axis=2)
    nonZeroFs = np.sum(linesUsedFSum>0,axis=1)
    validRoutes = np.equal(nonZeroFs,(transfers+1))
    linesUsedFSum = linesUsedFSum[validRoutes]
    wt = (3600/linesUsedFSum)/2
    wt[np.isinf(wt)] = np.nan # wherever frequency is zero (because transfer is not needed) change to nan
    return(wt,validRoutes)
    
def listStore(allStations,infoOrDT,wt,transfers,ivt,stationChanges,
                  noWalkLines,checkInStations,checkOutStations,validRoutes,
                  listOrigin,listDest,listTime,
                  listWt1,listWt2,listWt3,listWtTotal,
                  listTransfers,listIvt,listStationChanges,
                  listLines,listCheckIns,listCheckOuts,
                  listWt1Ratio,listWt2Ratio,listWt3Ratio,listWtRatio,
                  listIvtRatio,listTransfersDiff,listStationChangesDiff):
    listOrigin += allStations[infoOrDT[:,0].astype(int)].tolist()#infoOrDT[:,0].tolist()
    listDest += allStations[infoOrDT[:,1].astype(int)].tolist()#infoOrDT[:,1].tolist()
    listTime += infoOrDT[:,2].tolist()
    listWt1 += wt[:,0].tolist()
    listWt2 += wt[:,1].tolist()
    listWt3 += wt[:,2].tolist()
    listWtTotal += np.nansum(wt,axis=1).tolist() # nansum considers nans as zeros
    listTransfers += transfers[validRoutes].tolist()
    listIvt += ivt[validRoutes].tolist()
    listStationChanges += stationChanges[validRoutes].tolist()
    listLines += noWalkLines[validRoutes].tolist()
    listCheckIns += checkInStations[validRoutes].tolist()
    listCheckOuts += checkOutStations[validRoutes].tolist()
    
    listWt1Ratio += (wt[:,0]/np.nanmin(wt[:,0])).tolist() # nanmin not reqd because only non nan selected
    if any(~np.isnan(wt[:,1])):
        listWt2Ratio += (wt[:,1]/np.nanmin(wt[:,1])).tolist()
    else:
        listWt2Ratio += (wt[:,1]).tolist()
    
    if any(~np.isnan(wt[:,2])):
        listWt3Ratio += (wt[:,2]/np.nanmin(wt[:,2])).tolist()
    else:
        listWt3Ratio += (wt[:,2]).tolist()
    
    listWtRatio += (np.nansum(wt,axis=1)/np.min(np.nansum(wt,axis=1))).tolist() # nansum considers nans as zeros
    listIvtRatio += (ivt[validRoutes]/np.min(ivt[validRoutes])).tolist()
    listTransfersDiff += (transfers[validRoutes]-np.min(transfers[validRoutes])).tolist()
    listStationChangesDiff += (stationChanges[validRoutes]-np.min(stationChanges[validRoutes])).tolist()
    return(listOrigin,listDest,listTime,
                listWt1,listWt2,listWt3,listWtTotal,
                listTransfers,listIvt,listStationChanges,
                listLines,listCheckIns,listCheckOuts,
                listWt1Ratio,listWt2Ratio,listWt3Ratio,listWtRatio,
                listIvtRatio,listTransfersDiff,listStationChangesDiff)

def getAttributes(day):
    #--------------------------------------------------------------------------
    # Get day dependent variables
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
    
    #obs['choiceId'] = indexStation.loc[obs['choice2'].values].values
    #routesFile = pd.read_csv('../gtfs/routes.txt')
    #busRoutes = indexLines.loc[routesFile.loc[routesFile['agency_id']=='HTMBUZZ']['route_id'].str[4:].astype(int).values].values
    #tramRoutes = indexLines.loc[routesFile.loc[routesFile['agency_id']=='HTM']['route_id'].str[4:].astype(int).values].values
    #sum(np.array([all(np.isin(np.where(stationLines[int(i),:]==1)[0],busRoutes)) for i in obs['choiceId'].values if ~np.isnan(i)]))
    
    # Select OD pairs for this day
    uniqueOd_day = uniqueOd.loc[uniqueOd['day']==day].copy()
    uniqueOd_day['orId'] = indexStation.loc[uniqueOd_day['origin'].values].values
    uniqueOd_day['destId'] = indexStation.loc[uniqueOd_day['dest'].values].values
    uniqueOd_day = uniqueOd_day.sort_values(['orId','destId'])
    uniqueOd_day = uniqueOd_day.loc[~(uniqueOd_day['orId'].isnull()|uniqueOd_day['destId'].isnull())]
    uniqueOd_day = uniqueOd_day.reset_index(drop=True)
    
    # Get headways for WT
    headway = pd.read_msgpack('../vars/headway'+str(day)+'.msg')
    headway = headway.set_index(['arrivalHour','lineId'])
    
    # Get IVT Pspace
    ivtList = pd.read_msgpack('../vars/ivt'+str(day)+'.msg')
    ivtLspace = np.zeros((numStations,numStations,numLines),dtype=int)
    for i in ivtList.index.values:
        if np.isnan(ivtList[i]):
            j = i
            continue
        else:
            ivtLspace[(allStations==j[1]),(allStations==i[1]),i[0]] = ivtList[i]
            j = i
    
    ivtPspace = np.zeros((numStations,numStations,numLines),dtype=int)
    for i in range(numLines):
        t_line = routeTopos.loc[tuple(allLines[i,:]),'parent_station']
        for j in range(routeSizes.iloc[i]-1):
            t_dist = 0
            for k in range(j+1,routeSizes.iloc[i]):
                if t_dist == 0:
                    t_dist += ivtLspace[(allStations==t_line.iloc[j]),(allStations==t_line.iloc[k]),i]
                else:
                    t_dist += ivtLspace[(allStations==t_line.iloc[k-1]),(allStations==t_line.iloc[k]),i]
                ivtPspace[(allStations==t_line.iloc[j]),(allStations==t_line.iloc[k]),i] = t_dist
                
    ivtMat = np.ma.masked_equal(ivtPspace,0).mean(axis=2)
    ivtMat.mask = np.ma.nomask
    
    
    # Get tram stop indices for walkable dest stops
    routesFile = pd.read_csv('../gtfs/routes.txt',usecols = ['route_id','agency_id'])
    tramLines = routesFile.loc[routesFile['agency_id']=='HTM','route_id'].str[4:].astype(int).values
    tramStIndex = routeTopos.loc[routeTopos.reset_index()['route_id'].isin(tramLines).values,
                                 'parentStationIndex'].unique()
    
    
    
    #--------------------------------------------------------------------------
    # Get Master Choice Set for OD
    names = ['origin','dest','arrivalHour',
             'wt1','wt2','wt3','wtTotal',
             'transfers','ivt','stationChanges',
             'lines','checkIns','checkOuts',
             'wt1Ratio','wt2Ratio','wt3Ratio','wtRatio',
             'ivtRatio','transfersDiff','stationChangesDiff'] # choice set column names
    
    # Initialize storage 
    lenInfo = 0
    listOrigin = []
    listDest = []
    listTime = []
    listWt1 = []
    listWt2 = []
    listWt3 = []
    listWtTotal = []
    listTransfers = []
    listIvt = []
    listStationChanges = []
    listLines = []
    listCheckIns = []
    listCheckOuts = []
    listWt1Ratio = []
    listWt2Ratio = []
    listWt3Ratio = []
    listWtRatio = []
    listIvtRatio = []
    listTransfersDiff = []
    listStationChangesDiff = []
    
    i = 0
    orId = uniqueOd_day.iloc[i]['orId'] # orIds are sorted; so import the first relevant files
    filenum = int(orId/maxRowsInFile)
    matPathStations,matPathLines = getMasterChoiceSetPart(filenum+1)
    odsize = len(uniqueOd_day)
    
    for i in range(len(uniqueOd_day)):
        print(str(i)+'/'+str(odsize))
    
        # Get or, dest, hr
        orId = uniqueOd_day.iloc[i]['orId']
        destId = uniqueOd_day.iloc[i]['destId']
        timeH = uniqueOd_day.iloc[i]['arrivalHour']
        
        
        
        # If the file doesn't contain this orId, import next file
        if orId >= (filenum+1)*maxRowsInFile:
            filenum += 1
            matPathStations,matPathLines = getMasterChoiceSetPart(filenum+1)
        
        # Get path alts: stations, lines
        odPathStations=()
        odPathLines=()
        walkableDestStops = np.where(lspaceWalk[int(destId),:])[0]
        for i in walkableDestStops:
            if i in tramStIndex:
                odPathStations = odPathStations + matPathStations.iloc[int(orId),int(i)]
                odPathLines = odPathLines + matPathLines.iloc[int(orId),int(i)]
        
        
#        odPathStations = matPathStations.iloc[int(orId),int(destId)]
#        odPathLines = matPathLines.iloc[int(orId),int(destId)]
        
        # Get ivt, wt
        (linesUsed,transfers,ivt,stationChanges) = intializing(odPathLines,odPathStations,maxTransfers,
                                                                numLines,walkLineId,ivtPspace)
        (noWalkLines,noWalkStations,checkInStations,checkOutStations) = noWalkStuff(
                odPathLines,odPathStations,allStations,allLines,walkLineId)
        wt,validRoutes = getActiveLinesWt(timeH,headway,linesUsed,transfers,numLines)
        
        # If this OD doesn't have any valid route: ignore
        if not validRoutes.any():
            continue
        
        t_len = np.sum(validRoutes)
        lenInfo += t_len
        infoOrDT = np.tile(np.array([orId,destId,timeH]),(t_len,1))
        
        (listOrigin,listDest,listTime,
         listWt1,listWt2,listWt3,listWtTotal,
         listTransfers,listIvt,listStationChanges,
         listLines,listCheckIns,listCheckOuts,
         listWt1Ratio,listWt2Ratio,listWt3Ratio,listWtRatio,
         listIvtRatio,listTransfersDiff,
         listStationChangesDiff) = listStore(allStations,infoOrDT,wt,transfers,ivt,stationChanges,
                               noWalkLines,checkInStations,checkOutStations,validRoutes,
                               listOrigin,listDest,listTime,
                               listWt1,listWt2,listWt3,listWtTotal,
                               listTransfers,listIvt,listStationChanges,
                               listLines,listCheckIns,listCheckOuts,
                               listWt1Ratio,listWt2Ratio,listWt3Ratio,listWtRatio,
                               listIvtRatio,listTransfersDiff,listStationChangesDiff)
    
    #--------------------------------------------------------------------------
    # Store
    attributes = pd.DataFrame(data=[listOrigin,listDest,listTime,
                                    listWt1,listWt2,listWt3,listWtTotal,
                                    listTransfers,listIvt,listStationChanges,
                                    listLines,listCheckIns,listCheckOuts,
                                    listWt1Ratio,listWt2Ratio,listWt3Ratio,listWtRatio,
                                    listIvtRatio,listTransfersDiff,listStationChangesDiff])
    attributes = attributes.transpose()
    attributes.columns = names
    attributes.to_hdf('../vars/attributes/attributes1'+str(day)+'.h5',key='df',mode='w')

def getBestRoutes(days):
    bestRoutes = []
    for day in days:
        att = pd.read_hdf('../vars/attributes/attributes1'+str(day)+'.h5')
        att['day'] = day
        indicators = ['wtRatio','ivtRatio','transfersDiff']
        att['perf'] = att[indicators].sum(axis=1)
        idx = att.groupby(['origin','dest','arrivalHour']).apply(lambda att:att.perf.idxmin()) # index of minimum for each group
        bestRoutes += [att.loc[idx,['origin','dest','arrivalHour','day',
                                    'wtTotal','ivt','transfers','lines']].reset_index(drop=True)]
    pd.concat(bestRoutes).to_hdf('../vars/attributes/bestRoutes1.h5',key='df',mode='w')


#%% Main
if __name__ == '__main__':
    # Access station: Alternatives and choices
    alts = pd.read_excel('../Data/20190715 DT Tramtrips - onlytramparents.xlsx')
    alts = alts.rename(columns={'Userid':'persId',
                                'potential start station':'origin','endstation parent':'dest',
                                'hour':'arrivalHour','day':'day'})
    alts = alts.loc[~alts['persId'].isin([306,1306])] # (sometimes) NaNs in the hours column of these two users
#    alts = alts.loc[alts['origin'].isin(allStations)]
    uniqueOd = alts.drop_duplicates(['origin','dest','arrivalHour','day'])[['origin','dest','arrivalHour','day']]
    
    # Settings
    maxRowsInFile = 100
    maxTransfers = 2
    
    for day in uniqueOd['day'].unique():
        print(day)
        getAttributes(day)
    
    getBestRoutes(uniqueOd['day'].unique())
#    getBestRoutes([0])

#%%
#if odPathStations:
#    minTransfers += [min(np.array([len(k) for k in odPathStations])-1)]
#else:
#    minTransfers += [np.nan]
#
#minTransfers = np.array(minTransfers)


