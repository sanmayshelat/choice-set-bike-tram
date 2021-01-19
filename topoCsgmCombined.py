#%% Author and Title Information
# Title: Topological route choice generation: Choice set generation
# Author: Sanmay Shelat
# Institute: Dept. of Transport and Planning, Delft University of Technology

#%% Imports
import time
start_time = time.time()
import pandas as pd
import numpy as np
from myGlobalFxns import flattenD2
from multiprocessing import Pool
import pickle


#%% Definitions
# Initialize stop visited/traversed, lines used
def booleanInit(routeInfo,numStations,numLines):
    boolVisitedStations = np.zeros((numStations,)).astype(bool) # converting to boolean arrays for easy indexing
    boolVisitedStations[routeInfo[0]] = True
    boolUsedLines = np.zeros((numLines+1,)).astype(bool) # numlines+1 required to allow walking to be marked true in the next line
    boolUsedLines[flattenD2(routeInfo[1])] = True
    boolUsedLines=boolUsedLines[:-1] # remove the last one which was there only to allow walking to work in previous line
    boolTraversedStations = np.zeros((numStations,)).astype(bool)
    boolTraversedStations[flattenD2(routeInfo[2])] = True
    return(boolVisitedStations,boolUsedLines,boolTraversedStations)

# Logic 1: No stop repetitions
    # Passengers should not transfer at stops they have already visited
def noVisitingVisited(boolVisitedStations,cxnsMask): # boolVisitedStations is boolean array
    cxnsMask[boolVisitedStations] = False
    return(cxnsMask)

# Logic 2: No circular routes
    # Passengers should not visit stations they have already traversed through
def noVisitingTraversed(boolTraversedStations,cxnsMask): # boolTraversedStations is boolean array
    cxnsMask[boolTraversedStations] = False
    return(cxnsMask)

# Get traversed stations, split paths, get records
    # Get traversed stations from pre-stored array; based on this split the cnxs
    # So cxns can be repeated if there are more than one paths between the origin and cxn
    # Also create path records now because the paths need to be stored
def getTraversedNodes(pspaceTraversalsNumpyList,numLinesUniquePspace,pspaceTraversals,pspaceLines,orderedPreviousRouteInfo,orIndex,cxns): # also need to return used lines,used stations for path record
    cxnsNew = np.repeat(cxns,numLinesUniquePspace[orIndex][cxns])
    traversedStationsNumpy = np.vstack([pspaceTraversalsNumpyList[orIndex][i] for i in cxns])
    
    recordVisitedStations = np.column_stack((np.tile(orderedPreviousRouteInfo[0],(len(cxnsNew),1)),cxnsNew)).tolist()
    recordUsedLines = [orderedPreviousRouteInfo[1] + [j] for i in cxns for j in pspaceLines[orIndex][i]]
    recordTraversedStations = [orderedPreviousRouteInfo[2] + [j] for i in cxns for j in pspaceTraversals[orIndex][i]]
    return(cxnsNew,traversedStationsNumpy,recordVisitedStations,recordUsedLines,recordTraversedStations)

# Logic 3: No transfers to common lines
    # Should not transfer to one of the lines that connected the previous origin to this origin.
    # This allows passengers to shift their waiting time; that is, if two lines take them from 
    # previous origin to this origin, they can take the first one and then wait for the other line
    # at a farther station
    # Should also not transfer to a common line. Passengers should not be allowed to get down of one
    # of the lines connecting their previous origin to this origin, only to board a line that
    # is common (connects the same sequence of stops) as the others.
    # Delete records too, for the deleted cxns 
def noCommonLineTransfers(orderedPreviousLines,traversedStations,recordVisitedStations,recordUsedLines,recordTraversedStations,cxns,walkLineId): # visitedStations is ordered 'list' of stations visited
    t_mask = []
    if orderedPreviousLines[-1]==[walkLineId]:
        t_commonLinesPreviousTwo = orderedPreviousLines[-2]
        for i in range(len(cxns)):
            t = np.array([j in recordUsedLines[i][-1] for j in t_commonLinesPreviousTwo])
            t_mask.append(~np.all(t))
            recordUsedLines[i][-3] = np.array(recordUsedLines[i][-3])[~t].tolist()
    else:
        t_commonLinesPreviousTwo = orderedPreviousLines[-1]
        for i in range(len(cxns)):
            t = np.array([j in recordUsedLines[i][-1] for j in t_commonLinesPreviousTwo])
            t_mask.append(~np.all(t))
            recordUsedLines[i][-2] = np.array(recordUsedLines[i][-2])[~t].tolist()

    traversedStations = traversedStations[t_mask]
    cxns = cxns[t_mask]
    recordVisitedStations = np.array(recordVisitedStations)[t_mask].tolist() # convert to array for faster mask-based filtering
    recordUsedLines = np.array(recordUsedLines)[t_mask].tolist()
    recordTraversedStations = np.array(recordTraversedStations)[t_mask].tolist()
    return(cxns,traversedStations,recordVisitedStations,recordUsedLines,recordTraversedStations)

# Logic 4: No traversing through traversed + visited nodes
    # Passengers should not traverse through nodes that they have already visited or 
    # traversed through
    # Delete records too, for the deleted cxns    
def noTraversingVisitedTraversed(boolVisitedStations,boolTraversedStations,traversedStations,recordVisitedStations,recordUsedLines,recordTraversedStations,cxns): # boolVisitedStations, boolTraversedStations, traversedStations are boolean arrays
    # traversedStations should not overlap with boolVisitedStations,boolTraversedStations
    t_mask = ~np.any((traversedStations&(boolVisitedStations|boolTraversedStations)),axis=1) # traversedStations,boolVisitedStations,boolTraversedStations automatically upgraded to int
    
    traversedStations = traversedStations[t_mask]
    cxns = cxns[t_mask]
    recordVisitedStations = np.array(recordVisitedStations)[t_mask].tolist() # convert to array for faster mask-based filtering
    recordUsedLines = np.array(recordUsedLines)[t_mask].tolist()
    recordTraversedStations = np.array(recordTraversedStations)[t_mask].tolist()
    return(cxns,traversedStations,recordVisitedStations,recordUsedLines,recordTraversedStations)

# Logic 5: Don't walk to station with subset of lines
def noWalkingToSameLines(stationLines,previousStop,walkableStops):
    walkableStops[walkableStops] = np.any((stationLines[previousStop,:]-stationLines[walkableStops,:])<0,axis=1)
    return(walkableStops)

# Get logical walkable transfers
def getLogicalWalkableTransfers(stationLines,thisCxnVisitedStations,thisCxnTraversedStations,lspaceWalk,cxn):
    walkCxnsMask = lspaceWalk[cxn,:].copy() # ??? Is the copy needed? lspaceWalk seems to change without it
    walkCxnsMask = noVisitingVisited(thisCxnVisitedStations,walkCxnsMask) # cannot visit visited stops
    walkCxnsMask = noVisitingTraversed(flattenD2(thisCxnTraversedStations),walkCxnsMask) # cannot visit traversed stops
    walkCxnsMask = noWalkingToSameLines(stationLines,cxn,walkCxnsMask) # don't walk to stops with a subset of lines
    return(walkCxnsMask)

# Add walked transfer path records
def addWalkedPathRecords(walkLineId,walkCxnsMask,
                         thisCxnVisitedStations,thisCxnUsedRouteDirs,thisCxnTraversedStations,
                         walkCxns,walkStations,walkRouteDirs,walkTraversed):
    walkCxns += np.where(walkCxnsMask)[0].tolist()
    walkStations += [thisCxnVisitedStations+[j] for j in np.where(walkCxnsMask)[0].tolist()]
    walkRouteDirs += [thisCxnUsedRouteDirs+[[walkLineId]] for _ in range(np.sum(walkCxnsMask))]
    walkTraversed += [thisCxnTraversedStations] * np.sum(walkCxnsMask)
    return(walkCxns,walkStations,walkRouteDirs,walkTraversed)

# Initialization of final matrices
def initializeMats(allStations,numStations):
    matPathStations = pd.DataFrame(columns=allStations,index=allStations,data=[[[] for _ in range(numStations)] for _ in range(numStations)])
    matPathLines = pd.DataFrame(columns=allStations,index=allStations,data=[[[] for _ in range(numStations)] for _ in range(numStations)])
    matPathTraversed = pd.DataFrame(columns=allStations,index=allStations,data=[[[] for _ in range(numStations)] for _ in range(numStations)])
    return(matPathStations,matPathLines,matPathTraversed)







def main(day):
    #%% Import basic variables
    with open('../vars/routeExtraction'+str(day)+'.pkl','rb') as f:# so that extract routes doesn't have to run every time
        (allStations,numStations,allLines,numLines,stationLines,
         pspaceTraversals,pspaceTraversalsNumpyList,pspaceLines,
         numLinesPspace,numLinesUniquePspace,numpyPspace,
         walkLineId,lspaceWalk,transferStationsMask,routeSizes,routeTopos) = pickle.load(f)
    
    #%% Parameters
    pMaxTransfers = 2
    maxRowsInFile = 100
    
    #%% Initialize
    matPathStations,matPathLines,matPathTraversed = initializeMats(allStations,numStations)
    
    #%% For every origin
    for originIndex in range(numStations):
        #originIndex = 0
        print(originIndex,numStations)
        
        # Prevent low memory issues by managing variable size
        if (originIndex%maxRowsInFile==0) & (originIndex>0):
            # Save to file
            matPathLines.to_msgpack('../vars/topoGtfs/'+str(day)+'/matPathLines'+str(int(originIndex/maxRowsInFile))+'.msg')
            matPathStations.to_msgpack('../vars/topoGtfs/'+str(day)+'/matPathStations'+str(int(originIndex/maxRowsInFile))+'.msg')
            matPathTraversed.to_msgpack('../vars/topoGtfs/'+str(day)+'/matPathTraversed'+str(int(originIndex/maxRowsInFile))+'.msg')
            # Initialize
            matPathStations,matPathLines,matPathTraversed = initializeMats()
            print('File saved')
        
        
        level = 0
        # Initialize route information
        recordVisitedStations = [originIndex] # keep as list becuase will increase in length iteratively
        recordUsedLines = []
        recordTraversedStations = []
        routeInfo = [(recordVisitedStations,recordUsedLines,recordTraversedStations)]
        
        #%% For each level
        while level <= pMaxTransfers:
            level += 1
            # Make space for new route information
            prevRouteInfo = routeInfo
            routeInfo = []
            
            #%% For every origin in this level
            for orNum in range(len(prevRouteInfo)):
                #orNum = 5
                orIndex = prevRouteInfo[orNum][0][-1] # assign previous station for easy access
                boolVisitedStations,boolUsedLines,boolTraversedStations = booleanInit(prevRouteInfo[orNum],numStations,numLines) # converting to boolean arrays for easy indexing
                cxnsMask = numLinesPspace[orIndex,:]>0
                #print(len(np.argwhere(cxnsMask)))
                if np.sum(cxnsMask)==0: continue
                #%% Logic checks 1,2
                if level > 1:
                    cxnsMask = noVisitingVisited(boolVisitedStations,cxnsMask)
                    #print(len(np.argwhere(cxnsMask)))
                    cxnsMask = noVisitingTraversed(boolTraversedStations,cxnsMask)
                    #print(len(np.argwhere(cxnsMask)))
                    if np.sum(cxnsMask)==0: continue
                #%% Get traversed stations, split paths, get records
                cxns, = np.where(cxnsMask)
                (cxns,traversedStations,
                 recordVisitedStations,recordUsedLines,recordTraversedStations) = getTraversedNodes(
                         pspaceTraversalsNumpyList,numLinesUniquePspace,pspaceTraversals,pspaceLines,
                         prevRouteInfo[orNum],orIndex,cxns)
                
                #%% Logic checks 3,4
                if level > 1:
                    (cxns,traversedStations,
                     recordVisitedStations,recordUsedLines,recordTraversedStations) = noCommonLineTransfers(
                             prevRouteInfo[orNum][1],traversedStations,
                             recordVisitedStations,recordUsedLines,recordTraversedStations,cxns,walkLineId)
                    (cxns,traversedStations,
                     recordVisitedStations,recordUsedLines,recordTraversedStations) = noTraversingVisitedTraversed(
                             boolVisitedStations,boolTraversedStations,traversedStations,
                             recordVisitedStations,recordUsedLines,recordTraversedStations,cxns)
                    #print(len(np.argwhere(cxnsMask)))
                    if np.sum(cxnsMask)==0: continue                
                
                #%% Get walkable transfer path records (uses stop based path records)
                walkCxns = []
                walkStations = []
                walkRouteDirs = []
                walkTraversed = []
                for i in range(len(cxns)):
                    walkCxnsMask = getLogicalWalkableTransfers(stationLines,recordVisitedStations[i],recordTraversedStations[i],lspaceWalk,cxns[i])
                    walkCxns,walkStations,walkRouteDirs,walkTraversed = addWalkedPathRecords(walkLineId,walkCxnsMask,
                                                                        recordVisitedStations[i],recordUsedLines[i],recordTraversedStations[i],
                                                                        walkCxns,walkStations,walkRouteDirs,walkTraversed)
                    
                #%% Assign paths to origin-connections
                for i in range(len(cxns)):
                    # Assign paths to origin-connections
                    matPathStations.iloc[originIndex,cxns[i]].append(recordVisitedStations[i])
                    matPathLines.iloc[originIndex,cxns[i]].append(recordUsedLines[i])
                    matPathTraversed.iloc[originIndex,cxns[i]].append(recordTraversedStations[i])
                                
                #%% Keep only connections which are transfer stations
                t_mask = transferStationsMask[cxns]
                recordVisitedStations = np.array(recordVisitedStations)[t_mask].tolist() # convert to array for faster mask-based filtering
                recordUsedLines = np.array(recordUsedLines)[t_mask].tolist()
                recordTraversedStations = np.array(recordTraversedStations)[t_mask].tolist()
                #print(len(np.argwhere(cxnsMask)))
    
                #%% Merge walked connections
                recordVisitedStations += walkStations
                recordUsedLines += walkRouteDirs
                recordTraversedStations += walkTraversed
    
                
                #%% Add route info
                routeInfo += zip(recordVisitedStations,recordUsedLines,recordTraversedStations)
                
    # Last save to file
    matPathLines.to_msgpack('../vars/topoGtfs/'+str(day)+'/matPathLines'+str(int(originIndex/maxRowsInFile)+1)+'.msg')
    matPathStations.to_msgpack('../vars/topoGtfs/'+str(day)+'/matPathStations'+str(int(originIndex/maxRowsInFile)+1)+'.msg')
    matPathTraversed.to_msgpack('../vars/topoGtfs/'+str(day)+'/matPathTraversed'+str(int(originIndex/maxRowsInFile)+1)+'.msg')
    #%% Timer
    print("--- %s seconds ---" % (time.time() - start_time))


#%% Main
if __name__ == '__main__':
#    for day in range(7):
#        main(day)
    with Pool(7) as pool:
        pool.map(main,range(7))

