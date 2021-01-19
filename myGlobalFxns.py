# Global function definitions
# Generic functions useful in multiple places
import numpy as np
from math import sin,cos,sqrt,atan2,radians

def with_next(iterable): 
    #Yield (current, next_item) tuples for each item in iterable.
    iterator = iter(iterable)
    current = next(iterator)
    for next_item in iterator:
        yield current, next_item
        current = next_item

def flattenD2(var): # flattener for nesting depth of 2
    return([i for x in var for i in x])

#!!! from https://stackoverflow.com/a/19412565/9995225
def latLonDist(lat1,lon1,lat2,lon2): # calculate lat long dist (lat,lons in degrees)
    R = 6373.0*1000 #earth radius in m
    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat1)
    lon2 = radians(lon1)
    difflon = lon2 - lon1
    difflat = lat2 - lat1
    
    a = sin(difflat / 2)**2 + cos(lat1) * cos(lat2) * sin(difflon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    return(R*c)
    
def latLonDistArray(lat,lon): # give list/array of lat,lon in degrees get ndarray of distance in metres
    R = 6373.0*1000 #earth radius in m
    t_len = len(lat)
    latr = np.tile(np.radians(lat),(t_len,1))
    lonr = np.tile(np.radians(lon),(t_len,1))
    difflat = latr-latr.T
    difflon = lonr-lonr.T
    a = (np.sin(difflat/2))**2 + (np.cos(latr))*(np.cos(latr.T))*((np.sin(difflon/2))**2)
    c = 2*(np.arctan2(np.sqrt(a),np.sqrt(1-a)))
    return(R*c)

def uniqueElementIds(alist):
    uniqueList = []
    uniqueIds = []
    alistIds = list(range(len(alist)))
    for i,x in enumerate(alist):
        checkNotUnique = [(j,y) for j,y in enumerate(uniqueList) if np.all(x==y)]
        if checkNotUnique: # i.e. x is NOT unique
            alistIds[i] = checkNotUnique[0][0]
        else:
            uniqueList.append(x)
            uniqueIds.append(i)
    return(alistIds,uniqueIds)

def checkIntersection(list1,list2):
    return(any(x in max(list1,list2,key=len) for x in min(list1,list2,key=len)))