### Prepare Alternatives

#%% Import
import pandas as pd
import numpy as np
import itertools

sdfsdf
#%% Definitions
def getAllAlts():
#    obs = pd.read_excel('../Data/20190429 DT halte choice2.xlsx')
#    obs = obs.rename(columns={'Persid':'persId','parentopstap (new)':'choice'})
    alts = pd.read_excel('../Data/20190715 DT Tramtrips - onlytramparents.xlsx')
#    alts = alts.rename(columns={'Userid':'persId',
#                                'potential start station':'origin','endstation parent':'dest',
#                                'accessdistance':'dist',
#                                'hour':'arrivalHour','day':'day'})
    alts = alts.rename(columns={'Userid':'persId',
                                'potential start station':'origin','endstation parent':'dest',
                                'chosenstartstation parent':'chosenOr',
                                'accessdistance':'dist',
                                'hour':'arrivalHour','day':'day'})
    bestRoutes = pd.read_hdf('../vars/attributes/bestRoutes1.h5').sort_index()
    alts = pd.merge(alts,bestRoutes,on=['origin','dest','arrivalHour','day'],how='left') #merging left only because current dataset has some faulty origins,destinations that are not in routeTopos
    
    alts = alts.loc[~alts['transfers'].isnull()] # could also just merge 'inner' in above line
    alts['tt'] = alts['wtTotal'] + alts['ivt']
    alts['dist'] = alts['dist']*1000
    
    # Get ratios/diff
    altGroups = alts.groupby(['persId'])
    alts['distRatio'] = altGroups['dist'].apply(lambda x: x/x.min())
    alts['ttRatio'] = altGroups['tt'].apply(lambda x: x/x.min())
    alts['transfersDiff'] = altGroups['transfers'].apply(lambda x: x-x.min())
    
    alts['distDiff'] = altGroups['dist'].apply(lambda x: x-x.min())
    alts['ttDiff'] = altGroups['tt'].apply(lambda x: x-x.min())
    
    
    # Get chosen alternative
#    alts = pd.merge(alts,obs[['persId','choice']],how='left',
#                    left_on=['persId','origin'],right_on=['persId','choice'])
#    alts.loc[~np.isnan(alts['choice']),'choice'] = 1
#    alts.loc[np.isnan(alts['choice']),'choice'] = 0
    alts.loc[alts['origin']==alts['chosenOr'],'choice'] = 1
    alts.loc[alts['origin']!=alts['chosenOr'],'choice'] = 0
    alts = alts.loc[alts['persId'].isin(alts.loc[alts['choice']==1,'persId'])] # remove people whose observed choice has been removed
    
    # Remove dom
    alts = alts.loc[~getDominatedAlts(alts)]
    notWeirdPeople = alts.loc[alts['choice']==1,'persId'].values
    alts = alts.loc[alts['persId'].isin(notWeirdPeople)]
    
    # Final prep
    alts['size'] = alts.groupby('persId')['persId'].transform('size')
    alts = alts[['persId','size','origin','dest','choice',
                 'ivt','wtTotal',
                 'dist','tt','transfers',
                 'distRatio', 'ttRatio', 'transfersDiff',
                 'distDiff', 'ttDiff','lines']]
    return(alts)
    
def getDominated(persAlts):
    a = persAlts.values
    dom = []
    for i in range(len(persAlts)):
        check = np.all(a[i,:]>=a,axis=1)
        check[i] = False # for the same index will always return true for soft dom; so set to false
        dom += [any(check)]
    return(dom)

def getDominatedAlts(alts):
    dominated = []
    for pers in alts['persId'].unique():
        persAlts = alts.loc[alts['persId']==pers,['dist','tt','transfers']]
        dominated += getDominated(persAlts)
    return(np.array(dominated))

def thresholdFilter(altSet,prevAtt,minIndicator,attributeTypes):
    for t in range(len(prevAtt)):
        altSet = altSet.loc[altSet[attributeTypes[prevAtt[t]]]<=minIndicator[t]]
    return(altSet)

def thresholdFilter2(altSet,att,minIndicator,attributeTypes):
    for t in att:
        altSet = altSet.loc[altSet[attributeTypes[t]]<=minIndicator[t]]
    return(altSet)

def indicatorCalc(altSet,threshold,att,attributeTypes):
    # In choice set or not (less is in)
    choiceSet = altSet[attributeTypes[att]]<=threshold
    notChoiceSet = altSet[attributeTypes[att]]>threshold
    # Observed or not
    notObserved = altSet['choice']==0
    observed = altSet['choice']==1
    
    truePositives = sum(choiceSet & observed)
    truePositivesPlusFalseNegatives = sum(observed)
#    coverage = truePositives/truePositivesPlusFalseNegatives
#    print(coverage)
    
    trueNegativesPlusFalsePositives = sum(notObserved)
    trueNegatives = sum(notChoiceSet & notObserved)
#    efficiency = trueNegatives/trueNegativesPlusFalsePositives
#    print(efficiency)
    
#    indicator = coverage - efficiency
#    print(indicator)
    return(truePositives,truePositivesPlusFalseNegatives,
           trueNegativesPlusFalsePositives,trueNegatives)

def indicatorCalc2(altSet,threshold,att,attributeTypes):
    # In choice set or not (less is in)
    newChoiceSet = altSet[attributeTypes[att]]<=threshold
    personsInCS = newChoiceSet.loc[newChoiceSet['choice']==1]['persId'].values
    newChoiceSet = newChoiceSet.loc[newChoiceSet['persId'].isin(personsInCS)]
    
    choiceSet = newChoiceSet[attributeTypes[att]]<=threshold
    notChoiceSet = altSet[attributeTypes[att]]>threshold
    # Observed or not
    notObserved = altSet['choice']==0
    observed = altSet['choice']==1
    
    truePositives = sum(choiceSet & observed)
    truePositivesPlusFalseNegatives = sum(observed)
#    coverage = truePositives/truePositivesPlusFalseNegatives
#    print(coverage)
    
    trueNegativesPlusFalsePositives = sum(notObserved)
    trueNegatives = sum(notChoiceSet & notObserved)
#    efficiency = trueNegatives/trueNegativesPlusFalsePositives
#    print(efficiency)
    
#    indicator = coverage - efficiency
#    print(indicator)
    return(truePositives,truePositivesPlusFalseNegatives,
           trueNegativesPlusFalsePositives,trueNegatives)


def finalIndicatorCalc(alts,choiceSet):
    truePositives = sum(choiceSet['choice']==1)
    truePositivesPlusFalseNegatives = sum(alts['choice']==1)
    coverage = truePositives/truePositivesPlusFalseNegatives
    trueNegativesPlusFalsePositives = sum(alts['choice']==0)
    trueNegatives = sum((~alts.index.isin(choiceSet.index)) & (alts['choice']==0))
    efficiency = trueNegatives/trueNegativesPlusFalsePositives
    indicator = np.abs(coverage - efficiency)
    return(indicator)


def eba(multiplier):
    # Data Prep------------------------------------------------------------------------------
    alts = getAllAlts()
#    weirdPeople = 
#    alts = alts.loc[(alts['distRatio']<4)&(alts['choice']==1)])
    # Settings-----------------------------------------------------------------------
#    attributeTypes = ['distRatio', 'ttRatio', 'transfersDiff']
#    thresholdRange = [np.array(range(1000,10010,10))/1000,np.array(range(1000,10010,10))/1000,
#                      np.array(range(min(alts['transfersDiff']),max(alts['transfersDiff'])+1))]
    attributeTypes = ['distDiff', 'ttDiff', 'transfersDiff']
    thresholdRange = [np.array(range(500,2000)),
                      np.array(range(100,1500)),
                      np.array(range(int(min(alts[attributeTypes[2]])),int(max(alts[attributeTypes[2]])+1)))]
#    thresholdRange = [np.array(range(int(min(alts[attributeTypes[0]])),int(max(alts[attributeTypes[0]])+1))),
#                      np.array(range(int(min(alts[attributeTypes[1]])),int(max(alts[attributeTypes[1]])+1))),
#                      np.array(range(int(min(alts[attributeTypes[2]])),int(max(alts[attributeTypes[2]])+1)))]
    # Indicator calculation---------------------------------------------------------
    attributeCombination =  list(itertools.permutations(range(len(attributeTypes)))) # all possible combinations of attributes
    index = pd.MultiIndex.from_tuples(attributeCombination)
    indicatorStore = pd.DataFrame(columns=attributeTypes,index=index)
    minIndicatorStore = pd.DataFrame(columns=attributeTypes,index=index)
    truePositives = pd.DataFrame(columns=attributeTypes,index=index)
    truePositivesPlusFalseNegatives = pd.DataFrame(columns=attributeTypes,index=index)
    trueNegativesPlusFalsePositives = pd.DataFrame(columns=attributeTypes,index=index)
    trueNegatives = pd.DataFrame(columns=attributeTypes,index=index)
    
    for attributeSet in attributeCombination:
        print(attributeSet)
        prevAtt = []
        minIndicator = []
        for att in attributeSet:
            print(attributeTypes[att])
            temp_size = (len(thresholdRange[att]),)
            temp_indicator = np.zeros(temp_size)
            temp_tp = np.zeros(temp_size)
            temp_tpfn = np.zeros(temp_size)
            temp_tnfp = np.zeros(temp_size)
            temp_tn = np.zeros(temp_size)
            
            altSet = alts
            if prevAtt: # filter according to thresholds of previous attributes
                altSet = thresholdFilter(altSet,prevAtt,minIndicator,attributeTypes)
            if len(altSet)==0:
                continue
            
            for iThreshold in range(len(thresholdRange[att])):
                threshold = thresholdRange[att][iThreshold]
                (temp_tp[iThreshold],temp_tpfn[iThreshold],
                 temp_tnfp[iThreshold],temp_tn[iThreshold]) = indicatorCalc(altSet,threshold,att,attributeTypes)
            
            truePositives.loc[attributeSet,attributeTypes[att]] = temp_tp
            truePositivesPlusFalseNegatives.loc[attributeSet,attributeTypes[att]] = temp_tpfn
            trueNegativesPlusFalsePositives.loc[attributeSet,attributeTypes[att]] = temp_tnfp
            trueNegatives.loc[attributeSet,attributeTypes[att]] = temp_tn
            
            temp_coverage = temp_tp/temp_tpfn
            temp_efficiency = temp_tn/temp_tnfp
            temp_indicator = temp_coverage/multiplier - temp_efficiency
            indicatorStore.loc[attributeSet,attributeTypes[att]] = temp_indicator
            
            temp_minIndicator = np.argmin(np.absolute(temp_indicator),axis=0)
            
            prevAtt += [att]
            if len([temp_minIndicator])>1:
                minIndicator += [thresholdRange[att][temp_minIndicator[0]]]
                print(attributeTypes[att]+' '+attributeSet+' '+len(temp_minIndicator))
            else:
                minIndicator += [thresholdRange[att][temp_minIndicator]]
        
        for i in range(len(attributeSet)):
            minIndicatorStore.loc[attributeSet][attributeSet[i]] = minIndicator[i]
    
    # Final indicators -------------------------------------------------------------
    attributeCombination =  list(itertools.permutations(range(len(attributeTypes)))) # all possible combinations of attributes
    index = pd.MultiIndex.from_tuples(attributeCombination)
    smallestIndicators = pd.DataFrame(columns=attributeTypes,index=index)
    for attributeSet in attributeCombination:
        for att in attributeTypes:
            smallestIndicators.loc[attributeSet,att] = np.nanmin(
                    np.absolute(indicatorStore.loc[attributeSet,att]))
    
    finalPerformance = pd.DataFrame(data=np.log(np.prod(smallestIndicators.values,axis=1).tolist()),
                                    index=smallestIndicators.index).sort_values([0])
    
    finalIndicatorStore = pd.DataFrame(columns=['indicator'],index=index)
    finalWithObs = pd.DataFrame(columns=['withObs'],index=index)
    for i in attributeCombination:
        thresholds = minIndicatorStore.loc[i]
        choiceSet = thresholdFilter2(alts,i,thresholds.tolist(),attributeTypes)
        finalIndicatorStore.loc[i] = finalIndicatorCalc(alts,choiceSet)
        finalWithObs.loc[i] = len(set(choiceSet.loc[choiceSet['choice']==1]['persId']))
        
    finalIndicatorStore = finalIndicatorStore.sort_values('indicator')
    finalWithObs = finalWithObs.sort_values('withObs')
#    ranking = (1,2,0)
    
    # Ranking selection-----------------------------------------------------------------------
    ranking = finalPerformance.idxmin().values[0]
    thresholds = minIndicatorStore.loc[ranking]
    
    # Choice set analysis---------------------------------------------------------------------
    
    choiceSet = thresholdFilter2(alts,list(ranking),thresholds.tolist(),attributeTypes)
    withObs = set(choiceSet.loc[choiceSet['choice']==1]['persId'])
    choiceSetSize = choiceSet.groupby(['persId']).size()
    withChoice = set(choiceSetSize.loc[choiceSetSize>1].index)
    withoutChoice = set(choiceSetSize.loc[choiceSetSize==1].index)
    finalChoiceSet = choiceSet.loc[choiceSet['persId'].isin(withObs & withChoice)]
    finalChoiceSetSize = finalChoiceSet.groupby(['persId']).size()
    a = pd.read_excel('../Data/walking-cyclingpeople.xlsx')
    withBicycle = set(a.loc[a['Bicycle']==1]['persId'])
    
    len(withChoice & withObs & withBicycle)
    len(withObs)
    len(withChoice)
    len(withoutChoice)
    len(withChoice & withObs)
    len(withoutChoice & withObs)
    # Return---------------------------------------------------------------------------------
    return(finalPerformance,finalIndicatorStore,finalWithObs,minIndicatorStore)


#%% Main
if __name__=='__main__':
    multiplier = 1
    
    m = []
    finalPerformance = []
    finalIndicatorStore = []
    finalWithObs = []
    minIndicatorStore = []
    for i in range(0,11):
        multiplier = 1+(i*0.1)
        print(multiplier)
        (a,b,c,d) = eba(multiplier)
        m += [multiplier]
        finalPerformance += [a]
        finalIndicatorStore += [b]
        finalWithObs += [c]
        minIndicatorStore += [d]
        
#
#    names = ['multiplier','hasObs','canChoose','isEligible','ranking',
#             'distRatio','ttRatio','transferDiff']
#    sensitivity = pd.DataFrame(data=[m,withObs,withChoice,eligiblePers,ranking,dist,tt,tra])
#    sensitivity = sensitivity.transpose()
#    sensitivity.columns = names
#    sensitivity.to_hdf('../vars/sensitivity.h5',key='df',mode='w')
#    sensitivity.to_excel('../vars/sensitivity.xlsx')

#%% Export final choice set for Danique
alts['eligibleAlt'] = 0
alts.loc[finalChoiceSet.index,'eligibleAlt'] = 1
alts['eligiblePerson'] = alts['persId'].isin(withObs & withChoice).astype(int)

a = pd.read_excel('../Data/walking-cyclingpeople.xlsx')
a1 = a.loc[a['Bicycle']==1]['persId']
alts['usesBicycle'] = 0
alts.loc[alts['persId'].isin(a1),'usesBicycle'] = 1

alts.to_excel('../Data/finalAlts5.xlsx')

finalChoiceSet.to_excel('../Data/justChoiceSet5.xlsx')

#%% Histogram of max acceptable access distances according to EBA analysis
a2 = finalChoiceSet.groupby(['persId'])['dist'].min()+thresholds.iloc[0]
print(a2.median())
print(a2.quantile(0.9))
import matplotlib.pyplot as plt
plt.hist(a2.values,bins=30,color='orange')
plt.xlabel('Maximum access distance (m)')
plt.ylabel('Frequency')
plt.savefig('../figures/minPlusThresholdAccessDist.png', bbox_inches='tight',dpi=300)

a3 = alts.loc[alts['persId'].isin(alts.loc[alts['choice']==1]['persId'])].groupby('persId').size()
print(a3.median())
print(a3.quantile(0.9))
plt.hist(a3.values,bins=np.arange(0-0.5,20+0.5,1),color='blue')
plt.xlabel('Choice set size')
plt.ylabel('Frequency')
plt.savefig('../figures/cssMaster.png', bbox_inches='tight',dpi=300)

a4 = finalChoiceSet.groupby(['persId']).size()
print(a4.median())
print(a4.quantile(0.9))
a4.name = 'css'
plt.hist(a4.values,bins=np.arange(0-0.5,a4.max()+0.5,1),color='blue')
plt.xlabel('Choice set size')
plt.ylabel('Frequency')
plt.savefig('../figures/cssEBA.png', bbox_inches='tight',dpi=300)

a5 = finalChoiceSet.groupby(['persId'])['dist'].max()
print(a5.median())
print(a5.quantile(0.9))
plt.hist(a5.values,bins=30,color='red')
plt.xlabel('Maximum access distance (m)')
plt.ylabel('Frequency')
plt.savefig('../figures/maxInCSAccessDist.png', bbox_inches='tight',dpi=300)

a6 = finalChoiceSet.loc[finalChoiceSet['choice']==1]['dist']
print(a6.median())
print(a6.quantile(0.9))
plt.hist(a6.values,bins=30,color='cyan')
plt.xlabel('Maximum access distance (m)')
plt.ylabel('Frequency')
plt.savefig('../figures/csAccessDist.png', bbox_inches='tight',dpi=300)


#%% Plot choice set size by location
csSize = alts.drop_duplicates(['persId','lat-home','lon-home'])[['persId','lat-home','lon-home']].reset_index(drop=True)
csSize = pd.merge(csSize,a4,left_on='persId',right_index=True)
csSize.to_csv('../GIS/csSize1.csv')

#%%
cond = ((alts['distRatio']<=a['distRatio'].max()) &
 (alts['ttRatio']<=a['ttRatio'].max()) &
 (alts['transfersDiff']<=a['transfersDiff'].max()))

#%% Match tram type to lines
tramTypes = pd.read_excel('../Data/tramTypes.xlsx')
finalChoiceSet['lines2'] = np.array([[k for j in i for k in j] for i in finalChoiceSet['lines']])
finalChoiceSet['tramTypeGTL'] = [len(set(i)&set(tramTypes.loc[tramTypes['tramType']==1,'tramLine']))>0 for i in finalChoiceSet['lines2']]
finalChoiceSet['tramTypeAvenio'] = [len(set(i)&set(tramTypes.loc[tramTypes['tramType']==2,'tramLine']))>0 for i in finalChoiceSet['lines2']]
finalChoiceSet['tramTypeRand'] = [len(set(i)&set(tramTypes.loc[tramTypes['tramType']==3,'tramLine']))>0 for i in finalChoiceSet['lines2']]
finalChoiceSet.to_excel('../Data/justChoiceSet5_withTramTypes.xlsx')







