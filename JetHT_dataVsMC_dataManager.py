#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 15:42:04 2021

@author: si_sutantawibul1
"""

import pickle
#from EffInfoABC import EfficiencyInfoABC,EfficiencyInfoAB,EfficiencyInfoC
from ScaleFactor import ScaleFactor
import uproot
import pandas as pd
from analib import PhysObj, Event
import sys
import os
import numpy as np
from dataManager import getSubJetData, getnSVCounts, trainVars
from htoaaRootFilesLoc import TTJetsPaths, WJetsPaths, bEnrPaths, BGenPaths
from htoaaRootFilesLoc import ZJetsPaths, ParkedDataPaths, JetHTPaths, ggHPaths, MuonEGPaths
import copy as cp
#from dataVsMC_DataManager import getdR, getdRCount, getMaxPt, getNthPt, BGenDict, bEnrDict, ptkeys, ipkeys, npvsGkeys, ParkedDataDict
import dataVsMC_DataManager as DM
import matplotlib.pyplot as plt
import json

np.seterr(divide='ignore', invalid='ignore')
np.set_printoptions(suppress=True)

trigTensor = pickle.load(open('TrigTensor.p','rb'))
pathlist = ['ABC', 'AB', 'C']
taglist = ['BGen', 'bEnr', 'Parked', 'ggH', 'JetHT', 'ZJets', 'WJets', 'TTJets', 'QCDInc']

plotVars = ['FatJet_pt', 'FatJet_eta', 'FatJet_mass', 'FatJet_btagCSVV2', 'FatJet_btagDeepB',
            'FatJet_msoftdrop', 'FatJet_btagDDBvL', 'FatJet_deepTagMD_H4qvsQCD', 'FatJet_n2b1',
            'SubJet_mass(1)','SubJet_mass(2)','SubJet_tau1(1)','FatJet_nSV', 'PV_npvs', 'PV_npvsGood',
            'BDTScore', 'LHE_HT', 'FatJet_nSV']


def getMaxPtDf(filepath, ev, MC, path, tag, events):
    jets = ev.objs['jets']
    other = ev.objs['other']
    if jets.FatJet_pt.empty:
        return pd.DataFrame()
    maxPtData = DM.getMaxPt(jets, 'FatJet_pt')
    if MC:
        maxPtData = maxPtData.assign(LHE_HT=other['LHE_HT'].to_numpy())
    maxPtData = maxPtData.assign(PV_npvs=other['PV_npvs'].to_numpy())
    maxPtData = maxPtData.assign(PV_npvsGood=other['PV_npvsGood'].to_numpy())
    if 'AB' != path:
        ak4Jets = ev.objs['ak4Jets']
        maxPtData = maxPtData.assign(Jet_pt1=DM.getNthPt(n=1, physobj=ak4Jets,
                                                  sortBy='Jet_pt',
                                                  extractCol='Jet_pt'))
        maxPtData = maxPtData.assign(Jet_pt2=DM.getNthPt(n=2, physobj=ak4Jets,
                                                  sortBy='Jet_pt',
                                                  extractCol='Jet_pt'))
        maxPtData = maxPtData.assign(Jet_btagDeepB1=DM.getNthPt(n=1,
                                                         physobj=ak4Jets,
                                                         sortBy='Jet_pt',
                                                         extractCol='Jet_btagDeepB'))
        maxPtData = maxPtData.assign(Jet_btagDeepB2=DM.getNthPt(n=2,
                                                         physobj=ak4Jets,
                                                         sortBy='Jet_pt',
                                                         extractCol='Jet_btagDeepB'))
        btags = pd.DataFrame(maxPtData.Jet_btagDeepB1, columns=['Jet_btagDeepB1'])
        btags = btags.assign(Jet_btagDeepB2=maxPtData.Jet_btagDeepB2)
        maxPtData = maxPtData.assign(Jet_btagDeepB1=btags.max(axis=1))
        maxPtData = maxPtData.assign(Jet_btagDeepB2=btags.min(axis=1))


    maxPtData['final_weights'] = 1
    if 'ggH'==tag:
        maxPtData['LHE_weights'] = 0.0046788#1
        wgt = 3.9 - 0.4*np.log2(maxPtData.FatJet_pt)
        wgt[wgt<0.1] = 0.1
        maxPtData['ggH_weights'] = wgt
        maxPtData['final_weights'] = (maxPtData['LHE_weights'] *
                                      maxPtData['ggH_weights'])
    elif 'BGen'==tag:
        maxPtData['LHE_weights'] = DM.BGenDict[filepath]
        # wgt = 4.346 - 0.356*np.log2(maxPtData.LHE_HT)
        # wgt[wgt<0.1] = 0.1
        # maxPtData['QCD_correction'] = wgt
        # Xsec_wgt = 10.78
        maxPtData = maxPtData.assign(final_weights =
                                     maxPtData['LHE_weights'])
                                     #maxPtData['QCD_correction']
                                     #*Xsec_wgt)
    elif 'bEnr' == tag:
        maxPtData['LHE_weights'] = DM.bEnrDict[filepath]
        # wgt = 4.346 - 0.356*np.log2(maxPtData.LHE_HT)
        # wgt[wgt<0.1] = 0.1
        # maxPtData['QCD_correction'] = wgt
        # Xsec_wgt = 4.1
        maxPtData = maxPtData.assign(final_weights=
                                     maxPtData['LHE_weights'])
                                     #maxPtData['QCD_correction']
                                     #*Xsec_wgt)
    elif 'QCDInc'==tag:
        maxPtData['LHE_weights'] =DM.QCDIncDict[filepath]
        maxPtData = maxPtData.assign(final_weights=
                                     maxPtData['LHE_weights'])
    elif 'WJets'==tag:
        maxPtData['LHE_weights'] = DM.WJetsDict[filepath]
        maxPtData = maxPtData.assign(final_weights = maxPtData['LHE_weights'])
    elif 'ZJets'==tag:
        maxPtData['LHE_weights'] = DM.ZJetsDict[filepath]
        maxPtData = maxPtData.assign(final_weights = maxPtData['LHE_weights'])
    elif 'TTJets'==tag:
        maxPtData['LHE_weights'] = DM.TTJetsDict[filepath]
        maxPtData = maxPtData.assign(final_weights=maxPtData['LHE_weights'])
        ## added this for checks
        #wgt = 4.346 - 0.356*np.log2(maxPtData.LHE_HT)
        #wgt[wgt<0.1] = 0.1
        #maxPtData = maxPtData.assign(final_weights=maxPtData['LHE_weights']*wgt)

    if MC:
        #bins = pd.IntervalIndex.from_breaks(breaks=trigTensor['meta'][path], closed='left')
        if 'C'==path:
            bins = trigTensor['meta']['CX']
            ##!!! TODO: swap ak4 jets values here too don't forget to do that
            trigWgt = pd.cut(x=maxPtData.Jet_btagDeepB2, bins=bins, right=False, include_lowest=True,
                             labels=trigTensor['CX']).astype(float)
        else:
            bins = trigTensor['meta'][path]
            trigWgt = pd.cut(x=maxPtData.FatJet_pt, bins=bins, right=False, include_lowest=True,
                             labels=trigTensor[path]).astype(float)
        maxPtData = maxPtData.assign(trigWeight = trigWgt)
        maxPtData = maxPtData.assign(final_weights = maxPtData['final_weights']*
                                     maxPtData['trigWeight']*54.54)
        ## added ask bhoff if this is for all of just qcd
        # wgt = 4.346 - 0.356*np.log2(maxPtData.LHE_HT)
        # wgt[wgt<0.1] = 0.1
        # maxPtData['QCD_correction'] = wgt
        # maxPtData = maxPtData.assign(final_weights=maxPtData['final_weights']*
        #                              maxPtData['QCD_correction'])


    # else:
    #     maxPtData['event'] = other['event'].values.astype(int)
    #     maxPtData['run'] = other['run'].values.astype(int)
    #     maxPtData['luminosityBlock'] = other['luminosityBlock'].values.astype(int)

    maxPtData['event'] = other['event'].values.astype(int)
    maxPtData['run'] = other['run'].values.astype(int)
    maxPtData['luminosityBlock'] = other['luminosityBlock'].values.astype(int)
    maxPtData['FatJet_nSV'] = getnSVCounts(jets, events)
    maxPtData['FatJet_eta'] = maxPtData['FatJet_eta'].abs()
    #maxPtData['EventNum'] = jets.FatJet_pt.index
    return maxPtData


def process(filepath, MC, tag):
    print(filepath)

    if type(MC) != bool:
        print('MC needs to be set true/false')
        sys.exit()

    if tag not in taglist:
        print('check tags')
        sys.exit()

    fileName, fileExtension = os.path.splitext(filepath)
    if fileExtension != '.root':
        print('this program only supports .root  files')
        sys.exit()

    jets = PhysObj('jets')
    other = PhysObj('other')
    ak4Jets = PhysObj('ak4Jets')
    trig = PhysObj('trig')
    muon = PhysObj('muons')

    f = uproot.open(fileName + '.root')
    events = f.get('Events')

    jets['FatJet_pt'] = pd.DataFrame(events.array('FatJet_pt'))
    jets['FatJet_eta'] = pd.DataFrame(events.array('FatJet_eta'))
    jets['FatJet_mass'] = pd.DataFrame(events.array('FatJet_mass'))
    jets['FatJet_btagCSVV2'] = pd.DataFrame(events.array('FatJet_btagCSVV2'))
    jets['FatJet_btagDeepB'] = pd.DataFrame(events.array('FatJet_btagDeepB'))
    jets['FatJet_msoftdrop'] = pd.DataFrame(events.array('FatJet_msoftdrop'))
    jets['FatJet_btagDDBvL'] = pd.DataFrame(events.array('FatJet_btagDDBvL'))
    jets['FatJet_deepTagMD_H4qvsQCD'] = pd.DataFrame(events.array('FatJet_deepTagMD_H4qvsQCD'))
    jets['FatJet_phi'] = pd.DataFrame(events.array('FatJet_phi'))
    jets['FatJet_n2b1'] = pd.DataFrame(events.array('FatJet_n2b1'))
    jets['SubJet_mass(1)'] = getSubJetData(1,'SubJet_mass', events)
    jets['SubJet_mass(2)'] = getSubJetData(2, 'SubJet_mass', events)
    jets['SubJet_tau1(1)'] = getSubJetData(1, 'SubJet_tau1', events)

    ak4Jets['Jet_pt'] = pd.DataFrame(events.array('Jet_pt'))
    ak4Jets['Jet_eta'] = pd.DataFrame(events.array('Jet_eta'))
    ak4Jets['Jet_puId'] = pd.DataFrame(events.array('Jet_puId'))
    ak4Jets['Jet_phi'] = pd.DataFrame(events.array('Jet_phi'))
    ak4Jets['Jet_btagDeepB'] = pd.DataFrame(events.array('Jet_btagDeepB'))

    if MC:
        other['LHE_HT'] = pd.DataFrame(events.array('LHE_HT')).astype(np.float64)
    other['PV_npvs'] = pd.DataFrame(events.array('PV_npvs'))
    other['PV_npvsGood'] = pd.DataFrame(events.array('PV_npvsGood'))
    other['run'] = pd.DataFrame(events.array('run').astype(int))
    other['event'] = pd.DataFrame(events.array('event').astype(int))
    other['luminosityBlock'] = pd.DataFrame(events.array('luminosityBlock').astype(int))

    

    ## triggers
    trigABC1 = events.array('L1_SingleJet180') & (events.array('HLT_AK8PFJet500') |
                                               events.array('HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4'))
    trigABC2 = (events.array('L1_DoubleJet112er2p3_dEta_Max1p6') |
             events.array('L1_DoubleJet150er2p5')) & events.array('HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71')
    trig['trigABC'] = pd.DataFrame((trigABC1 | trigABC2).astype(bool))

    trigAB = events.array('L1_SingleJet180') & (events.array('HLT_AK8PFJet500') |
                                                events.array('HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4'))
    trig['trigAB'] = pd.DataFrame(trigAB.astype(bool))

    trigC = (events.array('L1_DoubleJet112er2p3_dEta_Max1p6') |
             events.array('L1_DoubleJet150er2p5')) & events.array('HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71')
    trig['trigC'] = pd.DataFrame(trigC.astype(bool))


    ev = Event(jets, other, trig)
    
    
    
    
    ii =104512
     
    for obj in ev.objs:
        for key in ev.objs[obj]:
            print(key)
            print(ev.objs[obj][key].loc[ii])
            print('--------------------------')
    
    

    ## select base on genparticles 
    if MC:
        gen = PhysObj('gen')
        gen['pdgId'] = pd.DataFrame(events.array('GenPart_pdgId'))
        gen['gen_pt'] = pd.DataFrame(events.array('GenPart_pt'))
        gen.cut((gen.pdgId==5)&(gen.gen_pt>15))
        if ('BGen'==tag) or ('bEnr'==tag):
            ev.register(gen)
            ev.sync()
            ev.objs.pop('gen')
        elif 'QCDInc'==tag:
            gen.cut((gen.pdgId==5)&(gen.gen_pt>15))
            jets.FatJet_pt.drop(gen.pdgId.index, inplace=True)
            ev.sync()

    jets.cut(jets['FatJet_pt'] >= 250) #250)
    jets.cut(jets['FatJet_eta'].abs() < 2.4)
    jets.cut(jets['FatJet_btagDDBvL'] > 0.8)
    jets.cut(jets['FatJet_btagDeepB'] > 0.4184)
    jets.cut(jets['FatJet_msoftdrop'] > 90)
    jets.cut(jets['FatJet_msoftdrop'] < 200)
    jets.cut(jets['FatJet_mass'] > 90)
    jets.cut(jets['FatJet_mass'] < 200)
    other.cut(other['PV_npvsGood'] >= 1)


    # ak4Jets.cut(ak4Jets['Jet_pt']  > 30)
    # ak4Jets.cut(ak4Jets['Jet_eta'].abs() < 2.4)
    # ak4Jets.cut(ak4Jets['Jet_puId'] >= 1)

    ev.sync()

    ## golden json cuts
    if False == MC:

        print('we in the GJSON selection')

        jdata = json.load(open('C2018.json'))
        dtev = PhysObj('event')
        dtev.run = pd.DataFrame(events.array('run'))
        dtev.lb = pd.DataFrame(events.array('luminosityBlock'))
        ev.register(dtev)
        ev.sync()
        ## Only events that are in keys are kept
        dtev.cut(dtev.run.isin(jdata.keys())==True)
        for elem in dtev:
            dtev[elem]=dtev[elem].astype(int)
        def fun(r,lb):
            return any([lb in range(a[0],a[1]+1) for a in jdata[str(r)] ])
        
        truthframe = pd.DataFrame([fun(r,lb) for r,lb in zip(dtev.run[0].values,dtev.lb[0].values)],index=dtev.run.index)
        dtev.cut(truthframe == True)
        ev.sync()


    ev.sync()

    # make C events
    Cjets = jets.deepcopy()
    Cak4Jets = ak4Jets.deepcopy()
    Cother = other.deepcopy()
    Ctrig = trig.deepcopy()

    Cev = Event(Cjets, Cother, Cak4Jets, Ctrig)

    #process further for ABC/AB
    jets.cut(jets['FatJet_pt'] >= 400)

    #make ABC/AB events
    ABjets = jets.deepcopy()
    ABak4Jets = ak4Jets.deepcopy()
    ABother = other.deepcopy()
    ABtrig = trig.deepcopy()
    ABev = Event(ABjets, ABother, ABtrig)

    ABCjets = jets
    ABCak4Jets = ak4Jets
    ABCother = other
    ABCtrig = trig
    ABCev = Event(ABCjets, ABCak4Jets, ABCother, ABCtrig)

    del jets, ak4Jets, other, trig

    ## perform ABC cuts
    ABCak4Jets.cut(ABCak4Jets['Jet_btagDeepB'] > 0.4184)
    ABCak4Jets.cut(ABCak4Jets['Jet_pt'] > 140)
    ABCak4Jets.cut(ABCak4Jets['Jet_eta'].abs() < 2.4)
    ABCak4Jets.cut(ABCak4Jets['Jet_puId'] >= 1)
    ABCev.sync()
    ABCak4Jets['dR'] = DM.getdR(objName='Jet', events=events, fatJetPhysObj=ABCjets, jetPhysObj=ABCak4Jets)
    ABCother['JetFatJet_dRCnt'] = DM.getdRCount(ABCak4Jets['dR'])
    ABCak4Jets.cut(ABCak4Jets['dR'] < 0.8)
    ABCother.cut(ABCother['JetFatJet_dRCnt'] >= 2)
    # ABCtrig.cut(ABCtrig['trigABC']==True)
    ABCev.sync()

    ## perform AB selection (drop events in ABC from B)
    ABjets.FatJet_pt.drop(ABCjets.FatJet_pt.index, inplace=True)
    ABtrig.cut(ABtrig['trigAB']==True)
    ABev.sync()

    ##since AB events that passes ABC trigger can survive, do ABCtrigger after
    ## event removal from AB
    ABCtrig.cut(ABCtrig['trigABC']==True)
    ABCev.sync()

    ## perform C cuts
    ## cut all but the fattest jet

    Cjets.cut(Cjets.FatJet_pt.rank(method='max', axis=1, ascending=False)==1)
    Cjets.cut(Cjets['FatJet_pt'] < 400)
    Cak4Jets.cut(Cak4Jets['Jet_btagDeepB'] > 0.4184)
    Cak4Jets.cut(Cak4Jets['Jet_pt'] > 140)
    Cak4Jets.cut(Cak4Jets['Jet_eta'].abs() < 2.4)
    Cak4Jets.cut(Cak4Jets['Jet_puId'] >= 1)
    Cjets.FatJet_pt.drop(ABCjets.FatJet_pt.index, inplace=True, errors='ignore')
    Cjets.FatJet_pt.drop(ABjets.FatJet_pt.index, inplace=True, errors='ignore')

    Cev.sync()

    Cak4Jets['dR'] = DM.getdR(objName='Jet', events=events, fatJetPhysObj=Cjets, jetPhysObj=Cak4Jets)
    Cother['JetFatJet_dRCnt'] = pd.DataFrame((Cak4Jets.dR<0.8).sum(axis=1))
    Cak4Jets.cut(Cak4Jets['dR'] < 0.8)
    Cother.cut(Cother['JetFatJet_dRCnt'] >= 2)
    Ctrig.cut(Ctrig['trigC']==True)
    Cev.sync()

    CXdf = pd.DataFrame(Cother['run'])
    # CXdf = CXdf.rename({0:'run'})
    CXdf = CXdf.assign(luminosityBlock = Cother['luminosityBlock'])
    CXdf = CXdf.assign(event = Cother['event'])

    ABCdf = pd.DataFrame(ABCother['run'])
    # ABCdf = ABCdf.rename({0:'run'})
    ABCdf = ABCdf.assign(luminosityBlock = ABCother['luminosityBlock'])
    ABCdf = ABCdf.assign(event = ABCother['event'])

    ABdf = pd.DataFrame(ABother['run'])
    #ABdf = ABdf.rename({0:'run'})
    ABdf = ABdf.assign(luminosityBlock = ABother['luminosityBlock'])
    ABdf = ABdf.assign(event = ABother['event'])

    pickle.dump(CXdf, open('JetHTTrigEff/frames/CXevents.pkl','wb'))
    pickle.dump(ABCdf, open('JetHTTrigEff/frames/ABCevents.pkl','wb'))
    pickle.dump(ABdf, open('JetHTTrigEff/frames/ABevents.pkl','wb'))

    ABCDf = getMaxPtDf(filepath=filepath, ev=ABCev, MC=MC, path='ABC', tag=tag, events=events)
    ABDf = getMaxPtDf(filepath=filepath, ev=ABev, MC=MC, path='AB', tag=tag, events=events)
    CDf = getMaxPtDf(filepath=filepath, ev=Cev, MC=MC, path='C', tag=tag, events=events)




    return {'ABC': ABCDf, 'AB': ABDf, 'C': CDf}



#%%
## function to get center of ins given binedges as np array
def getBinCenter(arr):
    arrCen = list()
    for i in range(len(arr)-1):
        arrCen.append((arr[i+1]+arr[i])/2)
    return arrCen

## function to add a BDTScore column to each of the background/signal/data DF
loadedModel = pickle.load(open('BDTModels/Htoaa_BDThigh disc.pkl', 'rb'))
def analyze(dataDf):
    prediction = loadedModel.predict_proba(dataDf[trainVars])
    dataDf = dataDf.assign(BDTScore=prediction[:,1])
    return dataDf


#======================== main ==============================
#def makefinaldf(ev, MC, tag, path):


ranges = {'FatJet_pt': (250,850),
          'FatJet_eta': (0,3),
          'FatJet_phi': (-3.2, 3.2),
          'FatJet_mass': (85, 205),
          'FatJet_btagCSVV2': (0,1.1),
          'FatJet_btagDeepB': (0.35, 1.05),
          'FatJet_msoftdrop': (85,205),
          'FatJet_btagDDBvL': (0.76, 1.02),
          'FatJet_deepTagMD_H4qvsQCD': (0,1),
          'PV_npvs': (0,80),
          'PV_npvsGood': (0,80),
          'LHE_HT': (0, 5000),
          'FatJet_n2b1': (0, 1),
          'SubJet_mass(1)': (-5,105),
          'SubJet_mass(2)': (-5,105),
          'SubJet_tau1(1)': (0,.5),
          'FatJet_nSV': (-1,11),
          'BDTScore': (0,1),
          }

nbins = {'FatJet_pt': 30,
         'FatJet_eta': 15,
         'FatJet_phi': 32,
         'FatJet_mass': 32,
         'FatJet_btagCSVV2': 22,
         'FatJet_btagDeepB': 14,
         'FatJet_msoftdrop': 32,
         'FatJet_btagDDBvL': 26,
         'FatJet_deepTagMD_H4qvsQCD': 20,
         'PV_npvs': 40,
         'PV_npvsGood': 40,
         'LHE_HT': 500,
         'FatJet_n2b1': 20,
         'SubJet_mass(1)': 22,
         'SubJet_mass(2)': 22,
         'SubJet_tau1(1)': 20,
         'FatJet_nSV': 12,
         'BDTScore': 10,
          }



pickledir = 'JetHTTrigEff/dataVsMC/pickles'
append_params = {'ignore_index':True, 'sort':False}
root=True
#%%
if root:
    ggH = process(filepath=DM.ggHPaths, MC=True, tag='ggH')
    for key in ggH:
        ggH[key] = analyze(ggH[key])
    pickle.dump(ggH, open(f'{pickledir}/ggH.pkl','wb'))
#%%
    BGen = {'ABC': pd.DataFrame(), 'AB': pd.DataFrame(), 'C': pd.DataFrame()}
    for fileName in DM.BGenPaths:
        tmp = process(filepath=fileName, MC=True, tag='BGen')
        BGen['ABC'] = BGen['ABC'].append(tmp['ABC'], **append_params)
        BGen['AB'] = BGen['AB'].append(tmp['AB'], **append_params)
        BGen['C'] = BGen['C'].append(tmp['C'], **append_params)
    for key in BGen:
        BGen[key] = analyze(BGen[key])
    pickle.dump(BGen, open(f'{pickledir}/BGen.pkl','wb'))

    bEnr = {'ABC': pd.DataFrame(), 'AB': pd.DataFrame(), 'C': pd.DataFrame()}
    for fileName in DM.bEnrPaths:
        tmp=process(filepath=fileName, MC=True, tag='bEnr')
        bEnr['ABC'] = bEnr['ABC'].append(tmp['ABC'], **append_params)
        bEnr['AB'] = bEnr['AB'].append(tmp['AB'], **append_params)
        bEnr['C'] = bEnr['C'].append(tmp['C'], **append_params)
    for key in bEnr:
        bEnr[key] = analyze(bEnr[key])
    pickle.dump(bEnr, open(f'{pickledir}/bEnr.pkl','wb'))
    #%%
    WJets = {'ABC': pd.DataFrame(), 'AB': pd.DataFrame(), 'C': pd.DataFrame()}
    for fileName in DM.WJetsPaths:
        tmp=process(filepath=fileName, MC=True, tag='WJets')
        WJets['ABC'] = WJets['ABC'].append(tmp['ABC'], **append_params)
        WJets['AB'] = WJets['AB'].append(tmp['AB'], **append_params)
        WJets['C'] = WJets['C'].append(tmp['C'], **append_params)
    for key in WJets:
        WJets[key] = analyze(WJets[key])
    pickle.dump(WJets, open(f'{pickledir}/WJets.pkl','wb'))
    
    ZJets = {'ABC': pd.DataFrame(), 'AB': pd.DataFrame(), 'C': pd.DataFrame()}
    for fileName in DM.ZJetsPaths:
        tmp=process(filepath=fileName, MC=True, tag='ZJets')
        ZJets['ABC'] = ZJets['ABC'].append(tmp['ABC'], **append_params)
        ZJets['AB'] = ZJets['AB'].append(tmp['AB'], **append_params)
        ZJets['C'] = ZJets['C'].append(tmp['C'], **append_params)
    for key in ZJets:
        ZJets[key] = analyze(ZJets[key])
    pickle.dump(ZJets, open(f'{pickledir}/ZJets.pkl','wb'))
#%%
    TTJets = process(filepath=DM.TTJetsPaths[0], MC=True, tag='TTJets')
    for key in TTJets:
        TTJets[key] = analyze(TTJets[key])
    pickle.dump(TTJets, open(f'{pickledir}/TTJets.pkl','wb'))

#%%
    QCDInc = {'ABC': pd.DataFrame(), 'AB': pd.DataFrame(), 'C': pd.DataFrame()}
    for fileName in DM.QCDIncPaths:
        tmp = process(filepath=fileName, MC=True, tag='QCDInc')
        QCDInc['ABC'] = QCDInc['ABC'].append(tmp['ABC'], **append_params)
        QCDInc['AB'] = QCDInc['AB'].append(tmp['AB'], **append_params)
        QCDInc['C'] = QCDInc['C'].append(tmp['C'], **append_params)
    for key in QCDInc:
        QCDInc[key] = analyze(QCDInc[key])
    pickle.dump(QCDInc, open(f'{pickledir}/QCDInc.pkl','wb'))

#%%
    
    JetHT = {'ABC': pd.DataFrame(), 'AB': pd.DataFrame(), 'C': pd.DataFrame()}
    for fileName in DM.JetHTPaths:
        tmp=process(filepath=fileName, MC=False, tag='JetHT')
        JetHT['ABC'] = JetHT['ABC'].append(tmp['ABC'], **append_params)
        JetHT['AB'] = JetHT['AB'].append(tmp['AB'], **append_params)
        JetHT['C'] = JetHT['C'].append(tmp['C'], **append_params)
    for key in JetHT:
        JetHT[key] = analyze(JetHT[key])
    pickle.dump(JetHT, open(f'{pickledir}/JetHT.pkl','wb'))
#%%
else:
    ggH = pickle.load(open(f'{pickledir}/ggH.pkl','rb'))
    BGen = pickle.load(open(f'{pickledir}/BGen.pkl','rb'))
    bEnr = pickle.load(open(f'{pickledir}/bEnr.pkl','rb'))
    WJets = pickle.load(open(f'{pickledir}/WJets.pkl','rb'))
    ZJets = pickle.load(open(f'{pickledir}/ZJets.pkl','rb'))
    TTJets = pickle.load(open(f'{pickledir}/TTJets.pkl','rb'))
    JetHT = pickle.load(open(f'{pickledir}/JetHT.pkl','rb'))
    QCDInc = pickle.load(open(f'{pickledir}/QCDInc.pkl','rb'))

bgs = {'WJets': WJets,
       'ZJets': ZJets,
       'QCDInc' :QCDInc,
       'TTJets': TTJets,
       'BGen': BGen,
       'bEnr': bEnr,
       }

wgt = pd.DataFrame()
for path in pathlist:
    for var in plotVars:
        nbin = nbins[var]
        range_local = ranges[var]

        fig, (ax0, ax1) = plt.subplots(figsize=(15,9),nrows=2, gridspec_kw={'height_ratios': [3, 1]},
                                       sharex=True)
        ax0.set_ylabel('events')
        ax1.set_ylabel('ratio')

        ## prep range for plot + prep bg for plot
        xmin = list()
        xmax = list()
        toplot = pd.DataFrame()
        toplotweights = pd.DataFrame()
        toplotlabel = list()
        for key, bgdict in bgs.items():
            bg = bgdict[path]
            toplot = pd.concat([toplot, bg[var]], ignore_index=False, axis=1)
            toplotweights = pd.concat([toplotweights, bg['final_weights']], ignore_index=False,
                      axis=1)
            toplotweights = toplotweights.rename({'final_weights': key})
            toplotlabel.append(f'{key} ({round(np.sum(bg.final_weights))})')

        if 'BDTScore' != var:
            bins = np.linspace(range_local[0], range_local[1], num=nbin)

        density = False
        hist_params = {'density': density, 'histtype': 'bar', 'range' : range_local, 'bins':nbin, 'stacked':True}

        ## making color palette for the QCD stakcs
        pal = ['#603514', '#b940f2','#ba6861', '#ec6f38', '#6acaf8', '#82f759', ]

        ## plotting background MC
        if 'BDTScore'==var:
            bgvals, bgbins, _ = ax0.hist(toplot.values.astype(float),
                                    weights=toplotweights.values.astype(float),
                                    label=toplotlabel,
                                    color=pal, **hist_params)
        else:
            bgvals, bgbins, _ = ax0.hist(np.clip(toplot.values.astype(float), bins[0],bins[-1]),
                                        weights=toplotweights.values.astype(float),
                                        label=toplotlabel,
                                        color=pal, **hist_params)

        ## plotting signal MC
        ## scaling GGH area to equal bg area
        ggHDf = ggH[path]
        scale = np.nansum(toplotweights)/np.sum(ggHDf.final_weights)
        gghweights = ggHDf['final_weights']*scale
        #ggHDf['final_weights'] = gghweights

        if 'BDTScore'==var:
            ax0.hist(ggHDf[var].values, label=f'GGH ({round(np.sum(ggHDf.final_weights))}x{round(scale)})', histtype='step',
                density=density, bins=nbin, linewidth=3, color='r',
                range=range_local, weights=gghweights)
        else:
            ax0.hist(np.clip(ggHDf[var].values, bins[0], bins[-1]), label=f'GGH ({round(np.sum(ggHDf.final_weights))}x{round(scale)})', histtype='step',
                     density=density, bins=nbin, linewidth=3, color='r',
                     range=range_local, weights=gghweights)

        # plotting JetHT
        JetHTDf = JetHT[path]
        x = getBinCenter(bgbins)
        if 'BDTScore'==var:
            # datavals, databins, _ = ax0.hist(JetHTDf[var].values,
            #                              label=f'JetHT ({round(np.sum(JetHTDf.final_weights))})', histtype='step',
            #                              density=density, bins=8, linewidth=3, color='k',
            #                              range=(0,0.8), weights=JetHTDf.final_weights)
            datavals, databins = np.histogram(JetHTDf[var].values, 
                                              density=density, bins=8, range=(0,0.8),
                                              weights=JetHTDf.final_weights.values)
            datavals = np.append(datavals, [0,0])
        elif 'LHE_HT' == var:
            continue
        else:
            # datavals, databins, _ = ax0.hist(np.clip(JetHTDf[var].values, bins[0], bins[-1]),
            #                              label=f'JetHT ({round(np.sum(JetHTDf.final_weights))})', histtype='step',
            #                              density=density, bins=nbin, linewidth=3, color='k',
            #                              range=range_local, weights=JetHTDf.final_weights)
            datavals, databins = np.histogram(np.clip(JetHTDf[var].values, bins[0], bins[-1]), 
                                              density=density, bins=nbin, range=range_local,
                                              weights=JetHTDf.final_weights.values)
        xerr = (x[1]-x[0])/2
        yerr = 1/np.sqrt(datavals)
        ax0.errorbar(x, datavals, xerr=xerr, yerr=yerr, fmt='ko')

        ## plotting ratio
        totalbgvals = bgvals[-1]
        y = np.divide((datavals-totalbgvals), totalbgvals, out = np.zeros_like(totalbgvals),
                      where=(totalbgvals>0))
        #x = getBinCenter(bgbins)
        ax1.grid()
        ymax = 0.5
        ax1.set_ylim(bottom=-ymax, top=ymax)
        #ax1.scatter(x,y, color='k')
        ax1.errorbar(x,y, xerr=xerr, yerr=yerr, fmt='ko')
        #ax1.set_xticklabels(xlabels)


        ax0.set_title(f'{var} {path}' )# + ' JetHT')
        ax0.legend(loc='best')#, frameon=True,bbox_to_anchor=(1, 0.5))
        ax0.grid()

        plt.savefig(f'JetHTTrigEff/dataVsMC/plots/{var}_{path}'  , bbox_inches = "tight")
        plt.show()
        plt.close