#htoaa analysis main code

import os
import sys
import json
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
import numpy as np
from copy import copy, deepcopy
#import uproot
import uproot3 as uproot

'''
H->aa->4b boosted analysis macro

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
* Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''

#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
from coffea.nanoevents.methods import nanoaod, vector
from coffea.analysis_tools import PackedSelection, Weights
#import hist
from coffea import hist
import awkward as ak
import uproot
#from dask.distributed import Client
from particle import Particle # For PDG particle listing https://github.com/scikit-hep/particle



from htoaa_Settings import *
from htoaa_CommonTools import GetDictFromJsonFile, calculate_lumiScale, setXRootDRedirector
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl
)

'''
from htoaa_Settings import *
from htoaa_NanoAODBranches import htoaa_nanoAODBranchesToRead
from htoaa_CommonTools import GetDictFromJsonFile, DfColLabel_convert_bytes_to_string, calculate_lumiScale
from htoaa_CommonTools import cut_ObjectMultiplicity, cut_ObjectPt, cut_ObjectEta, cut_ObjectPt_1
'''


# use GOldenJSON

 
printLevel = 2
nEventToReadInBatch = 1 # 0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze = 20 #-1 # 1000 # 100000 # -1
#pd.set_option('display.max_columns', None)

#print("".format())

sWeighted = "Wtd: "




    
# -----------------------------------------------------------------------------------
def get_GenPartDaughters(awkArray, index_GenPart):
    if printLevel >= 9:
        print(f"\n get_GenPartDaughters:: awkArray: {awkArray},   index_GenPart: {index_GenPart}"); sys.stdout.flush()
    
    return False

# -----------------------------------------------------------------------------------
def printVariable(sName, var):
    printInDetail=True
    #if nEventsToAnalyze == -1: printInDetail = False
    if str(type(var)) in ['numpy.ndarray', "<class 'numpy.ndarray'>"]: printInDetail = False # as gave error
    #print(f"printInDetail: {printInDetail} {sName} ({type(var)}) ({len(var)}): {var}")
    if not printInDetail:
        print(f"{sName} ({type(var)}) ({len(var)}): {var}")
    else:
        print(f"{sName} ({type(var)}) ({len(var)}): {var.to_list()}")

# -----------------------------------------------------------------------------------
class ObjectSelection:
    def __init__(self, era):
        self.era = era
        
        self.tagger_btagDeepB = 'DeepCSV'
        self.wp_btagDeepB = 'M'

        self.FatJetPtThsh  = 400 #170
        self.FatJetEtaThsh = 2.4

        self.FatJetMSoftDropThshLow  = 90
        self.FatJetMSoftDropThshHigh = 200

        self.FatJetParticleNetMD_Xbb_Thsh = 0.8
        
        self.JetPtThshForHT = 30.0
        self.JetEtaThshForHT = 2.4
        
        self.nFatJetMin = 1
        self.GenHTThsh  = 100.0
        self.LHEHTThsh  = 100.0
        

    def selectFatJets(self, events):
        '''
        Select FatJets
        '''
        '''
        jets.cut(jets['FatJet_pt'] > 170)
        jets.cut(jets['FatJet_eta'].abs() < 2.4)
        jets.cut(jets['FatJet_btagDDBvL'] > 0.8)
        jets.cut(jets['FatJet_btagDeepB'] > 0.4184)
        jets.cut(jets['FatJet_msoftdrop'] > 90)
        jets.cut(jets['FatJet_msoftdrop'] < 200)#<= 200)
        jets.cut(jets['FatJet_mass'] > 90)
        #jets.cut(jets['FatJet_mass'] <= 200)
        other.cut(other['PV_npvsGood'] >= 1)
        '''

        maskSelFatJets = (
            (events.FatJet.pt > self.FatJetPtThsh) &
            (abs(events.FatJet.eta) < self.FatJetEtaThsh) &
            (events.FatJet.btagDeepB > bTagWPs[self.era][self.tagger_btagDeepB][self.wp_btagDeepB])
        )
        if printLevel >= 15:
            #print(f"era: {self.era}, bTagWPs[self.era]: {bTagWPs[self.era]}")
            print(f"selectFatJets()::maskSelFatJets {len(maskSelFatJets)}: {maskSelFatJets.to_list()}")
        return events.FatJet[maskSelFatJets]


    def selectGenHiggs(self, events):
        maskGenHiggs = (
            (events.GenPart.pdgId  == 25) & # pdgId:: 25: H0
            (events.GenPart.status == 62)   # statu 62: outgoing subprocess particle with primordial kT included https://pythia.org/latest-manual/ParticleProperties.html
        )
        if printLevel >= 15:
            print(f"\n maskGenHiggs:  {maskGenHiggs.to_list()} ")
            print(f"\n events.GenPart[maskGenHiggs]:  {events.GenPart[maskGenHiggs].to_list()} ")
            print(f"\n events.GenPart[maskGenHiggs].mass:  {events.GenPart[maskGenHiggs].mass.to_list()} ")
        return events.GenPart[maskGenHiggs]

    def selectGenABoson(self, events):
        maskGenA = (
            (events.GenPart.pdgId == 36)
        )
        if printLevel >= 15:
            print(f"\n maskGenA:  {maskGenA.to_list()} ")
            print(f"\n events.GenPart[maskGenA]:  {events.GenPart[maskGenA].to_list()} ")
            print(f"\n events.GenPart[maskGenA].mass:  {events.GenPart[maskGenA].mass.to_list()} ")
        return events.GenPart[maskGenA]


    def GenHT(self, events):
        maskForGenHT = (
            (events.GenJet.pt > self.JetPtThshForHT) &
            (abs(events.GenJet.eta) < self.JetEtaThshForHT)
        )
        selGenJetPt = events.GenJet[maskForGenHT].pt
        GenHT = ak.sum(selGenJetPt, axis=-1)
        
        if printLevel >= 15:
            print(f"\nevents.GenJet.fields: {events.GenJet.fields}")
            print(f"\nevents.GenJet: {events.GenJet.to_list()}")
            print(f"\nevents.GenJet.pt ({len(events.GenJet.pt)}): {events.GenJet.pt.to_list()} ")
            print(f"\nmaskForGenHT: {maskForGenHT.to_list()} ")
            print(f"\nselGenJetPt: {selGenJetPt.to_list()} ")
            print(f"\nGenHT ({len(GenHT)}): {GenHT.to_list()} ")
            
        return GenHT




    
    
class HToAATo4bProcessor(processor.ProcessorABC):
    def __init__(self, datasetInfo={}):

        ak.behavior.update(nanoaod.behavior)

        self.datasetInfo = datasetInfo
        #self.isMC = isMC
        self.objectSelector = ObjectSelection(era=self.datasetInfo["era"])
        
        # set self.pdgId_BHadrons for 'QCD_bGenFilter' sample requirement
        self.pdgId_BHadrons = []
        bHadrons_ = Particle.findall(lambda p: p.pdgid.has_bottom) # Find all bottom hadrons
        print(f"List of B-hadrons for QCD B-GEN-filter ({len(bHadrons_)}):")
        print("%s %-20s %15s %15s" %(" "*4, "pdg name", "pdgId", "Mass in MeV"))
        for bHadron in bHadrons_:
            self.pdgId_BHadrons.append(bHadron.pdgid.abspid)
            print("%s %-20s %15d %15s" % (" "*4, str(bHadron), bHadron.pdgid.abspid, str(bHadron.mass)))
        print(f"self.pdgId_BHadrons ({len(self.pdgId_BHadrons)}): {self.pdgId_BHadrons}")
        self.pdgId_BHadrons = list(set(self.pdgId_BHadrons))
        print(f" after duplicate removal --> \nself.pdgId_BHadrons ({len(self.pdgId_BHadrons)}): {self.pdgId_BHadrons}")
        
        
        #dataset_axis = hist.axis.StrCategory(name="dataset", label="", categories=[], growth=True)
        #muon_axis = hist.axis.Regular(name="massT", label="Transverse Mass [GeV]", bins=50, start=15, stop=250)
        dataset_axis    = hist.Cat("dataset", "Dataset")
        systematic_axis = hist.Cat("systematic", "Systematic Uncertatinty")

        cutFlow_axis  = hist.Bin("CutFlow",   r"Cuts",            21, -0.5, 20.5)
        nObject_axis  = hist.Bin("nObject",   r"No. of object",   21, -0.5, 20.5)
        pt_axis       = hist.Bin("Pt",        r"$p_{T}$ [GeV]",   200, 0, 1000)
        ptLow_axis    = hist.Bin("PtLow",     r"$p_{T}$ [GeV]",   200, 0, 50)
        eta_axis      = hist.Bin("Eta",       r"$#eta$",          100, -5, 5)
        phi_axis      = hist.Bin("Phi",       r"$\phi$",          100, -3.14, 3.13)
        #mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       200, 0, 600)
        #mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       400, 0, 200)
        mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       300, 0, 300)
        mass_axis1    = hist.Bin("Mass1",     r"$m$ [GeV]",       300, 0, 300)
        mlScore_axis  = hist.Bin("MLScore",   r"ML score",        100, -1.1, 1.1)
        jetN2_axis    = hist.Bin("N2",        r"N2b1",            100, 0, 3)
        jetN3_axis    = hist.Bin("N3",        r"N3b1",            100, 0, 5)
        jetTau_axis   = hist.Bin("TauN",      r"TauN",            100, 0, 1)
        deltaR_axis   = hist.Bin("deltaR",    r"$delta$ r ",      500, 0, 5)
        HT_axis       = hist.Bin("HT",        r"HT",             2500, 0, 2500)
        PytPartStatus_axis = hist.Bin("PytPartStatus",        r"PytPartStatus",       421,-210.5, 210.5)
        boolean_axis  = hist.Bin("Boolean",   r"Boolean",           2,-0.5,1.5)
        pdgId_axis    = hist.Bin("PdgId",     r"PdgId",           101,-0.5,100.5)

        sXaxis      = 'xAxis'
        sXaxisLabel = 'xAxisLabel'
        sYaxis      = 'yAxis'
        sYaxisLabel = 'yAxisLabel'
        histos = OD([
            # ('histogram_name',  {sXaxis: hist.Bin() axis,  sXaxisLabel: "histogram axis label"})
            ('hCutFlow',                                  {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
            ('hCutFlowWeighted',                          {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
            
            ('nSelFatJet',                                {sXaxis: nObject_axis,    sXaxisLabel: 'No. of selected FatJets'}),
            ('hLeadingFatJetPt',                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
            ('hLeadingFatJetEta',                         {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
            ('hLeadingFatJetPhi',                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
            ('hLeadingFatJetMass',                        {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet) [GeV]"}),
            ('hLeadingFatJetMSoftDrop',                   {sXaxis: mass_axis,       sXaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
            
            ('hLeadingFatJetBtagDeepB',                   {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagDeepB"}),
            ('hLeadingFatJetBtagDDBvLV2',                 {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagDDBvLV2"}),

            ('hLeadingFatJetBtagDDCvBV2',                 {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagDDCvBV2"}),
            ('hLeadingFatJetBtagHbb',                     {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagHbb"}),
            ('hLeadingFatJetDeepTagMD_H4qvsQCD',          {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD"}),
            ('hLeadingFatJetDeepTagMD_HbbvsQCD',          {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD"}),
            ('hLeadingFatJetDeepTagMD_ZHbbvsQCD',         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD"}),
            ('hLeadingFatJetDeepTagMD_ZHccvsQCD',         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHccvsQCD"}),
            
            ('hLeadingFatJetDeepTagMD_ZbbvsQCD',          {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetdeepTagMD_ZbbvsQCD"}),
            ('hLeadingFatJetDeepTagMD_ZvsQCD',            {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZvsQCD"}),
            ('hLeadingFatJetDeepTagMD_bbvsLight',         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight"}),
            ('hLeadingFatJetDeepTagMD_ccvsLight',         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ccvsLight"}),
            ('hLeadingFatJetDeepTag_H',                   {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTag_H"}),
            ('hLeadingFatJetDeepTag_QCD',                 {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTag_QCD"}),
            ('hLeadingFatJetDeepTag_QCDothers',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTag_QCDothers"}),
            
            ('hLeadingFatJetN2b1',                        {sXaxis: jetN2_axis,      sXaxisLabel: r"LeadingFatJetn2b1"}),
            ('hLeadingFatJetN3b1',                        {sXaxis: jetN3_axis,      sXaxisLabel: r"LeadingFatJetn3b1"}),
            ('hLeadingFatJetTau1',                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau1"}),
            ('hLeadingFatJetTau2',                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau2"}),
            ('hLeadingFatJetTau3',                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau3"}),
            ('hLeadingFatJetTau4',                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau4"}),

            ('hLeadingFatJetTau4by3',                     {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau4by3"}),
            ('hLeadingFatJetTau3by2',                     {sXaxis: jetTau_axis,     sXaxisLabel: r"hLeadingFatJetTau3by2"}),
            ('hLeadingFatJetTau2by1',                     {sXaxis: jetTau_axis,     sXaxisLabel: r"hLeadingFatJetTau2by1"}),
            
            ('hLeadingFatJetNBHadrons',                   {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNBHadrons"}),
            ('hLeadingFatJetNCHadrons',                   {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNCHadrons"}),            
            ('hLeadingFatJetNConstituents',               {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNConstituents"}),
            
            ('hLeadingFatJetParticleNetMD_QCD',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_QCD"}),
            ('hLeadingFatJetParticleNetMD_Xbb',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xbb"}),
            ('hLeadingFatJetParticleNetMD_Xcc',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xcc"}),

            ('hLeadingFatJetParticleNetMD_Xqq',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xqq"}),
            ('hLeadingFatJetParticleNet_H4qvsQCD',        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_H4qvsQCD"}),
            ('hLeadingFatJetParticleNet_HbbvsQCD',        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_HbbvsQCD"}),
            ('hLeadingFatJetParticleNet_HccvsQCD',        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_HccvsQCD"}),
            ('hLeadingFatJetParticleNet_QCD',             {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_QCD"}),
            
            ('hLeadingFatJetParticleNet_mass',            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_mass"}),

            ('hGenHiggsPt_all',                           {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),
            ('hGenHiggsPt_sel',                           {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),
            ('hGenHiggsPt_sel_wGenCuts',                  {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),

            ('hGenHiggsMass_all',                         {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]"}),
            ('hMass_GenA_all',                            {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A) [GeV]"}),
            ('hMass_GenAApair_all',                       {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN HToAA) [GeV]"}),
            ('hMass_GenAToBBbarpair_all',                 {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN AToBB) [GeV]"}),
            ('hMass_Gen4BFromHToAA_all',                  {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN HTOAATo4B) [GeV]"}),
            ('hMass_GenAToBBbarpair_all_1',               {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN AToBB) [GeV]"}),
            ('hMass_Gen4BFromHToAA_all_1',                {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN HTOAATo4B) [GeV]"}),
            ('hDeltaR_GenH_GenB_max',                     {sXaxis: deltaR_axis,     sXaxisLabel: r"$Delta$r (GEN H, GEN B)_{max}"}),

            # QCD sample sticking
            ('hGenLHE_HT_all',                            {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
            ('hGenBquark_Status_all',                     {sXaxis: PytPartStatus_axis, sXaxisLabel: r"GEN Bquark Pythia status"}),
            ('hGenBquark_first_Status_all',               {sXaxis: PytPartStatus_axis, sXaxisLabel: r"GEN first Bquark Pythia status"}),
            ('hGenBquark_first_PdgId_all',               {sXaxis: pdgId_axis,       sXaxisLabel: r"GEN first Bquark pdgId"}),

            #('hGenBquark_first_isPrompt_all',             {sXaxis: boolean_axis,    sXaxisLabel: r"GEN first Bquark isPrompt"}),

            
            
            ('hGenBquark_leadingPt_all',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, leading pT) [GeV]"}),
            ('hGenBquark_subleadingPt_all',               {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, subleading pT) [GeV]"}),
            ('hGenBquark_thirdLeadingPt_all',             {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, third leading pT) [GeV]"}),
            ('hGenBquark_forthLeadingPt_all',             {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, forth leading pT) [GeV]"}),
            

            # 2-D distribution
            ('hMass_GenA1_vs_GenA2_all',                       {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A1) [GeV]",
                                                                sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A2) [GeV]"}),
            ('hMass_GenA1ToBBbar_vs_GenA2ToBBbar_all',         {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A1ToBBbar) [GeV]",
                                                                sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A2ToBBbar) [GeV]"}),
            ('hMass_GenAHeavy_vs_GenALight_all',               {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A heavy) [GeV]",
                                                                sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A light) [GeV]"}),
            ('hMass_GenH_vs_GenAHeavy_all',                    {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]",
                                                                sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A heavy) [GeV]"}),
            ('hMass_GenH_vs_GenALight_all',                    {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]",
                                                                sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A light) [GeV]"}),
            ('hMassGenH_vs_maxDRGenHGenB_all',                 {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]",
                                                                sYaxis: deltaR_axis,      sYaxisLabel: r"$Delta$r (GEN H, GEN B)_{max}"}),
            ('hMassGenAHeavy_vs_maxDRGenHGenB_all',            {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A heavy) [GeV]",
                                                                sYaxis: deltaR_axis,      sYaxisLabel: r"$Delta$r (GEN H, GEN B)_{max}"}),
            ('hMassGenALight_vs_maxDRGenHGenB_all',            {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A light) [GeV]",
                                                                sYaxis: deltaR_axis,      sYaxisLabel: r"$Delta$r (GEN H, GEN B)_{max}"}),

            
        ])

        for statusFlag_ in GENPART_STATUSFLAGS:
            histos['hGenBquark_first_%s_all' % (statusFlag_)] = {sXaxis: boolean_axis,    sXaxisLabel: r"GEN first Bquark %s"  % (statusFlag_)}
        

        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int)
        })

        for histName, histAttributes in histos.items():
            #hXaxis = histAttributes[sXaxis].copy()
            hXaxis = deepcopy(histAttributes[sXaxis])
            hXaxis.label = histAttributes[sXaxisLabel]

            if sYaxis not in histAttributes.keys():
                # TH1
                self._accumulator.add({
                    histName: hist.Hist(
                        "Counts",
                        dataset_axis,
                        hXaxis, #nObject_axis,
                        systematic_axis,
                    )
                })
            else:
                # TH2
                hYaxis = deepcopy(histAttributes[sYaxis])
                hYaxis.label = histAttributes[sYaxisLabel]
                
                self._accumulator.add({
                    histName: hist.Hist(
                        "Counts",
                        dataset_axis,
                        hXaxis, #nObject_axis,
                        hYaxis,
                        systematic_axis,
                    )
                })
                


                

        '''
        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int),
            #'hnFatJet_level0': hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of FatJet").Weight(),
            #'hnFatJet_level1': hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of FatJet").Weight(),
            'nSelFatJet': hist.Hist(
                "Counts",
                dataset_axis,
                nObject_axis,
                systematic_axis,
            ),
            'hLeadingFatJetPt': hist.Hist(
                "Counts",
                dataset_axis,
                pt_axis,
                systematic_axis,
            ),
            
        })

        self._accumulator.add({
             'hLeadingFatJetDeepTagMD_HbbvsQCD_1': hist.Hist(
                "Counts",
                dataset_axis,
                mlScore_axis,
                systematic_axis,
            ),
        })
        '''


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata["dataset"] # dataset label
        self.datasetInfo[dataset]['isSignal'] = False
        self.datasetInfo[dataset]['isQCD'] = False

        if printLevel >= 20:
            print(f"nEvents: {len(events)}")
        if printLevel >= 20:
            print(f"\n events.fields: {events.fields}")
            print(f"\n events.GenPart.fields: {events.GenPart.fields}")
            #print(f"\n events.HLT.fields: {events.HLT.fields}")
            #printVariable('\n events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4', events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4)
            #print(f"\n events.L1.fields: {events.L1.fields}")
            #printVariable('\n events.L1.SingleJet180', events.L1.SingleJet180)
            #print(f"\n events.FatJet.fields: {events.FatJet.fields}")
            #print(f"\n events.LHE.fields: {events.LHE.fields}")
            #print(f"\n events.LHE.HT: {events.LHE.HT.to_list()}")
            
        if nEventToReadInBatch == 1:
            print(f"\n (run:ls:event): {ak.zip([events.run, events.luminosityBlock, events.event])}")

        self.datasetInfo[dataset]['isSignal'] = False
        self.datasetInfo[dataset]['isQCD'] = False
        self.datasetInfo[dataset]['isQCDIncl'] = False
        self.datasetInfo[dataset]['isQCD_bEnrich'] = False
        self.datasetInfo[dataset]['isQCD_bGen'] = False
        
        if self.datasetInfo[dataset]['isMC']:
            self.datasetInfo[dataset]['isSignal']      = True if "HToAATo4B"   in dataset else False
            #self.datasetInfo[dataset]['isQCD']         = True if "QCD"         in dataset else False
            self.datasetInfo[dataset]['isQCDIncl']     = True if kQCDIncl      in dataset else False
            self.datasetInfo[dataset]['isQCD_bEnrich'] = True if kQCD_bEnrich  in dataset else False
            self.datasetInfo[dataset]['isQCD_bGen']    = True if kQCD_bGen     in dataset else False
            self.datasetInfo[dataset]['isQCD'] = self.datasetInfo[dataset]['isQCDIncl'] or self.datasetInfo[dataset]['isQCD_bEnrich'] or self.datasetInfo[dataset]['isQCD_bGen']
               
            
            output = self.accumulator.identity()
            systematics_shift = [None]
            for _syst in systematics_shift:
                output += self.process_shift(events, _syst)
        else:
            output = self.process_shift(events, None)

        return output


    
    def process_shift(self, events, shift_syst=None):
        
        output = self.accumulator.identity()
        dataset = events.metadata["dataset"] # dataset label
        #print(f"dataset: {dataset}")

        

        ones_list  = np.ones(len(events))
        trues_list = np.ones(len(events), dtype=bool)


          
        ##################
        # OBJECT SELECTION
        ##################


        # Gen-level selection ---------------------------------------------------------------------
        genHiggs = None
        genHT    = None
        if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isSignal']: 
            genHiggs  = self.objectSelector.selectGenHiggs(events)        
            genHT     = self.objectSelector.GenHT(events)

            # m(bbar from A) and m(4b from HToAA) ----------
            '''
            genHiggsCollection = events.GenPart[(
                (events.GenPart.pdgId  == 25) & # pdgId:: 25: H0
                (events.GenPart.status == 62)   # statu 62: outgoing subprocess particle with primordial kT included https://pythia.org/latest-manual/ParticleProperties.html
            )]
            '''

            genACollection = self.objectSelector.selectGenABoson(events)
            genA_First  = genACollection[:, 0]
            genA_Second = genACollection[:, 1]
            
            idxGenA_sortByMass = ak.argsort(genACollection.mass, axis=-1, ascending=False)
            # genACollection[idxGenA_sortByMass[0]] : Leading mass GenA
            # genACollection[idxGenA_sortByMass[1]] : Subleading mass GenA
            
            
            genBBar_pairs_all = ak.argcombinations(events.GenPart, 2, fields=['b', 'bbar'])
            genBBar_pairs = genBBar_pairs_all[(
                (abs(events.GenPart[genBBar_pairs_all['b'   ]].pdgId) == 5) &
                (abs(events.GenPart[genBBar_pairs_all['bbar']].pdgId) == 5) &
                ((events.GenPart[genBBar_pairs_all['b']].pdgId) == (-1*events.GenPart[genBBar_pairs_all['bbar']].pdgId)  ) &
                (events.GenPart[genBBar_pairs_all['b']].genPartIdxMother == events.GenPart[genBBar_pairs_all['bbar']].genPartIdxMother) &
                (events.GenPart[ events.GenPart[genBBar_pairs_all['b'   ]].genPartIdxMother ].pdgId == 36) &
                (events.GenPart[ events.GenPart[genBBar_pairs_all['bbar']].genPartIdxMother ].pdgId == 36) 
            )]

            # LorentVector of GenB quarks from HToAATo4b
            nEvents_11 = ak.num(events.GenPart[genBBar_pairs['b']][:, 0].pt, axis=0)
            mass_bQuark = 4.18
            #print(f"\n np.full(nEvents_11, mass_bQuark): {np.full(nEvents_11, mass_bQuark)}")

            # https://coffeateam.github.io/coffea/modules/coffea.nanoevents.methods.vector.html
            # bQuark from 1st A
            LVGenB_0 = ak.zip(
                {
                    "pt"  : events.GenPart[genBBar_pairs['b']][:, 0].pt,
                    "eta" : events.GenPart[genBBar_pairs['b']][:, 0].eta,
                    "phi" : events.GenPart[genBBar_pairs['b']][:, 0].phi,
                    "mass": np.full(nEvents_11, mass_bQuark),
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )

            # bbarQuark from 1st A
            LVGenBbar_0 = ak.zip(
                {
                    "pt"  : events.GenPart[genBBar_pairs['bbar']][:, 0].pt,
                    "eta" : events.GenPart[genBBar_pairs['bbar']][:, 0].eta,
                    "phi" : events.GenPart[genBBar_pairs['bbar']][:, 0].phi,
                    "mass": np.full(nEvents_11, mass_bQuark),
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )

            # bQuark from 2nd A
            LVGenB_1 = ak.zip(
                {
                    "pt"  : events.GenPart[genBBar_pairs['b']][:, 1].pt,
                    "eta" : events.GenPart[genBBar_pairs['b']][:, 1].eta,
                    "phi" : events.GenPart[genBBar_pairs['b']][:, 1].phi,
                    "mass": np.full(nEvents_11, mass_bQuark),
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )

            # bbarQuark from 2nd A
            LVGenBbar_1 = ak.zip(
                {
                    "pt"  : events.GenPart[genBBar_pairs['bbar']][:, 1].pt,
                    "eta" : events.GenPart[genBBar_pairs['bbar']][:, 1].eta,
                    "phi" : events.GenPart[genBBar_pairs['bbar']][:, 1].phi,
                    "mass": np.full(nEvents_11, mass_bQuark),
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )            

            dr_GenH_GenB = ak.concatenate([genHiggs.delta_r(LVGenB_0), genHiggs.delta_r(LVGenBbar_0), genHiggs.delta_r(LVGenB_1), genHiggs.delta_r(LVGenBbar_1)], axis=-1)
            max_dr_GenH_GenB = ak.max(dr_GenH_GenB, axis=-1)            

                
        # QCD MC ----------------------------------------------
        mask_genBHadrons_status2_eventwise   = None
        mask_genBQuarks_hardSctred_eventwise = None
        mask_QCD_stitch_eventwise            = None        
        if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isQCD'] :
            genBQuarks_QCD = events.GenPart[(
                (abs(events.GenPart.pdgId) == 5)
            )]
            mask_genBQuarks_wMotherIdx = (
                (abs(events.GenPart.pdgId) == 5) &
                (events.GenPart.genPartIdxMother >= 0)  
            )
            genBQuarks_pT = ak.sort(genBQuarks_QCD.pt, axis=-1, ascending=False)
            genBQuarks_first = ak.firsts(genBQuarks_QCD)
            #genBQuarks_first = ak.firsts(genBQuarks_QCD, axis=-1)
            mask_genBQuarks = (ak.count(genBQuarks_QCD.pdgId, axis=1) >= 1)

            
            if printLevel >= 4:
                #printVariable('\n ', ); sys.stdout.flush()
                printVariable('\n mask_genBQuarks', mask_genBQuarks); sys.stdout.flush()
                printVariable('\n genBQuarks_QCD[mask_genBQuarks]', genBQuarks_QCD[mask_genBQuarks]); sys.stdout.flush()
                printVariable('\n genBQuarks_QCD.pdgId', genBQuarks_QCD.pdgId); sys.stdout.flush()
                printVariable('\n genBQuarks_first.pdgId', genBQuarks_first.pdgId); sys.stdout.flush()
                printVariable('\n abs(genBQuarks_first.pdgId)', abs(genBQuarks_first.pdgId)); sys.stdout.flush()
                #print(f"\n ak.drop_none(abs(genBQuarks_first.pdgId)): {ak.drop_none(abs(genBQuarks_first.pdgId))}"); sys.stdout.flush()
                print(f"\n ak.to_numpy(abs(genBQuarks_first.pdgId)) {type(ak.to_numpy(abs(genBQuarks_first.pdgId)))}: {ak.to_numpy(abs(genBQuarks_first.pdgId))}"); sys.stdout.flush()
                print(f"\n ak.to_numpy(abs(genBQuarks_first[mask_genBQuarks].pdgId)) {type(ak.to_numpy(abs(genBQuarks_first[mask_genBQuarks].pdgId)))}: {ak.to_numpy(abs(genBQuarks_first[mask_genBQuarks].pdgId))}"); sys.stdout.flush()
                #print(f"\n ak.flatten(abs(genBQuarks_first.pdgId)): {ak.flatten(abs(genBQuarks_first.pdgId))}"); sys.stdout.flush()
                
                
            if printLevel >= 4:
                #print(f"genBQuarks_QCD")
                #printVariable('\n ', ); sys.stdout.flush()
                #printVariable('\n mask: genBQuarks_QCD', (abs(events.GenPart.pdgId) == 5)); sys.stdout.flush()
                #printVariable('\n events.GenPart', events.GenPart); sys.stdout.flush()
                #printVariable('\n genBQuarks_QCD', genBQuarks_QCD); sys.stdout.flush()
                printVariable('\n genBQuarks_QCD.status', genBQuarks_QCD.status); sys.stdout.flush()
                printVariable('\n genBQuarks_QCD.statusFlags', genBQuarks_QCD.statusFlags); sys.stdout.flush()
                printVariable('\n genBQuarks_QCD.genPartIdxMother', genBQuarks_QCD.genPartIdxMother); sys.stdout.flush()
                #printVariable('\n genBQuarks_QCD.genPartIdxMotherG', genBQuarks_QCD.genPartIdxMotherG); sys.stdout.flush()
                #printVariable('\n genBQuarks_QCD.distinctParentIdxG', genBQuarks_QCD.distinctParentIdxG); sys.stdout.flush()
                #printVariable('\n genBQuarks_QCD.childrenIdxG', genBQuarks_QCD.childrenIdxG); sys.stdout.flush()
                #printVariable('\n genBQuarks_QCD.distinctChildrenIdxG', genBQuarks_QCD.distinctChildrenIdxG); sys.stdout.flush()
                #printVariable('\n genBQuarks_QCD', genBQuarks_QCD)
                
                #printVariable('\n ak.first(genBQuarks_QCD)', ak.firsts(genBQuarks_QCD)); sys.stdout.flush()

                # genBQuarks_QCD[:, 0] status, IdxMother

                for flag_ in GENPART_STATUSFLAGS:
                    printVariable('\n genBQuarks_first.hasFlags(%s))' % (flag_), genBQuarks_first.hasFlags(flag_)); sys.stdout.flush()
                printVariable('\n genBQuarks_first.parent', genBQuarks_first.parent); sys.stdout.flush()
                printVariable('\n genBQuarks_first.parent.pdgId', genBQuarks_first.parent.pdgId); sys.stdout.flush()
                printVariable('\n genBQuarks_first.parent.status', genBQuarks_first.parent.status); sys.stdout.flush()
                
                printVariable('\n genBQuarks_first.distinctParent', genBQuarks_first.distinctParent); sys.stdout.flush()
                printVariable('\n genBQuarks_first.distinctParent.pdgId', genBQuarks_first.distinctParent.pdgId); sys.stdout.flush()
                printVariable('\n genBQuarks_first.distinctParent.status', genBQuarks_first.distinctParent.status); sys.stdout.flush()

            
            if printLevel >= 4:
                #printVariable('\n ', ); sys.stdout.flush()
                printVariable('\n *** genBQuarks_first', genBQuarks_first); sys.stdout.flush()
                #printVariable('\n mask_genBQuarks_wMotherIdx', mask_genBQuarks_wMotherIdx); sys.stdout.flush()
                #printVariable('\n events.GenPart[mask_genBQuarks_wMotherIdx]', events.GenPart[mask_genBQuarks_wMotherIdx]); sys.stdout.flush()
                #printVariable('\n events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother', events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother); sys.stdout.flush()
                printVariable('\n ak.firsts( events.GenPart[ events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ] )', ak.firsts( events.GenPart[ events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ] ) ); sys.stdout.flush()
                printVariable('\n ak.firsts( events.GenPart[ events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ].pdgId )', ak.firsts( events.GenPart[ events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ].pdgId )); sys.stdout.flush()


                '''
                printVariable('\n ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother)', ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ) ); sys.stdout.flush()
                #printVariable('\n ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother, highlevel=False)', ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother, highlevel=False ) ); sys.stdout.flush()
                print(f"\n ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother, highlevel=False ): {ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother, highlevel=False )}")
                #print(f"\n ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother, behavior=dict ): {ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother, behavior=dict )}")
                print(f"\n ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother, behavior=ak.behavior ): {ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother, behavior=ak.behavior )}")
                
                #printVariable('\n events.GenPart[ ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ) ]', events.GenPart[ ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ) ]); sys.stdout.flush()
                printVariable('\n events.GenPart[ ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ) ].pdgId', events.GenPart[ ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ) ].pdgId); sys.stdout.flush()
                #printVariable('\n ak.count( events.GenPart[ ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ) ] )', ak.count( events.GenPart[ ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ) ] )); sys.stdout.flush()
                print(f"\n ak.count( events.GenPart[ ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ) ].pdgId, axis=-1 ): {  ak.count( events.GenPart[ ak.firsts( events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ) ].pdgId, axis=-1 ).to_list() }")
                '''

                

                
                print(f"\n ak.count(events.GenPart[ events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ].pdgId, axis=-1): {ak.count(events.GenPart[ events.GenPart[mask_genBQuarks_wMotherIdx].genPartIdxMother ].pdgId, axis=-1).to_list()}")
                print(f"\n ak.count(events.GenPart.pdgId, axis=-1): {ak.count(events.GenPart.pdgId, axis=-1).to_list()}")
                
                '''
                printVariable('\n events.GenPart[(genBQuarks_first.genPartIdxMother)]', events.GenPart[(genBQuarks_first.genPartIdxMother)]); sys.stdout.flush()
                printVariable('\n events.GenPart[(genBQuarks_first.genPartIdxMother)].pdgId', events.GenPart[(genBQuarks_first.genPartIdxMother)].pdgId); sys.stdout.flush()
                printVariable('\n events.GenPart[(genBQuarks_first.genPartIdxMother)].status', events.GenPart[(genBQuarks_first.genPartIdxMother)].status); sys.stdout.flush()
                '''


            '''    
            def bitwiseTransformer(layout, **kwarg):
                print(f"layout {type(layout)} : {layout}")

            #genPart_copy = ak.with_field(events.GenPart)
            #genPart_copy = events.GenPart.statusFlags
            genPart_copy = copy(events.GenPart)
            #genPart_copy['statusFlags_inBits'] = genPart_copy.statusFlags.apply(lambda x: bin(x), axis=1)
            #genPart_copy1 = ak.transform(genPart_copy.statusFlags)
            #ak.transform(bitwiseTransformer, genPart_copy.statusFlags)

            if printLevel >= 2:
                printVariable('\n genPart_copy', genPart_copy); sys.stdout.flush()
                printVariable('\n genPart_copy[statusFlags]', genPart_copy['statusFlags']); sys.stdout.flush()
                printVariable('\n events.GenPart.statusFlags', events.GenPart.statusFlags); sys.stdout.flush()
                #print(f"bin(genPart_copy.statusFlags): {bin(genPart_copy.statusFlags)}")

                #print(f"genBQuarks_first.statusFlags: {genBQuarks_first.statusFlags}")
                #print(f"np.binary_repr(genBQuarks_first.statusFlags): {np.binary_repr(genBQuarks_first.statusFlags)}")

                #printVariable('\n ', ); sys.stdout.flush()
                #printVariable('\n genBQuarks_first.statusFlags', genBQuarks_first.statusFlags); sys.stdout.flush()
                #printVariable('\n ak.drop_none(genBQuarks_first.statusFlags)', ak.drop_none(genBQuarks_first.statusFlags)); sys.stdout.flush()
                #printVariable('\n ak.fill_none(genBQuarks_first.statusFlags, value=0)', ak.fill_none(genBQuarks_first.statusFlags, value=0)); sys.stdout.flush()
                #printVariable('\n np.unpackbits(genBQuarks_first.statusFlags)', np.unpackbits(genBQuarks_first.statusFlags)); sys.stdout.flush()
                #print(f"\n np.unpackbits(ak.fill_none(genBQuarks_first.statusFlags, value=0)): {np.unpackbits(ak.fill_none(genBQuarks_first.statusFlags, value=0))}")

            '''

            #print(f"events.GenPart.hasFlags(['isHardProcess']): {events.GenPart.hasFlags(['isHardProcess'])}")
            
            if printLevel >= 4:
                mask_genBQuarks_first_isHardProcess = (
                    #(ak.num(genBQuarks_QCD.pdgId, axis=1) >=1) &
                    (genBQuarks_first.hasFlags('isHardProcess'))                    
                ) 
                mask_genBQuarks_first_isHardProcess = ak.fill_none(mask_genBQuarks_first_isHardProcess, value=False)
                #printVariable("\n mask_genBQuarks_first_isHardProcess", mask_genBQuarks_first_isHardProcess)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess]", genBQuarks_first[mask_genBQuarks_first_isHardProcess])
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].pdgId", genBQuarks_first[mask_genBQuarks_first_isHardProcess].pdgId)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].status", genBQuarks_first[mask_genBQuarks_first_isHardProcess].status)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].genPartIdxMother", genBQuarks_first[mask_genBQuarks_first_isHardProcess].genPartIdxMother)
                
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].parent", genBQuarks_first[mask_genBQuarks_first_isHardProcess].parent)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctParent", genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctParent)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].children", genBQuarks_first[mask_genBQuarks_first_isHardProcess].children)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctChildren", genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctChildren)
                #printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctChildrenDeep", genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctChildrenDeep)
                #printVariable("\n ", )
                #printVariable("\n ", )
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].parent.pdgId", genBQuarks_first[mask_genBQuarks_first_isHardProcess].parent.pdgId)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].parent.status", genBQuarks_first[mask_genBQuarks_first_isHardProcess].parent.status)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctParent.pdgId", genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctParent.pdgId)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctParent.status", genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctParent.status)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].children.pdgId", genBQuarks_first[mask_genBQuarks_first_isHardProcess].children.pdgId)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].children.status", genBQuarks_first[mask_genBQuarks_first_isHardProcess].children.status)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctChildren.pdgId", genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctChildren.pdgId)
                printVariable("\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctChildren.status", genBQuarks_first[mask_genBQuarks_first_isHardProcess].distinctChildren.status)
            
                for flag_ in GENPART_STATUSFLAGS:
                    printVariable('\n genBQuarks_first[mask_genBQuarks_first_isHardProcess].hasFlags(%s))' % (flag_), genBQuarks_first[mask_genBQuarks_first_isHardProcess].hasFlags(flag_)); sys.stdout.flush()

                printVariable("\n\n\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess]", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess])
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].pdgId", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].pdgId)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].status", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].status)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].genPartIdxMother", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].genPartIdxMother)

                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].parent", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].parent)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctParent", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctParent)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].children", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].children)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctChildren", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctChildren)
                #printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctChildrenDeep", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctChildrenDeep)
                
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].parent.pdgId", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].parent.pdgId)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].parent.status", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].parent.status)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctParent.pdgId", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctParent.pdgId)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctParent.status", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctParent.status)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].children.pdgId", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].children.pdgId)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].children.status", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].children.status)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctChildren.pdgId", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctChildren.pdgId)
                printVariable("\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctChildren.status", genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].distinctChildren.status)

                for flag_ in GENPART_STATUSFLAGS:
                    printVariable('\n genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].hasFlags(%s))' % (flag_), genBQuarks_QCD[mask_genBQuarks_first_isHardProcess].hasFlags(flag_)); sys.stdout.flush()
                


            if printLevel >= 2:
                printVariable('\n events.GenPart[(events.GenPart.status == 2)].pdgId', events.GenPart[(events.GenPart.status == 2)].pdgId); sys.stdout.flush()


            mask_genBHadrons_status2 = None
            #PDGID_BHadrons = [511, 521, 513, 523]
            for ipdgId_tmp in range(len(self.pdgId_BHadrons)):
                pdgId_tmp = self.pdgId_BHadrons[ipdgId_tmp]
                mask_genBHadrons_status2_i = (
                    (events.GenPart.status == 2) &
                    (abs(events.GenPart.pdgId) == pdgId_tmp)
                )

                if ipdgId_tmp == 0:
                    mask_genBHadrons_status2 = mask_genBHadrons_status2_i
                else:
                    mask_genBHadrons_status2 = mask_genBHadrons_status2 | mask_genBHadrons_status2_i

                if printLevel >= 100:
                    printVariable('\n\n %d \t mask_genBHadrons_status2_i' % pdgId_tmp, mask_genBHadrons_status2_i); sys.stdout.flush()
                    printVariable('\n     \t mask_genBHadrons_status2  ', mask_genBHadrons_status2); sys.stdout.flush()
                    
            genBHadrons_status2 = events.GenPart[mask_genBHadrons_status2]
            mask_genBHadrons_status2_eventwise = ak.any(mask_genBHadrons_status2, axis=1)
                    
            if printLevel >= 2:
                printVariable('\n genBHadrons_status2', genBHadrons_status2); sys.stdout.flush()
                printVariable(' genBHadrons_status2.pdgId', genBHadrons_status2.pdgId); sys.stdout.flush()
                printVariable(' ak.any(mask_genBHadrons_status2, axis=1)', ak.any(mask_genBHadrons_status2, axis=1)); sys.stdout.flush()

            mask_genBQuarks_hardSctred = (
                (abs(events.GenPart.pdgId) == 5) &
                (abs(events.GenPart.status) == 23)
            )
            genBQuarks_hardSctred = events.GenPart[mask_genBQuarks_hardSctred]
            mask_genBQuarks_hardSctred_eventwise = ak.any(mask_genBQuarks_hardSctred, axis=1)
            
            if printLevel >= 2:
                printVariable('\n genBQuarks_hardSctred', genBQuarks_hardSctred); sys.stdout.flush()
                printVariable(' genBQuarks_hardSctred.hasFlags("isHardProcess")', genBQuarks_hardSctred.hasFlags("isHardProcess")); sys.stdout.flush()
                printVariable(' genBQuarks_hardSctred.hasFlags("fromHardProcess")', genBQuarks_hardSctred.hasFlags("fromHardProcess")); sys.stdout.flush()
                printVariable(' ak.any(mask_genBQuarks_hardSctred, axis=1)', ak.any(mask_genBQuarks_hardSctred, axis=1)); sys.stdout.flush()


            if self.datasetInfo[dataset]['isQCD_bEnrich']:
                mask_QCD_stitch_eventwise = trues_list
            elif self.datasetInfo[dataset]['isQCD_bGen']:
                mask_QCD_stitch_eventwise = (
                    (mask_genBQuarks_hardSctred_eventwise == False)
                )
            elif self.datasetInfo[dataset]['isQCDIncl']:
                mask_QCD_stitch_eventwise = (
                    (mask_genBQuarks_hardSctred_eventwise == False) &
                    (mask_genBHadrons_status2_eventwise == False)
                )
                
            if printLevel >= 2:
                printVariable('\n mask_QCD_stitch_eventwise', mask_QCD_stitch_eventwise); sys.stdout.flush()
            

                    
        # Reco-level -----------------------------------------------------------------------------------
        # FatJet selection
        #selFatJet = self.objectSelector.selectFatJets(events)

        '''
        mask_FatJetPt = (events.FatJet.pt > self.objectSelector.FatJetPtThsh)
        mask_FatJetEta = (abs(events.FatJet.eta) < self.objectSelector.FatJetEtaThsh)
        mask_FatJetBtagDeepB = (events.FatJet.btagDeepB > bTagWPs[self.objectSelector.era][self.objectSelector.tagger_btagDeepB][self.objectSelector.wp_btagDeepB])
        mask_FatJetMSoftDrop = (
            (events.FatJet.msoftdrop > self.objectSelector.FatJetMSoftDropThshLow) &
            (events.FatJet.msoftdrop < self.objectSelector.FatJetMSoftDropThshHigh)
        )
        
        selFatJet = events.FatJet[(
            mask_FatJetPt &
            mask_FatJetEta &
            mask_FatJetBtagDeepB #&
            #mask_FatJetMSoftDrop
        )]
        '''

        
        ##################
        # EVENT VARIABLES
        ##################
        
        #leadingFatJet = ak.firsts(selFatJet)
        #leadingFatJet_asSingletons = ak.singletons(leadingFatJet)

        #leadingFatJet = ak.singletons( ak.firsts(events.FatJet) )
        leadingFatJet = ak.firsts(events.FatJet) # for e.g. [0.056304931640625, None, 0.12890625, 0.939453125, 0.0316162109375]
        leadingFatJet_asSingletons = ak.singletons(leadingFatJet) # for e.g. [[0.056304931640625], [], [0.12890625], [0.939453125], [0.0316162109375]]
        
        if printLevel >= 13:
            #printVariable("\n ", )
            printVariable("\n ak.firsts(events.FatJet)", ak.firsts(events.FatJet))
            printVariable("\n ak.singletons( ak.firsts(events.FatJet) )", ak.singletons( ak.firsts(events.FatJet) ))

            printVariable("\n leadingFatJet.btagDeepB", leadingFatJet.btagDeepB)
            printVariable("\n ak.singletons(leadingFatJet.btagDeepB)", ak.singletons(leadingFatJet.btagDeepB))
            printVariable("\n leadingFatJet.btagDeepB > bTagWPs[self.objectSelector.era][self.objectSelector.tagger_btagDeepB][self.objectSelector.wp_btagDeepB]", leadingFatJet.btagDeepB > bTagWPs[self.objectSelector.era][self.objectSelector.tagger_btagDeepB][self.objectSelector.wp_btagDeepB])

        #if printLevel >= 2:
        #    print(f" : {}")
            
        #####################
        # EVENT SELECTION
        #####################

        HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4" 
        
        # sel_names_all = dict of {"selection name" : [list of different cuts]}; for cut-flow table
        sel_names_all = OD([
            ("SR",                    [
                "nPV",
                "leadingFatJetPt",
                "leadingFatJetEta",
                #"leadingFatJetBtagDeepB",
                "leadingFatJetMSoftDrop",
                "leadingFatJetParticleNetMD_Xbb",
                "L1_SingleJet180",
                HLT_AK8PFJet330_name
            ]),
        ])
        if self.datasetInfo[dataset]['isQCD']: sel_names_all["SR"].append("QCDStitch")
        # reconstruction level cuts for cut-flow table. Order of cuts is IMPORTANT
        cuts_reco = ["dR_LeadingFatJet_GenB_0p8"] + sel_names_all["SR"] #.copy()

        
        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        selection = PackedSelection()

        if "nPV" in sel_names_all["SR"]:
            # nPVGood >= 1
            selection.add("nPV", events.PV.npvsGood >= 1)

        

        if "leadingFatJetPt" in sel_names_all["SR"]:
            # >=1 FatJet
            #selection.add("FatJetGet", ak.num(selFatJet) >= self.objectSelector.nFatJetMin)
            selection.add(
                "leadingFatJetPt",
                leadingFatJet.pt > self.objectSelector.FatJetPtThsh
            )


        if "leadingFatJetEta" in sel_names_all["SR"]:
            selection.add(
                "leadingFatJetEta",
                abs(leadingFatJet.eta) < self.objectSelector.FatJetEtaThsh
            )


        if "leadingFatJetBtagDeepB" in sel_names_all["SR"]:
            selection.add(
                "leadingFatJetBtagDeepB",
                leadingFatJet.btagDeepB > bTagWPs[self.objectSelector.era][self.objectSelector.tagger_btagDeepB][self.objectSelector.wp_btagDeepB]
            )

 
        if"leadingFatJetMSoftDrop"  in sel_names_all["SR"]:
           selection.add(
                "leadingFatJetMSoftDrop",
                (leadingFatJet.msoftdrop > self.objectSelector.FatJetMSoftDropThshLow) &
                (leadingFatJet.msoftdrop < self.objectSelector.FatJetMSoftDropThshHigh)
            )

 
        if "leadingFatJetParticleNetMD_Xbb" in sel_names_all["SR"]:
            selection.add(
                "leadingFatJetParticleNetMD_Xbb",
                leadingFatJet.particleNetMD_Xbb > self.objectSelector.FatJetParticleNetMD_Xbb_Thsh
            )
            
            
        if "L1_SingleJet180" in sel_names_all["SR"]:
            selection.add(
                "L1_SingleJet180",
                events.L1.SingleJet180 == True
            )

 
        if HLT_AK8PFJet330_name in sel_names_all["SR"]:
            # some files of Run2018A do not have HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 branch
            #HLT_AK8PFJet330_name = None
            if "AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4" in events.HLT.fields:
                #HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4"
                selection.add(
                    HLT_AK8PFJet330_name,
                    events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 == True
                )
            elif "AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02" in events.HLT.fields:
                #HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02"
                selection.add(
                    HLT_AK8PFJet330_name,
                    events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02 == True
                )
        

        if "QCDStitch" in sel_names_all["SR"]:
            selection.add(
                "QCDStitch",
                mask_QCD_stitch_eventwise == True
            )
            
            

            
        #sel_SR          = selection.all("nPV", "FatJetGet")
        sel_SR           = selection.all(* sel_names_all["SR"])
        sel_GenHToAATo4B = None

        if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isSignal']:
            # max. dR(sel_leadingFatJet, GEN 4B from H->aa)
            dr_LeadingFatJet_GenB = ak.concatenate([leadingFatJet_asSingletons.delta_r(LVGenB_0), leadingFatJet_asSingletons.delta_r(LVGenBbar_0), leadingFatJet_asSingletons.delta_r(LVGenB_1), leadingFatJet_asSingletons.delta_r(LVGenBbar_1)], axis=-1)
            max_dr_LeadingFatJet_GenB = ak.max(dr_LeadingFatJet_GenB, axis=-1)

            if printLevel >= 13:
                printVariable("\n leadingFatJet_asSingletons.delta_r(LVGenB_0)", leadingFatJet_asSingletons.delta_r(LVGenB_0))
                printVariable("\n dr_LeadingFatJet_GenB", dr_LeadingFatJet_GenB)
                printVariable("\n max_dr_LeadingFatJet_GenB", max_dr_LeadingFatJet_GenB)
                
            # GEN level selection
            selection.add("1GenHiggs", ak.num(genHiggs) == 1)
            selection.add("2GenA", ak.num(genACollection) == 2)
            selection.add("2GenAToBBbarPairs", ak.num(genBBar_pairs) == 2)
            selection.add("dR_GenH_GenB_0p8", max_dr_GenH_GenB < 0.8)
            selection.add("dR_LeadingFatJet_GenB_0p8", max_dr_LeadingFatJet_GenB < 0.8)

            # 
            sel_names_GEN = ["1GenHiggs", "2GenA", "2GenAToBBbarPairs", "dR_GenH_GenB_0p8"]
            sel_names_all.update( OD([
                ("GenHToAATo4B_1", ["1GenHiggs", "2GenA", "2GenAToBBbarPairs"]),
                ("GenHToAATo4B", [*sel_names_GEN]),
            ]) )
            for idx, cutName in enumerate(cuts_reco):
                if idx == 0:
                    sel_names_all.update( OD([
                        ("SR_%d" % (idx+1),  [*sel_names_GEN, cutName]),
                    ]) )
                else:
                    sel_names_all.update( OD([
                        ("SR_%d" % (idx+1),  [*sel_names_all["SR_%d" % (idx)], cutName]),
                    ]) ) 
            '''
            sel_names_all.update( OD([
                ("GenHToAATo4B_1", ["1GenHiggs", "2GenA", "2GenAToBBbarPairs"]),
                ("GenHToAATo4B", [*sel_names_GEN]),

                ("SR_1",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8"]),
                ("SR_2",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt"]),
                ("SR_3",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt", "FatJetEta"]),
                ("SR_4",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt", "FatJetEta", "FatJetBtagDeepB"]),
                ("SR_5",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt", "FatJetEta", "FatJetBtagDeepB", "FatJetMSoftDrop"]),
            ]) )
            '''
            
            sel_GenHToAATo4B = selection.all(* sel_names_all["GenHToAATo4B"])

        #print(f"\nsel_SR ({len(sel_SR)}): {sel_SR}   nEventsPass: {ak.sum(sel_SR, axis=0)}")
        #print(f"\nsel_SR_wGenCuts ({len(sel_SR_wGenCuts)}): {sel_SR_wGenCuts}   nEventsPass: {ak.sum(sel_SR_wGenCuts, axis=0)}")
        
        
        # useful debugger for selection efficiency
        if shift_syst is None and printLevel >= 5:
            print(dataset)
            for n in selection.names:
                print(
                    f"- Cut {n} pass {selection.all(n).sum()} of {len(events)} events"
                )
                print(f"selection {n} ({type(selection.all(n))}): {selection.all(n)}")
                #wgt1=np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
                #print(f"wgt1 ({len(wgt1)}): {wgt1}")
                

        if printLevel >= 15:
            print(f"selFatJet.fields: {selFatJet.fields}")
            print(f"sel_SR ({len(sel_SR)}): {sel_SR}")
            print(f"selFatJet.pt[sel_SR].to_list(): {selFatJet.pt[sel_SR].to_list()} ")

        if printLevel >= 30:
            #print(f" ({type()}) ({len()}): {} \n")
            print(f"leadingFatJet ({type(leadingFatJet)}) ({len(leadingFatJet)}): {leadingFatJet} \n")
            print(f"leadingFatJet_asSingletons ({type(leadingFatJet_asSingletons)}) ({len(leadingFatJet_asSingletons)}): {leadingFatJet_asSingletons} \n")
            print(f"sel_SR ({type(sel_SR)}) ({len(sel_SR)}): {sel_SR} \n")
            print(f"sel_GenHToAATo4B ({type(sel_GenHToAATo4B)}) ({len(sel_GenHToAATo4B)}): {sel_GenHToAATo4B} \n")
            print(f"leadingFatJet ({type(leadingFatJet)}) ({len(leadingFatJet)}): {leadingFatJet} \n")
            print(f"leadingFatJet.pt[sel_SR] ({type(leadingFatJet.pt[sel_SR])}) ({len(leadingFatJet.pt[sel_SR])}): {leadingFatJet.pt[sel_SR]} \n")
            #print(f" ({type()}) ({len()}): {} \n")
                    
        
        if printLevel >= 30:
            #printVariable("", )
            printVariable("\n events.FatJet.pt", events.FatJet.pt)
            printVariable("\n mask_FatJetPt", mask_FatJetPt)
            printVariable("\n events.FatJet[mask_FatJetPt].pt", events.FatJet[mask_FatJetPt].pt)            
            printVariable("\n ak.num(events.FatJet[mask_FatJetPt])", ak.num(events.FatJet[mask_FatJetPt]))
            printVariable("\n ak.num(events.FatJet[mask_FatJetPt]) >= self.objectSelector.nFatJetMin", ak.num(events.FatJet[mask_FatJetPt]) >= self.objectSelector.nFatJetMin)
            print(f'\n selection.all("FatJetPt") ({type(selection.all("FatJetPt"))}) ({len(selection.all("FatJetPt"))}): {selection.all("FatJetPt")}')

            printVariable("\n\n events.FatJet.eta", events.FatJet.eta)
            printVariable("\n mask_FatJetEta", mask_FatJetEta)
            printVariable("\n events.FatJet[mask_FatJetEta].eta", events.FatJet[mask_FatJetEta].eta)
            printVariable("\n ak.num(events.FatJet[mask_FatJetEta]) >= self.objectSelector.nFatJetMin", ak.num(events.FatJet[mask_FatJetEta]) >= self.objectSelector.nFatJetMin)
            #printVariable('\n selection.all("FatJetEta")', selection.all("FatJetEta"))            
            print(f'\n selection.all("FatJetEta") ({type(selection.all("FatJetEta"))}) ({len(selection.all("FatJetEta"))}): {selection.all("FatJetEta")}')
            

            printVariable("\n events.FatJet.btagDeepB", events.FatJet.btagDeepB)
            printVariable("\n mask_FatJetBtagDeepB", mask_FatJetBtagDeepB)
            printVariable("\n events.FatJet[mask_FatJetBtagDeepB].btagDeepB", events.FatJet[mask_FatJetBtagDeepB].btagDeepB)
            printVariable("\n ak.num(events.FatJet[mask_FatJetBtagDeepB])", ak.num(events.FatJet[mask_FatJetBtagDeepB]))
            printVariable("\n ak.num(events.FatJet[mask_FatJetBtagDeepB]) >= self.objectSelector.nFatJetMin ", ak.num(events.FatJet[mask_FatJetBtagDeepB]) >= self.objectSelector.nFatJetMin )
            print(f'\n selection.all("FatJetBtagDeepB") ({type(selection.all("FatJetBtagDeepB"))}) ({len(selection.all("FatJetBtagDeepB"))}): {selection.all("FatJetBtagDeepB")}')

            printVariable("\n events.FatJet.msoftdrop", events.FatJet.msoftdrop)
            printVariable("\n mask_FatJetMSoftDrop", mask_FatJetMSoftDrop)
            printVariable("\n events.FatJet[mask_FatJetMSoftDrop].msoftdrop", events.FatJet[mask_FatJetMSoftDrop].msoftdrop)
            printVariable("\n ak.num(events.FatJet[mask_FatJetMSoftDrop])", ak.num(events.FatJet[mask_FatJetMSoftDrop]))
            printVariable("\n ak.num(events.FatJet[mask_FatJetMSoftDrop]) >= self.objectSelector.nFatJetMin ", ak.num(events.FatJet[mask_FatJetMSoftDrop]) >= self.objectSelector.nFatJetMin )
            print(f'\n selection.all("FatJetMSoftDrop") ({type(selection.all("FatJetMSoftDrop"))}) ({len(selection.all("FatJetMSoftDrop"))}): {selection.all("FatJetMSoftDrop")}')


            print(f"\n\nsel_names_all: {sel_names_all}")
            print(f'\n selection.all(* sel_names_all["SR"]) ({type(selection.all(* sel_names_all["SR"]))}) ({len(selection.all(* sel_names_all["SR"]))}): {selection.all(* sel_names_all["SR"])}')
            
            print(f'\n selection.all(* sel_names_all["GenHToAATo4B"]) ({type(selection.all(* sel_names_all["GenHToAATo4B"]))}) ({len(selection.all(* sel_names_all["GenHToAATo4B"]))}): {selection.all(* sel_names_all["GenHToAATo4B"])}')
            
            printVariable("\n events.FatJet.pt[sel_SR]", events.FatJet.pt[sel_SR])
            printVariable("\n selFatJet.pt", selFatJet.pt)
            printVariable("\n selFatJet.pt[sel_SR]", selFatJet.pt[sel_SR])

            printVariable("\n selFatJet.pt[mask_FatJetPt]", selFatJet.pt[mask_FatJetPt])

            printVariable("\n mask_FatJetPt", mask_FatJetPt)
            printVariable("\n mask_FatJetEta", mask_FatJetEta)
            printVariable("\n mask_FatJetBtagDeepB", mask_FatJetBtagDeepB)
            printVariable("\n mask_FatJetMSoftDrop", mask_FatJetMSoftDrop)
            printVariable("\n mask_all &:", (
                mask_FatJetPt &
                mask_FatJetEta &
                mask_FatJetBtagDeepB &
                mask_FatJetMSoftDrop
            ))
            #printVariable("\n ", )

            
            
            
            '''
            printVariable("events.FatJet.btagDeepB", events.FatJet.btagDeepB)
            printVariable("", )
            printVariable("", )

            printVariable("mask_FatJetBtagDeepB", mask_FatJetBtagDeepB)
            printVariable("", )
            '''

                    
        


            
            


        
        ################
        # EVENT WEIGHTS
        ################

        # create a processor Weights object, with the same length as the number of events in the chunk
        weights     = Weights(len(events))
        weights_gen = Weights(len(events))
        weights_GenHToAATo4B = Weights(len(events))


        if self.datasetInfo[dataset]["isMC"]:
            weights.add(
                "lumiWeight",
                weight=np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
            )
            weights.add(
                "genWeight",
                weight=np.copysign(np.ones(len(events)), events.genWeight)
            )
            weights.add(
                "btagWeight",
                weight=(events.btagWeight.DeepCSVB)
            )
            
 
            weights_gen.add(
                "lumiWeight",
                weight=np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
            )
            weights_gen.add(
                "genWeight",
                weight=np.copysign(np.ones(len(events)), events.genWeight)
            )
            
            
            if printLevel >= 5:
                print(f"\nevents.genWeight ({events.genWeight.fields}) ({len(events.genWeight)}): {events.genWeight.to_list()}")
                genWgt1 = np.copysign(np.ones(len(events)), events.genWeight)
                print(f"genWgt1 ({len(genWgt1)}): {genWgt1}")
                
                print(f"\n\nevents.btagWeight ({events.btagWeight.fields}) ({len(events.btagWeight)}): {events.btagWeight.to_list()}")
                #print(f"events.nPSWeight {events.nPSWeight.to_list()}")
                print(f"\n\nevents.PSWeight ({events.PSWeight.fields}) ({len(events.PSWeight)}): {events.PSWeight.to_list()}")
                print(f"\n\nevents.LHEWeight ({events.LHEWeight.fields}) ({len(events.LHEWeight)}): {events.LHEWeight.to_list()}")
                print(f"\n\nevents.LHEScaleWeight ({events.LHEScaleWeight.fields}) ({len(events.LHEScaleWeight)}): {events.LHEScaleWeight.to_list()}")
                print(f"\n\nevents.LHEReweightingWeight ({events.LHEReweightingWeight.fields}) ({len(events.LHEReweightingWeight)}): {events.LHEReweightingWeight.to_list()}")

                

            if printLevel >= 13:                
                printVariable("\n weights.weight()", weights.weight())
                printVariable("\n weights_gen.weight()", weights_gen.weight())
                printVariable("\n weights_GenHToAATo4B.weight()", weights_GenHToAATo4B.weight())
                #printVariable("\n weights.weight().sum()", weights.weight().sum())
                print(f"weights.weight().sum(): {weights.weight().sum()}")
                printVariable("\n weights.weight()[sel_SR]", weights.weight()[sel_SR])
                print(f"weights.weight()[sel_SR].sum(): {weights.weight()[sel_SR].sum()}")

            if printLevel >= 30:
                print(f"leadingFatJet.tau4[sel_SR]): {leadingFatJet.tau4[sel_SR]}")
                print(f"leadingFatJet.tau3[sel_SR]): {leadingFatJet.tau3[sel_SR]}")
                print(f"divide                     : {np.divide(leadingFatJet.tau4[sel_SR], leadingFatJet.tau3[sel_SR])} \n\n\n")
        
        ###################
        # FILL HISTOGRAMS
        ###################

        systList = []
        if self.datasetInfo[dataset]['isMC']:
            if shift_syst is None:
                systList = [
                    "central"
                ]
            else:
                systList = [shift_syst]
        else:
            systList = ["noweight"]

            
        output['cutflow']['all events'] += len(events)
        output['cutflow'][sWeighted+'all events'] += weights.weight().sum() 
        #for n in selection.names:
        #    output['cutflow'][n] += selection.all(n).sum()

        for iSelection in sel_names_all.keys():
            iName = f"{iSelection}: {sel_names_all[iSelection]}"
            sel_i = selection.all(* sel_names_all[iSelection])
            output['cutflow'][iName] += sel_i.sum()
            output['cutflow'][sWeighted+iName] +=  weights.weight()[sel_i].sum()
            

        for syst in systList:

            # find the event weight to be used when filling the histograms
            weightSyst = syst
            
            # in the case of 'central', or the jet energy systematics, no weight systematic variation is used (weightSyst=None)
            if syst in ["central", "JERUp", "JERDown", "JESUp", "JESDown"]:
                weightSyst = None

            
            
            if syst == "noweight":
                evtWeight = np.ones(len(events))
            else:
                evtWeight     = weights.weight(weightSyst)
                evtWeight_gen = weights_gen.weight(weightSyst)


            if printLevel >= 15:
                print(f"ak.num(selFatJet[sel_SR]) ({len(ak.num(selFatJet[sel_SR]))}): {ak.num(selFatJet[sel_SR])}")
                print(f"ak.to_numpy(ak.num(selFatJet[sel_SR])) ({len(ak.to_numpy(ak.num(selFatJet[sel_SR])))}): {ak.to_numpy(ak.num(selFatJet[sel_SR]))}")
                print(f"selFatJet.pt[sel_SR] ({len(selFatJet.pt[sel_SR])}): {selFatJet.pt[sel_SR]}")
                print(f"ak.flatten(selFatJet.pt[sel_SR]) ({len(ak.flatten(selFatJet.pt[sel_SR]))}): {ak.flatten(selFatJet.pt[sel_SR])}")
                print(f"selFatJet.pt[sel_SR][:, 0] ({len(selFatJet.pt[sel_SR][:, 0])}): {selFatJet.pt[sel_SR][:, 0]}")
                #print(f"ak.flatten(selFatJet.pt[sel_SR][:, 0]) : {ak.flatten(selFatJet.pt[sel_SR][:, 0])}")
                #print(f"ak.flatten(selFatJet.pt[sel_SR][:, 0]) ({len(ak.flatten(selFatJet.pt[sel_SR][:, 0]))}): {ak.flatten(selFatJet.pt[sel_SR][:, 0])}")

                print(f"evtWeight ({len(evtWeight)}): {evtWeight}")
                print(f"evtWeight[sel_SR] ({len(evtWeight[sel_SR])}): {evtWeight[sel_SR]}")

            if printLevel >= 30:
                printVariable("\n weights.weight()", weights.weight())
                printVariable("\n weights.weight()[sel_SR]", weights.weight()[sel_SR])
                printVariable("\n selFatJet.pt", selFatJet.pt)
                printVariable("\n selFatJet.pt[sel_SR]", selFatJet.pt[sel_SR])
                printVariable("\n leadingFatJet.pt", leadingFatJet.pt)
                printVariable("\n leadingFatJet.pt[sel_SR]", leadingFatJet.pt[sel_SR])

            if printLevel >= 30:
                printVariable("\n ones_list", ones_list)
                printVariable("\n ones_list[sel_SR]", ones_list[sel_SR])

                printVariable("\n evtWeight", evtWeight)
                printVariable("\n evtWeight[sel_SR]", evtWeight[sel_SR])

                k = 3
                print(f"k: {k}")
                printVariable("\n ones_list * k", ones_list * k)
                printVariable("\n ones_list[sel_SR] * k", ones_list[sel_SR] * k)
                
                printVariable("\n evtWeight * k", evtWeight * k)
                printVariable("\n evtWeight[sel_SR]", evtWeight[sel_SR] * k)
                

            # all events
            iBin = 0
            output['hCutFlow'].fill(
                dataset=dataset,
                CutFlow=(ones_list * iBin),
                systematic=syst
            )
            output['hCutFlowWeighted'].fill(
                dataset=dataset,
                CutFlow=(ones_list * iBin),
                systematic=syst,
                weight=evtWeight
            )

            # QCD MC ----------------------------------------------
            if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isQCD'] :
                # genBHadrons_status2 events
                iBin = 1
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_genBHadrons_status2_eventwise] * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_genBHadrons_status2_eventwise] * iBin),
                    systematic=syst,
                    weight=evtWeight
                )

                # genBHadrons_status2 events
                iBin = 2
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_genBQuarks_hardSctred_eventwise] * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_genBQuarks_hardSctred_eventwise] * iBin),
                    systematic=syst,
                    weight=evtWeight
                )

                # QCD_stitch events
                iBin = 3
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_QCD_stitch_eventwise] * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_QCD_stitch_eventwise] * iBin),
                    systematic=syst,
                    weight=evtWeight
                )

                

            # events passing SR
            iBin = 4
            output['hCutFlow'].fill(
                dataset=dataset,
                CutFlow=(ones_list[sel_SR] * iBin),
                systematic=syst
            )
            output['hCutFlowWeighted'].fill(
                dataset=dataset,
                CutFlow=(ones_list[sel_SR] * iBin),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )

            
            
            '''
            output['nSelFatJet'].fill(
                dataset=dataset,
                nObject=ak.to_numpy(ak.num(selFatJet[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            '''
            
            output['hLeadingFatJetPt'].fill(
                dataset=dataset,
                #Pt=ak.flatten(selFatJet.pt[sel_SR][:, 0]),
                #Pt=(selFatJet.pt[sel_SR][:, 0]),
                Pt=(leadingFatJet.pt[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            
            output['hLeadingFatJetEta'].fill(
                dataset=dataset,
                Eta=(leadingFatJet.eta[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetPhi'].fill(
                dataset=dataset,
                Phi=(leadingFatJet.phi[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetMass'].fill(
                dataset=dataset,
                Mass=(leadingFatJet.mass[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetMSoftDrop'].fill(
                dataset=dataset,
                Mass=(leadingFatJet.msoftdrop[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetBtagDeepB'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.btagDeepB[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetBtagDDBvLV2'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.btagDDBvLV2[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetBtagDDCvBV2'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.btagDDCvBV2[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetBtagHbb'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.btagHbb[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_H4qvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_HbbvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )            
            output['hLeadingFatJetDeepTagMD_ZHbbvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_ZHccvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ZHccvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_ZbbvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ZbbvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_ZvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ZvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_bbvsLight'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_ccvsLight'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ccvsLight[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTag_H'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTag_H[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTag_QCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTag_QCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTag_QCDothers'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTag_QCDothers[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )

            
            output['hLeadingFatJetN2b1'].fill(
                dataset=dataset,
                N2=(leadingFatJet.n2b1[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetN3b1'].fill(
                dataset=dataset,
                N3=(leadingFatJet.n3b1[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            
            output['hLeadingFatJetTau1'].fill(
                dataset=dataset,
                TauN=(leadingFatJet.tau1[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau2'].fill(
                dataset=dataset,
                TauN=(leadingFatJet.tau2[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau3'].fill(
                dataset=dataset,
                TauN=(leadingFatJet.tau3[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau4'].fill(
                dataset=dataset,
                TauN=(leadingFatJet.tau4[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )

            output['hLeadingFatJetTau4by3'].fill(
                dataset=dataset,
                TauN=(np.divide(leadingFatJet.tau4[sel_SR], leadingFatJet.tau3[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau3by2'].fill(
                dataset=dataset,
                TauN=(np.divide(leadingFatJet.tau3[sel_SR], leadingFatJet.tau2[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau2by1'].fill(
                dataset=dataset,
                TauN=(np.divide(leadingFatJet.tau2[sel_SR], leadingFatJet.tau1[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            
            
            output['hLeadingFatJetNConstituents'].fill(
                dataset=dataset,
                nObject=(leadingFatJet.nConstituents[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            if self.datasetInfo[dataset]['isMC']:
                output['hLeadingFatJetNBHadrons'].fill(
                    dataset=dataset,
                    nObject=(leadingFatJet.nBHadrons[sel_SR]),
                    systematic=syst,
                    weight=evtWeight[sel_SR]
                )
                output['hLeadingFatJetNCHadrons'].fill(
                    dataset=dataset,
                    nObject=(leadingFatJet.nCHadrons[sel_SR]),
                    systematic=syst,
                    weight=evtWeight[sel_SR]
                )            
                
            
            output['hLeadingFatJetParticleNetMD_QCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNetMD_QCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNetMD_Xbb'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNetMD_Xbb[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNetMD_Xcc'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNetMD_Xcc[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNetMD_Xqq'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNetMD_Xqq[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNet_H4qvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNet_H4qvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNet_HbbvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNet_HbbvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )

            output['hLeadingFatJetParticleNet_HccvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNet_HccvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNet_QCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNet_QCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNet_mass'].fill(
                dataset=dataset,
                Mass=(leadingFatJet.particleNet_mass[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )


            if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isSignal']: 
                output['hGenHiggsPt_all'].fill(
                    dataset=dataset,
                    Pt=(ak.flatten(genHiggs.pt[sel_GenHToAATo4B])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hGenHiggsPt_sel'].fill(
                    dataset=dataset,
                    Pt=(ak.flatten(genHiggs.pt[sel_SR])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_SR]
                )
                output['hGenHiggsPt_sel_wGenCuts'].fill(
                    dataset=dataset,
                    Pt=(ak.flatten(genHiggs.pt[sel_GenHToAATo4B])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                
                # m(2b from ATo2B) and m(4b from HToAATo4b) --------------                   
                output['hGenHiggsMass_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                    
                output['hMass_GenA_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 0].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hMass_GenA_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                '''
                output['hMass_GenA_all'].fill(
                    dataset=dataset,
                    Pt=(ak.flatten(genA_Second.mass)),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                '''
                output['hMass_GenAApair_all'].fill(
                    dataset=dataset,
                    Mass=((genA_First + genA_Second).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hMass_GenAToBBbarpair_all'].fill(
                    dataset=dataset,
                    Mass=((events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 0][sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )                
                output['hMass_GenAToBBbarpair_all'].fill(
                    dataset=dataset,
                    Mass=((events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 1][sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hMass_Gen4BFromHToAA_all'].fill(
                    dataset=dataset,
                    Mass=((events.GenPart[genBBar_pairs['b']][:, 0] + events.GenPart[genBBar_pairs['bbar']][:, 0] + events.GenPart[genBBar_pairs['b']][:, 1] + events.GenPart[genBBar_pairs['bbar']][:, 1]).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )


                output['hMass_GenAToBBbarpair_all_1'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_0 + LVGenBbar_0).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )                
                output['hMass_GenAToBBbarpair_all_1'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_1 + LVGenBbar_1).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hMass_Gen4BFromHToAA_all_1'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_0 + LVGenBbar_0 + LVGenB_1 + LVGenBbar_1).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )


                output['hMass_GenA1_vs_GenA2_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 0].mass[sel_GenHToAATo4B]),
                    Mass1=(genACollection[:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                
                output['hMass_GenA1ToBBbar_vs_GenA2ToBBbar_all'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_0 + LVGenBbar_0).mass[sel_GenHToAATo4B]),
                    Mass1=((LVGenB_1 + LVGenBbar_1).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMass_GenH_vs_GenAHeavy_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    Mass1=(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                
                output['hMass_GenH_vs_GenALight_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    Mass1=(genACollection[idxGenA_sortByMass][:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMass_GenAHeavy_vs_GenALight_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B]),
                    Mass1=(genACollection[idxGenA_sortByMass][:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hDeltaR_GenH_GenB_max'].fill(
                    dataset=dataset,
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMassGenH_vs_maxDRGenHGenB_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMassGenAHeavy_vs_maxDRGenHGenB_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B]),
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMassGenALight_vs_maxDRGenHGenB_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[idxGenA_sortByMass][:, 1].mass[sel_GenHToAATo4B]),
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                
            # QCD MC ----------------------------------------------
            if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isQCD'] :
                output['hGenLHE_HT_all'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_HT_QCDStitch'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_QCD_stitch_eventwise]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_QCD_stitch_eventwise]
                )
                
                output['hGenBquark_Status_all'].fill(
                    dataset=dataset,
                    PytPartStatus=(ak.flatten(genBQuarks_QCD.status)),
                    systematic=syst
                )
                output['hGenBquark_first_Status_all'].fill(
                    dataset=dataset,
                    PytPartStatus=(ak.to_numpy(genBQuarks_first[mask_genBQuarks].status)),
                    systematic=syst
                )
                output['hGenBquark_first_PdgId_all'].fill(
                    dataset=dataset,
                    PdgId=(ak.to_numpy(abs(genBQuarks_first[mask_genBQuarks].pdgId))),
                    systematic=syst
                )

                for statusFlag_ in GENPART_STATUSFLAGS:
                    output['hGenBquark_first_%s_all' % (statusFlag_)].fill(
                        dataset=dataset,
                        Boolean=(ak.to_numpy(genBQuarks_first[mask_genBQuarks].hasFlags(statusFlag_))),
                        systematic=syst
                    )

                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 1)
                output['hGenBquark_leadingPt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[mask_tmp][:, 0]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 2)
                output['hGenBquark_subleadingPt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[mask_tmp][:, 1]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 3)
                output['hGenBquark_thirdLeadingPt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[mask_tmp][:, 2]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 4)
                output['hGenBquark_forthLeadingPt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[mask_tmp][:, 3]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 1)
                output['hGenBquark_leadingPt_QCDStitch'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_eventwise)][:, 0]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_eventwise)]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 2)
                output['hGenBquark_subleadingPt_QCDStitch'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_eventwise)][:, 1]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_eventwise)]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 3)
                output['hGenBquark_thirdLeadingPt_QCDStitch'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_eventwise)][:, 2]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_eventwise)]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 4)
                output['hGenBquark_forthLeadingPt_QCDStitch'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_eventwise)][:, 3]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_eventwise)]
                )

                
                if printLevel >= 100:                    
                    printVariable("\n genBQuarks_pT", genBQuarks_pT)
                    printVariable("\n ak.count(genBQuarks_pT, axis=-1)", ak.count(genBQuarks_pT, axis=-1))
                    printVariable("\n (ak.count(genBQuarks_pT, axis=-1) >= 2)", (ak.count(genBQuarks_pT, axis=-1) >= 2))
                    printVariable("\n evtWeight_gen", evtWeight_gen)
                    printVariable("\n evtWeight_gen[(ak.count(genBQuarks_pT, axis=-1) >= 2)]", evtWeight_gen[(ak.count(genBQuarks_pT, axis=-1) >= 2)])
                    printVariable("\n genBQuarks_pT[(ak.count(genBQuarks_pT, axis=-1) >= 2)]", genBQuarks_pT[(ak.count(genBQuarks_pT, axis=-1) >= 2)])
                    printVariable("\n genBQuarks_pT[(ak.count(genBQuarks_pT, axis=-1) >= 2)][:, 1]", genBQuarks_pT[(ak.count(genBQuarks_pT, axis=-1) >= 2)][:, 1])

                '''
                Notes: Check status of first occuriing b-quark in event, and check how it is produced.

                '''
                    


                
            '''
            output[''].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output[''].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output[''].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            '''

            

            '''
            output[''].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            '''

        
        return output


    def postprocess(self, accumulator):
        #pass
        return accumulator



    
def getLorentzVectorFromAwkArray(awkArray, ptObjectName, etaObjectName, phiObjectName, massObjectName):
    '''
    return ak.Array(
        awkArray,
        with_name='Momentum4D'
    )
    '''
    v1 = ak.zip( {
        'pt':   awkArray[ptObjectName],
        'eta':  awkArray[etaObjectName],
        'phi':  awkArray[phiObjectName],
        'mass': awkArray[massObjectName],
    })
    if printLevel >= 12:
        print(f"getLorentzVectorFromAwkArray(): v1 ({type(v1)}): {v1.to_list()}")

    return ak.Array(v1, with_name='Momentum4D')
                       

def printWithType(sX, X):
    #print(f"{sX} ({type(X)}): {X}")
    print(f"{sX} : {X}")


        
    
    
if __name__ == '__main__':
    print("htoaa_Analysis:: main: {}".format(sys.argv)); sys.stdout.flush()


    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    
    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))

    lumiScale = 1
    sInputFiles         = config["inputFiles"]
    sOutputFile         = config["outputFile"]
    sample_category     = config['sampleCategory']
    isMC                = config["isMC"]
    era                 = config['era']
    if isMC:
        luminosity          = Luminosities[era][0]
        sample_crossSection = config["crossSection"]
        sample_nEvents      = config["nEvents"]
        sample_sumEvents    = config["sumEvents"] if config["sumEvents"] != -1 else sample_nEvents
        if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
        lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)
    #branchesToRead = htoaa_nanoAODBranchesToRead
    #print("branchesToRead: {}".format(branchesToRead))

    print(f"isMC: {isMC}, lumiScale: {lumiScale}")
    sInputFiles_toUse = []
    for sInputFile in sInputFiles:
        if "*" in sInputFile:  sInputFiles_toUse.extend( glob.glob( sInputFile ) )
        else:                  sInputFiles_toUse.append( sInputFile )
    sInputFiles = sInputFiles_toUse
    for iFile in range(len(sInputFiles)):        
        if sInputFiles[iFile].startswith("/store/"): # LFN: Logical File Name
            #sInputFiles[iFile] = xrootd_redirectorName + sInputFiles[iFile]
            sInputFiles[iFile] = setXRootDRedirector(sInputFiles[iFile])
    print(f"sInputFiles ({len(sInputFiles)}) (type {type(sInputFiles)}):");
    for sInputFile in sInputFiles:
        print(f"\t{sInputFile}")
    sys.stdout.flush()


    startTime = time.time()
    
    tracemalloc.start()
    



    #client = Client("tls://localhost:8786")
    #executor = processor.DaskExecutor(client=client)
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    nWorkers= 4 if chunksize > 1 else 1
    print(f"nEventsToAnalyze: {nEventsToAnalyze},  nEventToReadInBatch: {nEventToReadInBatch}, chunksize: {chunksize},  maxchunks: {maxchunks},  nWorkers: {nWorkers}")
    run = processor.Runner(
        #executor=executor,
        executor=processor.FuturesExecutor(workers=nWorkers),
        schema=schemas.NanoAODSchema,
        savemetrics=True,
        chunksize=chunksize,  #3 ** 20,  ## Governs the number of times LeptonJetProcessor "process" is called
        maxchunks=maxchunks
    )

    output, metrics = run(
        fileset={sample_category: sInputFiles},
        #fileset={"QCD": ["/home/siddhesh/Work/CMS/htoaa/analysis/tmp/20BE2B12-EFF6-8645-AB7F-AFF6A624F816.root"]},
        treename="Events",
        processor_instance=HToAATo4bProcessor(
            datasetInfo={
                "era": era, 
                sample_category: {"isMC": isMC, "lumiScale": lumiScale}
            }
        )
    )

    print(f"metrics: {metrics}")


    if 'cutflow' in output.keys():
        print("Cutflow::")
        #for key, value in output['cutflow'].items():
        for key in output['cutflow'].keys():
            #print(key, value)
            if key.startswith(sWeighted): continue

            print("%10f\t%10d\t%s" % (output['cutflow'][sWeighted+key], output['cutflow'][key], key))


    
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        #sOutputFile = sOutputFile.replace('.root', '_wCoffea.root') # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        sDir1 = 'evt/%s' % (sample_category)

        
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                #print(f"key: {key},  value ({type(value)}): {value}")
                #if not (key.startswith('h') or key != 'cutflow'): continue
                if not isinstance(value, hist.Hist): continue
                #print(f"1: key {key}, value ({type(value)})     Hist: {type(hist.Hist)},    isinstance(value, hist.Hist): {isinstance(value, hist.Hist)}") # value: {value}")

                '''
                print(f"value.DEFAULT_DTYPE {value.DEFAULT_DTYPE}")
                print(f"value.fields ({type(value.fields)}): {value.fields}")
                print(f"value.label ({type(value.label)}): {value.label}")
                print(f"value.axes() ({type(value.axes())}): {value.axes()}")
                print(f"value.sparse_axes() ({type(value.sparse_axes())}): {value.sparse_axes()}")
                print(f"value.sparse_dim() {value.sparse_dim()}")
                for ax_label in value.fields:
                    print(f"ax_label ({type(ax_label)}): {ax_label}")
                    print(f"value.axis(ax_label) ({type(value.axis(ax_label))}): {value.axis(ax_label)}")
                    #printWithType('value.axis(ax_label).identifiers()', value.axis(ax_label).identifiers())
                    print(f"value.axis(ax_label).identifiers()  : {value.axis(ax_label).identifiers()}")

                for ax in value.sparse_axes():
                    printWithType('ax', ax)

                h1 = value.to_hist()
                print(f"h1 ({type(h1)}): h1")
                '''


                for _dataset in value.axis('dataset').identifiers():
                    #print(f"_dataset ({type(_dataset)}): {_dataset}")

                    for _syst in value.axis('systematic').identifiers():
                        #print(f"_syst ({type(_syst)}): {_syst}")

                        h1 = value.integrate('dataset',_dataset).integrate('systematic',_syst).to_hist()
                        #print(f"h1 ({type(h1)}): h1")

                        #fOut['%s/%s_%s_%s' % (sDir1, key, _dataset, _syst)] = h1
                        fOut['%s/%s_%s' % (sDir1, key, _syst)] = h1
                
                #fOut['%s%s' % (sDir1, key)] = value
                #fOut['%s%s' % (sDir1, key)] = hist.export1d(value)
                #fOut['%s%s' % (sDir1, key)] = h1

                
        

        #util.save(output, sOutputFile)
            
        
        print("Wrote to sOutputFile {}".format(sOutputFile))
   















    
    current_memory, peak_memory = tracemalloc.get_traced_memory() # https://medium.com/survata-engineering-blog/monitoring-memory-usage-of-a-running-python-program-49f027e3d1ba
    print(f"\n\nMemory usage:: current {current_memory / 10**6}MB;  peak {peak_memory / 10**6}MB")

    endTime = time.time()
    totalTime = endTime - startTime
    totalTime_hr  = int(totalTime/60/60)
    totalTime_min = totalTime - float(totalTime_hr * 60)
    totalTime_min = int(totalTime_min/60)
    totalTime_sec = totalTime - float(totalTime_hr * 60*60) - float(totalTime_min * 60)
    print(f"Total run time: {totalTime_hr}h {totalTime_min}m {totalTime_sec}s = {totalTime}sec ")
    