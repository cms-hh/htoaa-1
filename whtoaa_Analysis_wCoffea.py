#https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=HIG-18-026&tp=an&id=2158&ancode=HIG-18-026
#htoaa analysis main code
#/afs/cern.ch/work/s/snandan/anaconda3/envs/myEnv/lib/python3.10/site-packages/boost_histogram/_internal/hist.py
#https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Jet_Energy_Corrections_in_Run2
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution
#https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Description
#https://cms-opendata-guide.web.cern.ch/analysis/systematics/objectsuncertain/jetmetuncertain/
#https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
#https://twiki.cern.ch/twiki/bin/view/CMS/JetEnergyScale#Jet_Energy_Corrections_for_CMS_u
#https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC
import os
import sys
import json
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
import numpy as np
from copy import deepcopy
import uproot
from auxiliary import met_p4
from auxiliary import selectRunLuminosityBlock
from metfilter import apply_metfilters
np.set_printoptions(threshold=sys.maxsize)
import hist
from histograms import histograms, clip_value
from weights import event_weights#calculate_weight, transfer_factor
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
#from coffea import hist
import awkward as ak
import uproot


from htoaa_Settings import *#bTagWPs
from htoaa_CommonTools import GetDictFromJsonFile, calculate_lumiScale, setXRootDRedirector, xrdcpFile, get_lf
from objectselection import ObjectSelection
from genparticle import genparticle

nEventToReadInBatch = 0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze =  -1 #-1 # 1000 # 100000 # -1
sWeighted = "Wtd: "

class HToAATo4bProcessor(processor.ProcessorABC):
    def __init__(self, datasetInfo={}):

        ak.behavior.update(nanoaod.behavior)
        self.config = datasetInfo
        self.objectSelector = ObjectSelection(era=self.config["era"])

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata["dataset"] # dataset label
        if not self.config['isMC']:
            with open('Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt') as f:
                rl = json.load(f)
            events = events[selectRunLuminosityBlock(rl, events.run, events.luminosityBlock)]
        events = apply_metfilters(self.config['isMC'], events)
        # Gen-level selection ---------------------------------------------------------------------
        genHiggs = None
        genHT    = None
        
        if self.config['isSignal']:
            genHiggs  = self.objectSelector.selectGenHiggs(events)        
            genHT     = self.objectSelector.GenHT(events)
            # m(bbar from A) and m(4b from HToAA) ----------

            genACollection = self.objectSelector.selectGenABoson(events)
            genA_First  = genACollection[:, 0]
            genA_Second = genACollection[:, 1]
            
            idxGenA_sortByMass = ak.argsort(genACollection.mass, axis=-1, ascending=False)
            
        # QCD MC ----------------------------------------------
        if self.config['isQCD'] :
            genBQuarks_QCD = events.GenPart[(
                (abs(events.GenPart.pdgId) == 5)
            )]

        # Reco-level -----------------------------------------------------------------------------------
        
        ##################
        # EVENT VARIABLES
        ##################

        FatJet = self.objectSelector.selectFatJets(events)
        hto4b_vs_QCDsum = FatJet.particleNetMD_Hto4b_Haa4b/(FatJet.particleNetMD_Hto4b_Haa4b+FatJet.particleNetMD_Hto4b_QCD0b+FatJet.particleNetMD_Hto4b_QCD1b+FatJet.particleNetMD_Hto4b_QCD2b+FatJet.particleNetMD_Hto4b_QCD3b+FatJet.particleNetMD_Hto4b_QCD4b)
        FatJet['particleNetMD_Hto4b_vs_QCDsum'] = hto4b_vs_QCDsum
        FatJet = FatJet[ak.argsort(FatJet.particleNetMD_Hto4b_vs_QCDsum, axis=-1, ascending=False)]
        ak4Jet = self.objectSelector.selectak4Jets(events)
        preselele = self.objectSelector.selectElectrons(events)
        preselmu = self.objectSelector.selectMuons(events)
        vetolep = ak.with_name(
            ak.concatenate([preselele[(preselele.mvaTTH >= 0.3)], preselmu[(preselmu.mvaTTH >= 0.5)]], axis=1), "PtEtaPhiMCandidate"
            )
        vetolep = vetolep[ak.argsort(vetolep.pt, axis=-1, ascending=False)]
        elesize = ak.sum(ak.num(preselele))
        musize = ak.sum(ak.num(preselmu))
        is_triggered_1ele = events.HLT.Ele28_eta2p1_WPTight_Gsf_HT150 | events.HLT.Ele32_WPTight_Gsf
        is_triggered_1mu = events.HLT.IsoMu24 | events.HLT.IsoMu27
        sel_triggered_1ele = (self.config['use_triggers_1e']) & is_triggered_1ele
        sel_triggered_1mu = (self.config['use_triggers_1mu']) & is_triggered_1mu
        preselele = ak.with_field(preselele, ((preselele.mvaTTH >=0.3) & (preselele.pt > 33) & (preselele.miniPFRelIso_all < 0.15)), "tight_ele")
        preselmu = ak.with_field(preselmu, ((preselmu.mvaTTH >=0.5) & (preselmu.pt > 25) & (preselmu.miniPFRelIso_all < 0.15)), "tight_mu")
        preselele = ak.with_field(preselele, (preselele.mvaTTH <0.3), "fake_ele")
        preselmu = ak.with_field(preselmu, (preselmu.mvaTTH <0.5), "fake_mu")

        #https://github.com/HEP-KBFI/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L94 
        trig_obj_ele = (abs(events.TrigObj.id) == 11) & ( (events.TrigObj.filterBits & 2) == 2 )
        trig_obj_mu = (abs(events.TrigObj.id) == 13) & ( (events.TrigObj.filterBits & 8) == 8 )

        if self.config['isSignal']:
            events['GenPart'] = ak.with_field(events.GenPart,
                        (abs(events.GenPart.pdgId) == 5)
                      & (events.GenPart[events.GenPart.genPartIdxMother].pdgId == 36),
                        'bfroma'
              )
            events['GenPart'] = ak.with_field(events.GenPart,
                        (abs(events.GenPart.pdgId) == 36),
                        'aboson'
              )
            FatJet = FatJet[ak.all(FatJet.metric_table(events.GenPart[events.GenPart.bfroma]) <0.8, axis=-1)]
            '''selection.add("1GenHiggs", ak.num(genHiggs) == 1)
            selection.add("2GenA", ak.num(genACollection) == 2)
            selection.add("2GenAToBBbarPairs", ak.num(genBBar_pairs) == 2)'''
        FatJet = ak.with_field(FatJet, FatJet.msoftdrop >= 90, "msoftdrop_GE_90")
        FatJet = ak.with_field(FatJet, FatJet.msoftdrop < 90, "msoftdrop_L_90")
        evt_weights = event_weights(events, self.config)
        weights = evt_weights.calculate_weight()

        systList = ["central"]
        '''if self.config['isMC']:
        if dataset.startswith('TT'):
        systList.extend(["topt_reweigingUp", "topt_reweigingDown"])'''

        output = histograms(systList, self.config['msoftdrop'])
        output['cutflow'] = {}
        output['cutflow_weighted'] = {}
        events['GenPart'] = ak.with_field(events.GenPart, ak.local_index(events.GenPart, axis=1), 'idx')
        events['GenPart'] = ak.with_field(events.GenPart,
                        (abs(events.GenPart.pdgId)==5)
                      & (abs(events.GenPart[events.GenPart.genPartIdxMother].pdgId) == 6),
                      'bfromtop'
        )
        events['GenPart'] = ak.with_field(events.GenPart,
                        (abs(events.GenPart.pdgId)==24)
                        & (abs(events.GenPart[events.GenPart.genPartIdxMother].pdgId) == 6),
                        'wfromtop'
        )
        events['GenPart'] = ak.with_field(events.GenPart,
                        (abs(events.GenPart.pdgId) < 6)
                    & (abs(events.GenPart[events.GenPart.genPartIdxMother].pdgId) == 24),
                        'qfromW'
        )
        if self.config['lepton_selection'] == 'tight_lep':
            selele = preselele[preselele.tight_ele]
            selmu = preselmu[preselmu.tight_mu]
            sellep_all = ak.with_name(
                ak.concatenate([selmu, selele], axis=1), "PtEtaPhiMCandidate"
            )
        else:
            if 1:
                sellep_all = preselele[preselele.fake_ele]
                sellep_tight = preselele[preselele.tight_ele]
            else:
                sellep_all = preselmu[preselmu.fake_mu]
                sellep_tight = preselmu[preselmu.tight_mu]
        sellep_all = sellep_all[ak.argsort(sellep_all.pt, axis=-1, ascending=False)]
        sellep = sellep_all[:,0:1]
        vetolep = vetolep[ak.all(vetolep.metric_table(sellep) > 0.3, axis=-1)]
        '''print('******** ', sellep.pdgId[1910:1920].to_list())
        print('********veto ', vetolep[1910:1920].to_list())
        print('********pdg ', ak.num(vetolep.pdgId)[1910:1920].to_list())
        v,s = ak.unzip(ak.cartesian([vetolep, sellep], axis=1, nested=True))
        print('****v: ', v.pdgId[1910:1920])
        vetoidx = vetolep.metric_table(sellep, metric=lambda a,b: (abs(a.pdgId) == abs(b.pdgId)) & (a.idx == b.idx))
        print('****ve: ', vetoidx[1910:1920].to_list())
        print('****s: ', sellep[1910:1920].to_list())
        vmasked = v.mask[~vetoidx]
        #print('***8v: ', vmasked.fill_none[1910:1920].to_list())
        vetolep = vetolep[~ak.all(vetolep.metric_table(sellep, metric=lambda a,b: (abs(a.pdgId) == abs(b.pdgId)) & (a.idx == b.idx)), axis=-1)
                          ]
        print('********2nd ', vetolep.pdgId[20:50].to_list())'''
        sel_ele = abs(sellep.pdgId) == 11
        sel_mu = abs(sellep.pdgId) == 13
        trigObj_matched_ele = ak.any(events.TrigObj[trig_obj_ele].metric_table(sellep[sel_ele]) < 0.4, axis=-1)
        trigObj_matched_mu = ak.any(events.TrigObj[trig_obj_mu].metric_table(sellep[sel_mu]) < 0.4, axis=-1)
        for msoftdrop in self.config['msoftdrop']:#['msoftdrop_GE_90', 'msoftdrop_L_90']:
            if msoftdrop == 'msoftdrop_GE_90':
                '''if self.config['isMC']:
                    selFatJet = FatJet[#(FatJet.msoftdrop > 90)
                        #&
                        (FatJet.particleNetMD_Hto4b_vs_QCDsum >0.98)
                        &
                        (ak.all(FatJet.metric_table(sellep)>0.8, axis=-1))
                    ]
                else:'''
                selFatJet = FatJet[#(FatJet.msoftdrop > 90)
                    #&
                    (FatJet.particleNetMD_Hto4b_vs_QCDsum >0.98)
                    &
                    (ak.all(FatJet.metric_table(sellep)>0.8, axis=-1))
                ]
            else:
                assert(msoftdrop == 'msoftdrop_L_90')
                selFatJet = FatJet[(FatJet.msoftdrop < 90)
                    & (ak.all(FatJet.metric_table(sellep)>0.8, axis=-1))
                ]
            delta_wq = selFatJet[:,0:1].metric_table(events.GenPart[events.GenPart.qfromW])
            if not self.config['isSignal']:
                delta_b = selFatJet[:,0:1].metric_table(events.GenPart[events.GenPart.bfromtop])
            else:
                delta_b = selFatJet[:,0:1].metric_table(events.GenPart[events.GenPart.bfroma])
            ak4jet_cleaned_wrt_lepton = ak4Jet[ak.all(ak4Jet.metric_table(sellep)>0.4, axis=-1)]
            ak4jet_cleaned = ak4Jet[
                ak.all(ak4Jet.metric_table(selFatJet[:,0:1])>0.8, axis=-1)
            ]
            ak4jet_cleaned_opshemi = ak4jet_cleaned[
                (ak4jet_cleaned.pt > 30)
                & ((ak4jet_cleaned.puId & 7) == 7)
                & ak.all(abs(ak4jet_cleaned.metric_table(selFatJet[:,0:1]#)>0.8, axis=-1)
                ,metric=lambda a, b: a.delta_phi(b))) > (3.14159265/2.), axis=-1)
            ]
            nbtag_medium = ak.num(ak4jet_cleaned[ak4jet_cleaned.btagDeepFlavB > 0.2770])
            if not self.config['isSignal']:
                bmatched_ak4jet = ak4jet_cleaned[ak.any(ak4jet_cleaned.metric_table(events.GenPart[events.GenPart.bfromtop]) <0.4, axis=-1)]
            else:
                bmatched_ak4jet = ak4jet_cleaned[ak.any(ak4jet_cleaned.metric_table(events.GenPart[events.GenPart.bfroma]) <0.4, axis=-1)]
            ak4jet_notcleaned = ak4Jet[
                ak.all(ak4Jet.metric_table(selFatJet[:,0:1])<=0.8, axis=-1)
            ]
            nbtag_notcleaned_medium = ak.num(ak4jet_notcleaned[ak4jet_notcleaned.btagDeepFlavB > 0.2770])
            ht = ak.sum(ak4jet_cleaned.pt, axis=1)
            
            selection = PackedSelection()
            selection.add("nPV", events.PV.npvsGood >= 1)
            if self.config['lepton_selection'] == 'tight_lep':
                selection.add( self.config['lepton_selection'], (ak.num(sellep_all) == 1) & (ak.num(sellep) == 1))
            else:
                selection.add( self.config['lepton_selection'], (ak.num(sellep_all) > 0)\
                    & (ak.num(sellep_tight) == 0)\
                )
            if msoftdrop == 'msoftdrop_GE_90':
                selection.add('FatJet', ak.num(selFatJet) > 0)
            else:
                selection.add(msoftdrop, (ak.num(selFatJet) > 0)\
                  & (ak.num(FatJet[FatJet.msoftdrop > 90]) == 0))
            selection.add('veto lepton', ak.num(vetolep) == 0)
            selection.add('ak4jet with nbtag_medium==0', nbtag_medium == 0)
            selection.add("met", events.MET.pt >= 20)
            selection.add("triggered", sel_triggered_1ele | sel_triggered_1mu)
            selection.add("trigObj_matched_lep",\
               (ak.any(trigObj_matched_ele, axis=-1) & (sel_triggered_1ele))\
                | (ak.any(trigObj_matched_mu, axis=-1) & sel_triggered_1mu ))
            if not self.config['isMC']:
                selection.add("dataset_check", (sel_triggered_1ele & is_triggered_1mu)==0)
            f0 = np.ones(len(events))
            if not self.config['isSignal']\
               and self.config['lepton_selection'] == 'tight_lep'\
               and msoftdrop == 'msoftdrop_L_90'\
               and len(events):
                f0 = evt_weights.transfer_factor(ht, f0)

            allcuts = selection.names#commoncuts + [lepton] + [msoftdrop]
            sel_reg = selection.all(*allcuts)
            ###################
            # FILL HISTOGRAMS
            ###################
            metp4 = met_p4(events)
            dphi_lm = sellep.metric_table(events.MET, metric=lambda a, b: a.delta_phi(b))
            dphi_lm_jet = (sellep+metp4).metric_table(selFatJet[:,0:1], metric=lambda a, b: a.delta_phi(b))
            sign_deta = sellep.metric_table(selFatJet[:,0:1],metric=lambda a, b: (a.eta - b.eta) * (1-2*(b.eta<0)))
            for syst in systList:
                weightSyst = syst
                if syst in ["central", "JERUp", "JERDown", "JESUp", "JESDown"]:
                    weightSyst = None
                if not self.config['isMC']:
                    evtWeight = f0 if self.config['lepton_selection'] == 'tight_lep'\
                                and msoftdrop == 'msoftdrop_L_90'\
                                else np.ones(len(events))
                else:
                    if self.config['lepton_selection'] == 'tight_lep'\
                       and msoftdrop == 'msoftdrop_L_90':
                        evtWeight = weights.weight(weightSyst)*f0
                    else:
                        evtWeight = weights.weight(weightSyst)
                    evtWeight_gen = weights.partial_weight(['lumiWeight', 'genWeight'])
                    
                cuts = []
                if syst == 'central':
                    for ibin, cut in enumerate(allcuts):
                        cuts.append(cut)
                        n = selection.all(*cuts).sum()
                        if msoftdrop not in output['cutflow'].keys():
                            output['cutflow'][msoftdrop] = {}
                            output['cutflow_weighted'][msoftdrop] = {}
                        output['cutflow'][msoftdrop][cut] = n
                        output['cutflow_weighted'][msoftdrop][cut] = np.sum(evtWeight[selection.all(*cuts)])
                        print('cut: ', cut, '\t', n)
                        n = np.ones(n)
                        output['hCutFlow'].fill(
                            jet=msoftdrop,
                            cutflow = (n*ibin),
                            systematic='central',
                            weight=evtWeight[selection.all(*cuts)]
                        )
                output['hLeadingFatJetPt'].fill(
                    jet=msoftdrop,
                    hLeadingFatJetPt=clip_value(
                        ak.flatten(selFatJet.pt[sel_reg][:, 0:1]),
                        min(output['hLeadingFatJetPt'].axes[2].edges),
                        max(output['hLeadingFatJetPt'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                #print(ak.sum(ak.num(ak.flatten(delta_wq[sel_reg][:,0:1],axis=0))), '\tweight ', ak.num(evtWeight[sel_reg],axis=0))
                output['dq_1'].fill(
                    jet=msoftdrop,
                    dq_1=clip_value(
                        ak.flatten(delta_wq[sel_reg],axis=None),#[:,:, 0:1], axis=None),
                        min(output['dq_1'].axes[2].edges),
                        max(output['dq_1'].axes[2].edges)
                    ),
                    systematic=syst
                    #weight=evtWeight[sel_reg]
                )
                output['db_1'].fill(
                    jet=msoftdrop,
                    db_1=clip_value(
                        ak.flatten(delta_b[sel_reg],axis=None),#[:,:, 0:1], axis=None),
                        min(output['dq_1'].axes[2].edges),
                        max(output['dq_1'].axes[2].edges)
                    ),
                    systematic=syst
                    #weight=evtWeight[sel_reg]
                )
                output['mass_vs_score'].fill(
                    jet=msoftdrop,
                    mass=clip_value(
                        ak.flatten(FatJet.particleNet_massH_Hto4b_v0[:, 0:1],axis=None),#[:,:, 0:1], axis=None),
                        min(output['mass_vs_score'].axes[2].edges),
                        max(output['mass_vs_score'].axes[2].edges)
                    ),
                    score=clip_value(
                        ak.flatten(FatJet.particleNetMD_Hto4b_vs_QCDsum[:, 0:1],axis=None),#[:,:, 0:1], axis=None),
                        min(output['mass_vs_score'].axes[3].edges),
                        max(output['mass_vs_score'].axes[3].edges)
                    ),
                    systematic=syst
                    #weight=evtWeight[sel_reg]
                )
                output['mass_50_110_score'].fill(
                    jet=msoftdrop,
                    mass_50_110_score=clip_value(
                        ak.flatten(FatJet[(FatJet.particleNet_massH_Hto4b_v0 > 50)
                                          & (FatJet.particleNet_massH_Hto4b_v0 < 110)].particleNetMD_Hto4b_vs_QCDsum
                                          ,axis=None),
                        min(output['mass_50_110_score'].axes[2].edges),
                        max(output['mass_50_110_score'].axes[2].edges)
                    ),
                    systematic=syst
                    #weight=evtWeight
                )
                output['mass_110_150_score'].fill(
                    jet=msoftdrop,
                    mass_110_150_score=clip_value(
			ak.flatten(FatJet[(FatJet.particleNet_massH_Hto4b_v0 >= 110)
                                          & (FatJet.particleNet_massH_Hto4b_v0 < 150)].particleNetMD_Hto4b_vs_QCDsum
                                          ,axis=None),
                        min(output['mass_50_110_score'].axes[2].edges),
                        max(output['mass_50_110_score'].axes[2].edges)
                    ),
                    systematic=syst
                    #weight=evtWeight
                )
                output['mass_150_200_score'].fill(
                    jet=msoftdrop,
                    mass_150_200_score=clip_value(
			ak.flatten(FatJet[(FatJet.particleNet_massH_Hto4b_v0 >= 150)
                                          & (FatJet.particleNet_massH_Hto4b_v0 < 200)].particleNetMD_Hto4b_vs_QCDsum
                                          ,axis=None),
                        min(output['mass_50_110_score'].axes[2].edges),
                        max(output['mass_50_110_score'].axes[2].edges)
                    ),
                    systematic=syst
                    #weight=evtWeight
                )
                output['dphi_lm'].fill(
                    jet=msoftdrop,
                    dphi_lm=clip_value(
                        ak.flatten(dphi_lm[sel_reg],axis=None),
                        min(output['dphi_lm'].axes[2].edges),
                        max(output['dphi_lm'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['dphi_lm_jet'].fill(
                    jet=msoftdrop,
                    dphi_lm_jet=clip_value(
                        ak.flatten(dphi_lm_jet[sel_reg],axis=None),
                        min(output['dphi_lm_jet'].axes[2].edges),
                        max(output['dphi_lm_jet'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['sign_deta'].fill(
                    jet=msoftdrop,
                    sign_deta=clip_value(
                        ak.flatten(sign_deta[sel_reg],axis=None),
                        min(output['sign_deta'].axes[2].edges),
                        max(output['sign_deta'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['wq'].fill(
                    jet=msoftdrop,
                    wq=clip_value(
                        ak.flatten(
                            abs(events.GenPart[events.GenPart.qfromW].pdgId)
                            [sel_reg],axis=None),
                        min(output['wq'].axes[2].edges),
                        max(output['wq'].axes[2].edges)
                    ),
                    systematic=syst
                    #weight=evtWeight[sel_reg]
                )
                output['hLeadingFatJetMSoftDrop'].fill(
                    jet=msoftdrop,
                    hLeadingFatJetMSoftDrop=clip_value(
                        ak.flatten(selFatJet.msoftdrop[sel_reg][:, 0:1]),
                        min(output['hLeadingFatJetMSoftDrop'].axes[2].edges),
                        max(output['hLeadingFatJetMSoftDrop'].axes[2].edges)
                        ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_vs_QCDsum'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_vs_QCDsum=clip_value(
                        ak.flatten(selFatJet.particleNetMD_Hto4b_vs_QCDsum[sel_reg][:, 0:1]),
                        min(output['FatJet_particleNetMD_Hto4b_vs_QCDsum'].axes[2].edges),
                        max(output['FatJet_particleNetMD_Hto4b_vs_QCDsum'].axes[2].edges)
                        ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNet_massH_Hto4b_v0'].fill(
                    jet=msoftdrop,
                    FatJet_particleNet_massH_Hto4b_v0=clip_value(
                        ak.flatten(selFatJet.particleNet_massH_Hto4b_v0[sel_reg][:, 0:1]),
                        min(output['hLeadingFatJetMSoftDrop'].axes[2].edges),
                        max(output['hLeadingFatJetMSoftDrop'].axes[2].edges)
                        ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['HT'].fill(
                    jet=msoftdrop,
                    HT=clip_value(
                        ht[sel_reg],
                        min(output['HT'].axes[2].edges),
                        max(output['HT'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                if self.config['isMC'] and syst == 'central' and 'LHE' in events.fields:
                    output['LHE_HT_gen'].fill(
                        jet=msoftdrop,
                        LHE_HT_gen=clip_value(
                            events.LHE.HT[sel_reg],
                            min(output['LHE_HT_gen'].axes[2].edges),
                            max(output['LHE_HT_gen'].axes[2].edges)
                        ),
                        systematic=syst,
                        weight=evtWeight_gen[sel_reg]
                    )
                    output['LHE_Njet_gen'].fill(
                        jet=msoftdrop,
                        LHE_Njet_gen=(events.LHE.Njets[sel_reg]),
                        systematic=syst,
                        weight=evtWeight_gen[sel_reg]
                    )
                    output['genweight'].fill(
                        jet=msoftdrop,
                        genWeight=(events.genWeight),
                        systematic=syst
                    )
                fatjet_sortedpnmd = selFatJet[:,0:1]
                output['hLeadingFatJetParticleNetMD_Xbb'].fill(
                    jet=msoftdrop,
                    hLeadingFatJetParticleNetMD_Xbb=(ak.flatten(fatjet_sortedpnmd.particleNetMD_Xbb[sel_reg][:, 0:1])),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                if syst == 'central':
                    cut_womet = ['nPV']#, msoftdrop]#, 'triggered']
                    ele = events.Electron#[events.Electron.mvaTTH > -0.6]
                    mu = events.Muon#[events.Muon.mvaTTH > -0.6]
                    mva_lep = ak.with_name(ak.concatenate([ele,mu], axis=1), "PtEtaPhiMCandidate")
                    mva_lep = mva_lep[ak.argsort(mva_lep.pt, ascending=False)][:,0:1]
                    output['hmet'].fill(
                        jet=msoftdrop,
                        met=np.minimum(
                            events.MET.pt[selection.all(*cut_womet)],
                            max(output['hmet'].axes[2].edges)-1
                        ),
                        systematic=syst,
                        weight=evtWeight[selection.all(*cut_womet)]
                    )
                    selection_nocut = PackedSelection()
                    selection_nocut.add("lep", ak.num(mva_lep) >= 1)
                    selection_nocut.add("fatjet", ak.num(events.FatJet) >= 1)
                    output['hmet_vs_lepmva'].fill(
                        jet=msoftdrop,
                        met=events.MET.pt[selection_nocut.all(*["lep"])],
                        lepmva=(ak.flatten(mva_lep.mvaTTH[selection_nocut.all(*["lep"])])),
                        systematic=syst,
                        weight=evtWeight[selection_nocut.all(*['lep'])]
                    )
                    output['vetolep_dr_b'].fill(
                        jet=msoftdrop,
                        vetolep_dr_b=(ak.flatten(vetolep[:,1:2].metric_table(events.GenPart[events.GenPart.bfromtop]),axis=None)),
                        systematic=syst
                    )
                    output['vetolep_btag'].fill(
                        jet=msoftdrop,
                        vetolep_btag=(ak.fill_none(ak.flatten(vetolep[:,1:2].ak4jet.btagDeepFlavB), -1)),
                        systematic=syst
                    )
                    output['lepiso'].fill(
                        jet=msoftdrop,
                        lepiso=(ak.flatten(sellep.miniPFRelIso_all[sel_reg])),
                        systematic=syst
                    )
                    output['vetolepiso'].fill(
                        jet=msoftdrop,
                        vetolepiso=(ak.flatten(vetolep.miniPFRelIso_all[:, 1:][sel_reg])),
                        systematic=syst
                    )
                    output['vetoleppt'].fill(
                        jet=msoftdrop,
                        vetoleppt=(ak.flatten(vetolep.pt[:, 1:][sel_reg])),
                        systematic=syst
                    )
                    output['hmet_vs_lepiso'].fill(
                        jet=msoftdrop,
                        met=events.MET.pt[selection_nocut.all(*["lep"])],
                        lepiso=(ak.flatten(mva_lep.miniPFRelIso_all[selection_nocut.all(*["lep"])])),
                        systematic=syst,
                        weight=evtWeight[selection_nocut.all(*['lep'])]
                    )
                    output['hsoftdrop_vs_lepiso'].fill(
                        jet=msoftdrop,
                        softdrop=(ak.flatten(events.FatJet[:,0:1].msoftdrop[selection_nocut.all(*["lep", 'fatjet'])])),
                        lepiso=(ak.flatten(mva_lep.miniPFRelIso_all[selection_nocut.all(*["lep",'fatjet'])])),
                        systematic=syst,
                        weight=evtWeight[selection_nocut.all(*['lep','fatjet'])]
                    )
                    output['hsoftdrop_vs_lepmva'].fill(
                        jet=msoftdrop,
                        softdrop=(ak.flatten(events.FatJet[:,0:1].msoftdrop[selection_nocut.all(*["lep", 'fatjet'])])),
                        lepmva=(ak.flatten(mva_lep.mvaTTH[selection_nocut.all(*["lep",'fatjet'])])),
                        systematic=syst,
                        weight=evtWeight[selection_nocut.all(*['lep','fatjet'])]
                    )
                    output['leptonsize'].fill(
                        jet=msoftdrop,
                        leptonsize=(np.concatenate((np.ones(elesize), np.ones(musize)*2))),
                        systematic=syst
                    )
                    output['ak4jetsize'].fill(
                        jet=msoftdrop,
                        ak4jetsize=(ak.num(ak4jet_cleaned)[sel_reg]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['ak4jetsize_mediumbtag'].fill(
                        jet=msoftdrop,
                        ak4jetsize_mediumbtag=(nbtag_medium[sel_reg]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['dr_lep_ak4jet'].fill(
                        jet=msoftdrop,
                        dr_lep_ak4jet=clip_value(
                            ak.flatten(ak4jet_cleaned[:,0:1].metric_table(sellep)[sel_reg]
                            ,axis=None
                            ),
                            min(output['dr_lep_ak4jet'].axes[2].edges),
                            max(output['dr_lep_ak4jet'].axes[2].edges)
                        ),
                        systematic=syst
                        #weight=evtWeight[sel_reg]
                    )
                    output['tau4_2'].fill(
                        jet=msoftdrop,
                        tau4_2=ak.flatten((fatjet_sortedpnmd.tau4/fatjet_sortedpnmd.tau2)[sel_reg][:,0:1]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['tau4_3'].fill(
                        jet=msoftdrop,
                        tau4_3=ak.flatten((fatjet_sortedpnmd.tau4/fatjet_sortedpnmd.tau3)[sel_reg][:,0:1]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['ak4_btag'].fill(
                        jet=msoftdrop,
                        ak4_btag=ak.flatten(bmatched_ak4jet.btagDeepFlavB[sel_reg]),
                        systematic=syst
                        #weight=evtWeight[sel_reg]
                    )
                output['ak8jetsize'].fill(
                    jet=msoftdrop,
                    ak8jetsize=(ak.num(selFatJet)[sel_reg]),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_Haa3b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_Haa3b=(ak.flatten(selFatJet.particleNetMD_Hto4b_Haa3b[:,0:1][sel_reg], axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_Haa4b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_Haa4b=(ak.flatten(selFatJet.particleNetMD_Hto4b_Haa4b[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_Haa34b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_Haa34b=(ak.flatten((selFatJet.particleNetMD_Hto4b_Haa3b+selFatJet.particleNetMD_Hto4b_Haa4b)[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_Haa2b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_Haa2b=(ak.flatten(selFatJet.particleNetMD_Hto4b_Haa2b[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_Haa01b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_Haa01b=(ak.flatten(selFatJet.particleNetMD_Hto4b_Haa01b[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_Haa01b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_Haa01b=(ak.flatten(selFatJet.particleNetMD_Hto4b_Haa01b[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_QCD4b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_QCD4b=(ak.flatten(selFatJet.particleNetMD_Hto4b_QCD4b[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_QCD3b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_QCD3b=(ak.flatten(selFatJet.particleNetMD_Hto4b_QCD3b[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_QCD2b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_QCD2b=(ak.flatten(selFatJet.particleNetMD_Hto4b_QCD2b[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_QCD1b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_QCD1b=(ak.flatten(selFatJet.particleNetMD_Hto4b_QCD1b[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNetMD_Hto4b_QCD0b'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto4b_QCD0b=(ak.flatten(selFatJet.particleNetMD_Hto4b_QCD0b[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                hto34b_vs_QCDsum = (selFatJet.particleNetMD_Hto4b_Haa4b+selFatJet.particleNetMD_Hto4b_Haa3b)/(selFatJet.particleNetMD_Hto4b_Haa4b+selFatJet.particleNetMD_Hto4b_Haa3b+selFatJet.particleNetMD_Hto4b_QCD0b+selFatJet.particleNetMD_Hto4b_QCD1b+selFatJet.particleNetMD_Hto4b_QCD2b+selFatJet.particleNetMD_Hto4b_QCD3b+selFatJet.particleNetMD_Hto4b_QCD4b)
                output['FatJet_particleNetMD_Hto34b_vs_QCDsum'].fill(
                    jet=msoftdrop,
                    FatJet_particleNetMD_Hto34b_vs_QCDsum=(ak.flatten(hto34b_vs_QCDsum[:,0:1][sel_reg],axis=None)),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                if self.config['isSignal']:
                    output['apt'].fill(
                        jet=msoftdrop,
                        pt=(ak.flatten(events.GenPart[events.GenPart.aboson].pt[sel_reg],axis=None)),
                        systematic=syst
                    )
                output['FatJet_particleNet_massA_Hto4b_v1'].fill(
                    jet=msoftdrop,
                    FatJet_particleNet_massA_Hto4b_v1=clip_value(
                        ak.flatten(selFatJet.particleNet_massA_Hto4b_v1[:,0:1][sel_reg],axis=None),
                        min(output['FatJet_particleNet_massA_Hto4b_v1'].axes[2].edges),
                        max(output['FatJet_particleNet_massA_Hto4b_v1'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNet_massA_Hto4b_v2'].fill(
                    jet=msoftdrop,
                    FatJet_particleNet_massA_Hto4b_v2=clip_value(
                        ak.flatten(selFatJet.particleNet_massA_Hto4b_v2[:,0:1][sel_reg],axis=None),
                        min(output['FatJet_particleNet_massA_Hto4b_v2'].axes[2].edges),
                        max(output['FatJet_particleNet_massA_Hto4b_v2'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNet_massA_Hto4b_v3'].fill(
                    jet=msoftdrop,
                    FatJet_particleNet_massA_Hto4b_v3=clip_value(
                        ak.flatten(selFatJet.particleNet_massA_Hto4b_v3[:,0:1][sel_reg],axis=None),
                        min(output['FatJet_particleNet_massA_Hto4b_v3'].axes[2].edges),
                        max(output['FatJet_particleNet_massA_Hto4b_v3'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNet_massA_Hto4b_v0'].fill(
                    jet=msoftdrop,
                    FatJet_particleNet_massA_Hto4b_v0=clip_value(
                        ak.flatten(selFatJet.particleNet_massA_Hto4b_v0[:,0:1][sel_reg],axis=None),
                        min(output['FatJet_particleNet_massA_Hto4b_v0'].axes[2].edges),
                        max(output['FatJet_particleNet_massA_Hto4b_v0'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNet_massH_Hto4b_v3'].fill(
                    jet=msoftdrop,
                    FatJet_particleNet_massH_Hto4b_v3=clip_value(
                        ak.flatten(selFatJet.particleNet_massH_Hto4b_v3[:,0:1][sel_reg],axis=None),
                        min(output['FatJet_particleNet_massH_Hto4b_v3'].axes[2].edges),
                        max(output['FatJet_particleNet_massH_Hto4b_v3'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNet_massH_Hto4b_v2'].fill(
                    jet=msoftdrop,
                    FatJet_particleNet_massH_Hto4b_v2=clip_value(
                        ak.flatten(selFatJet.particleNet_massH_Hto4b_v2[:,0:1][sel_reg],axis=None),
                        min(output['FatJet_particleNet_massH_Hto4b_v2'].axes[2].edges),
                        max(output['FatJet_particleNet_massH_Hto4b_v2'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['FatJet_particleNet_massH_Hto4b_v1'].fill(
                    jet=msoftdrop,
                    FatJet_particleNet_massH_Hto4b_v1=clip_value(
                        ak.flatten(selFatJet.particleNet_massH_Hto4b_v1[:,0:1][sel_reg],axis=None),
                        min(output['FatJet_particleNet_massH_Hto4b_v1'].axes[2].edges),
                        max(output['FatJet_particleNet_massH_Hto4b_v1'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
        return output

    def postprocess(self, accumulator):
        #pass
        return accumulator

if __name__ == '__main__':
    print("htoaa_Analysis:: main: {}".format(sys.argv)); sys.stdout.flush()

    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))

    lumiScale = 1
    sInputFiles         = config.pop("inputFiles")
    sOutputFile         = config.pop("outputFile")
    sample_category     = config['sampleCategory']
    process_name        = config.pop('process_name')
    if config['isMC']:
        luminosity          = Luminosities[config['era']][0]
        sample_crossSection = config.pop("crossSection")
        sample_nEvents      = config.pop("nEvents")
        sample_sumEvents    = config.pop("sumEvents") if config["sumEvents"] != -1 else sample_nEvents
        if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
        lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)
    config['lumiScale'] = lumiScale
    config['isSignal'] = 'SUSY' in process_name
    config['isQCD'] = 'QCD'in process_name
    sInputFiles = get_lf(sInputFiles, process_name)
    sys.stdout.flush()
    startTime = time.time()
    tracemalloc.start()
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    print(f"nEventsToAnalyze: {nEventsToAnalyze},  nEventToReadInBatch: {nEventToReadInBatch}, chunksize: {chunksize},  maxchunks: {maxchunks}")
    run = processor.Runner(
        #executor=executor,
        executor=processor.FuturesExecutor(workers=4),
        schema=schemas.NanoAODSchema,
        savemetrics=True,
        chunksize=chunksize,  #3 ** 20,  ## Governs the number of times LeptonJetProcessor "process" is called
        maxchunks=maxchunks
    )
    output, metrics = run(
        fileset={sample_category: sInputFiles},
        treename="Events",
        processor_instance=HToAATo4bProcessor(
            datasetInfo=config
        )
    )
    if 'cutflow' in output.keys():
        print("Cutflow::")
        for msoftdrop in config['msoftdrop']:
            for key in output['cutflow'][msoftdrop].keys():
                print("%10f\t%10d\t%s" % (output['cutflow_weighted'][msoftdrop][key], 
                        output['cutflow'][msoftdrop][key],
                        key)
                  )
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        sDir1 = f'evt/{config["lepton_selection"]}'# % (sample_category)
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                if not isinstance(value, hist.Hist): continue
                for jet in range(0, len(value.axes[0])):
                    jetname = value.axes[0][jet]
                    for syst in  range(0, len(value.axes[1])):
                        systname = value.axes[1][syst]
                        histname = key + '_' + systname if 'central' not in systname else key
                        if 'gen' not in histname:
                            if value.ndim == 3:
                                fOut[f'{sDir1}/{jetname}/{sample_category}/{histname}'] = value[jet,syst,:]
                            else:
                                fOut[f'{sDir1}/{jetname}/{sample_category}/{histname}'] = value[jet,syst,:,:]
                        else:
                            fOut[f'{sDir1}/{jetname}/{sample_category}/gen/{histname}'] = value[jet,syst,:]
    current_memory, peak_memory = tracemalloc.get_traced_memory() # https://medium.com/survata-engineering-blog/monitoring-memory-usage-of-a-running-python-program-49f027e3d1ba
    print(f"\n\nMemory usage:: current {current_memory / 10**6}MB;  peak {peak_memory / 10**6}MB")

    endTime = time.time()
    totalTime = endTime - startTime
    totalTime_hr  = int(totalTime/60/60)
    totalTime_min = totalTime - float(totalTime_hr * 60)
    totalTime_min = int(totalTime_min/60)
    totalTime_sec = totalTime - float(totalTime_hr * 60*60) - float(totalTime_min * 60)
    #for f in sInputFiles:
    #os.system(f'rm {f}')
    print(f"Total run time: {totalTime_hr}h {totalTime_min}m {totalTime_sec}s = {totalTime}sec ")
