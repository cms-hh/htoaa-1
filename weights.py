from ht_njets_binning import htbin, njetbin
from auxiliary import top_pT_reweighting
import json
import numpy as np
import awkward as ak
import uproot
from coffea.analysis_tools import Weights

class event_weights():
    def __init__(self, events, config):
        self.events = events
        self.nevents = len(events)
        self.config = config
        self.pu_hardmax = 3
        self.pu_maxshift = 0.0025
        stitchingfile = 'stitchinginfo.json'
        with open(stitchingfile) as fSamplesInfo:
            self.stitchinginfo = json.load(fSamplesInfo)

    def calculate_genweight(self):
        genweight_unique = set(self.events.genWeight)
        count = {}
        for gu in genweight_unique:
            count[gu] = ak.sum(ak.where(self.events.genWeight==gu, 1,0))
        count = sorted([(k,v) for k,v in count.items()], key=lambda kv: kv[1], reverse=True)
        assert(count[0][0] >= 0)
        wsign = np.copysign(np.ones(self.nevents), self.events.genWeight)
        self.events.genWeight = ak.where(
            abs(self.events.genWeight) > 3*count[0][0],\
            3*count[0][0],
            abs(self.events.genWeight)
        )
        self.events.genWeight = self.events.genWeight*wsign

    def calculate_stitchingweight(self):

        for idx_pt, (k_pt,v_pt) in enumerate(htbin.items()):
            pt_min = v_pt[0]
            if len(v_pt) > 1:
                pt_max = v_pt[1]
            else:
                pt_max = 10000000000
            for idx_njet, (k_njet,v_njet) in enumerate(njetbin.items()):
                njet_min = v_njet
                self.stitch_weight = ak.where(
                    (self.events.LHE.HT >= pt_min)
                    & (self.events.LHE.HT < pt_max)
                    & (self.events.LHE.Njets == njet_min),
                    59.83*1000*self.stitchinginfo[k_pt][k_njet]['xs']/self.stitchinginfo[k_pt][k_njet]['nevent']\
                    if self.stitchinginfo[k_pt][k_njet]['nevent'] else self.stitch_weight,
                    self.stitch_weight
                )

    def calculate_puweight(self):
        #https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#Recommended_cross_section
        pu_data = uproot.open('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root')['pileup'].to_hist()
        pu_mc = uproot.open(f'/afs/cern.ch/work/s/snandan/public/myforkhtoaa/htoaa-1/countpileup/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12/output/analyze_htoaa_SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_0_0.root')[f'evt/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12']
        pu_w = np.divide(pu_data.values()/pu_data.values().sum(), pu_mc.values()/pu_mc.values().sum(), where=pu_mc.values()>0, out=np.ones_like(pu_data.values()))
        clipped = pu_w.copy()
        while True:
            clipped= np.minimum(clipped, self.pu_hardmax)
            sf = 1/(1+(np.sum(clipped)-np.sum(pu_w))/np.sum(pu_w))
            clipped *= sf
            if sf <= (1+self.pu_maxshift):
                break
        clipped_pu_w = clipped
        self.pileupWeight = np.ones(self.nevents)
        for npu in range(0, len(pu_data.values()+1)):
            self.pileupWeight = ak.where(
                np.logical_and(
                    np.greater_equal(self.events.Pileup.nTrueInt, npu),
                    np.less(self.events.Pileup.nTrueInt, npu+1)
                ),
                clipped_pu_w[npu],
                self.pileupWeight
            )

    def calculate_weight(self):
        weights = Weights(self.nevents, storeIndividual=True)
        self.stitch_weight = np.full(self.nevents, self.config["lumiScale"])
        if self.config['isMC']:
            self.calculate_genweight()
            self.calculate_puweight()
            if self.config["applystitching"]:
                self.calculate_stitchingweight()
                weights.add(
                    "lumiWeight",
                    weight=self.stitch_weight
                )
            else:
                weights.add(
                    "lumiWeight",
                    weight=np.full(self.nevents, self.config["lumiScale"])
                )
            weights.add(
                "genWeight",
                weight=self.events.genWeight
            )
            weights.add(
                "puWeight",
                weight=self.pileupWeight
            )
            if self.config['sampleCategory'].startswith('TT'):
                topt_reweiging = ak.to_numpy(top_pT_reweighting(self.events.GenPart)).squeeze()
                weights.add(
                    "topt_reweiging",
                    weight=topt_reweiging,
                    weightUp=topt_reweiging**2,  # w**2
                    weightDown=ak.ones_like(topt_reweiging)  # w = 1
                )
            
        return weights

    def transfer_factor(self, ht, f0):
    
        '''prev_v = ''
        for idx, (v,w) in enumerate(tw.keys()):
        if idx == 0:
        prev_v = v
        f0 = ak.where((ht < v), w , f0)
        elif idx == len(tw.keys())-2:
        f0 = ak.where((ht >= v), tf[], f0)'''
        with open('transfer_factor.json') as f:
            tf = json.load(f)
    
        f0 = ak.where((ht < 150), tf['150'] , f0)
        f0 = ak.where((ht >= 150) & (ht < 200), tf['200'], f0)
        f0 = ak.where((ht >= 200) & (ht < 225), tf['225'], f0)
        f0 = ak.where((ht >= 225) & (ht < 250), tf['250'], f0)
        f0 = ak.where((ht >= 250) & (ht < 275), tf['275'], f0)
        f0 = ak.where((ht >= 275) & (ht < 300), tf['300'], f0)
        f0 = ak.where((ht >= 300) & (ht < 325), tf['325'], f0)
        f0 = ak.where((ht >= 325) & (ht < 350), tf['350'], f0)
        f0 = ak.where((ht >= 350) & (ht < 375), tf['375'], f0)
        f0 = ak.where((ht >= 375) & (ht < 400), tf['400'], f0)
        f0 = ak.where((ht >= 400) & (ht < 425), tf['425'], f0)
        f0 = ak.where((ht >= 425) & (ht < 450), tf['450'], f0)
        f0 = ak.where((ht >= 450) & (ht < 475), tf['475'], f0)
        f0 = ak.where((ht >= 475) & (ht < 500), tf['500'], f0)
        f0 = ak.where((ht >= 500) & (ht < 550), tf['550'], f0)
        f0 = ak.where((ht >= 550) & (ht < 600), tf['600'], f0)
        f0 = ak.where((ht >= 600), tf['1000'], f0)
        
        return f0
