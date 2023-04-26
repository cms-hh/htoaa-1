import awkward as ak
import sys

class ObjectSelection:
    def __init__(self, era):
        self.era = era

        self.tagger_btagDeepB = 'DeepCSV'
        self.wp_btagDeepB = 'M'

        self.FatJetPtThsh  = 170
        self.FatJetEtaThsh = 2.4

        self.FatJetMSoftDropThshLow  = 90
        self.FatJetMSoftDropThshHigh = 200

        self.FatJetParticleNetMD_Xbb_Thsh = 0.8

        self.JetPtThshForHT = 30.0
        self.JetEtaThshForHT = 2.4

        self.nFatJetMin = 1
        self.GenHTThsh  = 100.0
        self.LHEHTThsh  = 100.0

        self.leadingelePtThsh = 26
        self.leadingmuPtThsh = 26
        self.leadingeleEtaThsh = 2.3
        self.leadingmuEtaThsh = 2.4


    def selectGenHiggs(self, events):
        maskGenHiggs = (
            (events.GenPart.pdgId  == 25) & # pdgId:: 25: H0
            (events.GenPart.status == 62)   # statu 62: outgoing subprocess particle with primordial kT included https://pythia.org/latest-manual/ParticleProperties.html           
        )
        return events.GenPart[maskGenHiggs]

    def selectGenABoson(self, events):
        maskGenA = (
            (events.GenPart.pdgId == 36)
        )
        return events.GenPart[maskGenA]

    def GenHT(self, events):
        maskForGenHT = (
            (events.GenJet.pt > self.JetPtThshForHT) &
            (abs(events.GenJet.eta) < self.JetEtaThshForHT)
        )
        selGenJetPt = events.GenJet[maskForGenHT].pt
        GenHT = ak.sum(selGenJetPt, axis=-1)
        return GenHT

    def selectFatJets(self, events):
        mask = ((events.FatJet.pt >= self.FatJetPtThsh)
                & (abs(events.FatJet.eta) <= self.FatJetEtaThsh)
                & (events.FatJet.msoftdrop >= self.FatJetMSoftDropThshLow)
                )
        return events.FatJet[mask]

    def selectMuons(self, events):
        mask = ((events.Muon.pt >5) &
                (abs(events.Muon.eta) < 2.4) &
                (events.Muon.dxy < 0.05) &
                (events.Muon.dz < 0.1 ) &
                (events.Muon.miniPFRelIso_all < 0.4) &
                (events.Muon.sip3d <8 ) &
                (events.Muon.mvaTTH >= 0.5)
        )

        return events.Muon[mask]

    def selectElectrons(self, events):
        mask = ((events.Electron.pt >7) &
                (abs(events.Electron.eta) < 2.3) &
                (events.Electron.dxy < 0.05) &
                (events.Electron.dz < 0.1 ) &
                (events.Electron.miniPFRelIso_all < 0.4) &
                (events.Electron.sip3d <8 ) &
                (events.Electron.mvaTTH >= 0.3)
        )

        return events.Electron[mask]
