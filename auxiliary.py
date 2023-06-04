import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector

def met_p4(events):
    v= ak.zip(
        {
            "pt"  : events.MET.pt,
            "eta" : np.zeros(len(events)),
            "phi" : events.MET.phi,
            "mass": np.zeros(len(events)),
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior,
    )
    return v

def top_pT_sf_formula(pt):
    return np.exp(
        -2.02274e-01 + 1.09734e-04 * pt + -1.30088e-07 * pt ** 2 + (5.83494e01 / (pt + 1.96252e02))
    )

def top_pT_reweighting(gen):
    """
    Apply this SF only to TTbar datasets!
    Documentation:
        - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting 
        - https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf
        - https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf
    """
    top = gen[(gen.pdgId == 6) & gen.hasFlags(["isLastCopy"])]
    anti_top = gen[(gen.pdgId == -6) & gen.hasFlags(["isLastCopy"])]
    return np.sqrt(top_pT_sf_formula(top.pt) * top_pT_sf_formula(anti_top.pt))
