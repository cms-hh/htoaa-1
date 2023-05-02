import awkward as ak
from coffea.nanoevents.methods import vector

def genparticle(events, part, idx):
    v= ak.zip(
        {
            "pt"  : events.GenPart[part][:, idx].pt,
            "eta" : events.GenPart[part][:, idx].eta,
            "phi" : events.GenPart[part][:, idx].phi,
            "mass": events.GenPart[part][:, idx].mass,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior,
    )
    return v