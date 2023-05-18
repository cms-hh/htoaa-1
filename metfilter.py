import awkward as ak
from functools import reduce
from operator import and_, or_

#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
metfilters = {}
metfilters['goodVertices'] = {'MC': True, 'Data': True}
metfilters['globalSuperTightHalo2016Filter'] = {'MC': True, 'Data': True}
metfilters['HBHENoiseFilter'] = {'MC': True, 'Data': True}
metfilters['HBHENoiseIsoFilter'] = {'MC': True, 'Data': True}
metfilters['EcalDeadCellTriggerPrimitiveFilter'] = {'MC': True, 'Data': True}
metfilters['BadPFMuonFilter'] = {'MC': True, 'Data': True}
metfilters['eeBadScFilter'] = {'MC': False, 'Data': True}

def get_metfilters(isMC):
    mc_filters = []
    data_filters = []
    for f, info in metfilters.items():
        if info['MC']:
            mc_filters.append(f)
        else:
            data_filters.append(f)
    return mc_filters if isMC else data_filters

def apply_metfilters(isMC, events):
    mfs = get_metfilters(isMC)
    mask = reduce(and_, (getattr(events.Flag, str(mf)) for mf in mfs))
    return events[mask]
