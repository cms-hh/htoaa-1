from collections import OrderedDict as OD
import json
import ROOT as r
import sys

samplesList = 'Samples_2018UL.json'
with open(samplesList) as fSamplesInfo:
        samplesInfo = json.load(fSamplesInfo)

def get_xs(proc):
        for sample_info in samplesInfo:
                if samplesInfo[sample_info]['process_name'] == proc:
                        return samplesInfo[sample_info]['cross_section']
        print(f'{proc} not found')
        assert(0)

def check_sample(proc, samples_sum):
        for sample_sum in samples_sum:
                if proc in sample_sum[0]:
                        print(f'{proc} found on {sample_sum[0]}')
                        return (sample_sum[0], sample_sum[1])
        return ([], 0)

def calculate_samples_sum_event(samplesInfo):
    samples_sum = []
    for idx1, sample_info in enumerate(samplesInfo):
        #samplesInfo[sample_info]['sumEvents'] = get_histvalue(samplesInfo[sample_info]['process_name'])                                                                                                                                                     
        proc = samplesInfo[sample_info]['process_name']
        if 'QCD' in proc: continue
        proc_1_mod = proc.lower()
        sample_exist = check_sample(proc, samples_sum)
        if len(sample_exist[0]): continue
        sample_sum = []
        gensum = 0
        for idx2, sample_info_2 in enumerate(samplesInfo):
                if not (idx2 > idx1): continue
                add_sample = False
                proc_2 = samplesInfo[sample_info_2]['process_name']
                proc_2_mod = proc_2.lower()
                if proc_1_mod == proc_2_mod:
                        add_sample = True
                else:
                        proc_2_mod = proc_2_mod.split('_psweights')[0]
                        proc_2_mod = proc_2_mod.split('_ext')[0]
                        proce_1_mod = proc_1_mod.split('_psweights')[0]
                        proce_1_mod = proc_1_mod.split('_ext')[0]
                        if proc_2_mod == proc_1_mod:
                                add_sample = True
                if add_sample:
                        assert(samplesInfo[sample_info_2]['cross_section'] == samplesInfo[sample_info]['cross_section'])
                        if len(sample_sum) == 0:
                                sample_sum.append(proc)
                                gensum += samplesInfo[sample_info]['sumEvents']
                        sample_sum.append(proc_2)
                        gensum += samplesInfo[sample_info_2]['sumEvents']
        if len(sample_sum): samples_sum.append([sample_sum, gensum])
    return samples_sum

samples_sum_events = calculate_samples_sum_event(samplesInfo)

pt = OD([
    ('HT-0To70', [0, 70]),
    ('HT-70To100', [70,100]),
    ('HT-100To200', [100,200]),
    ('HT-200To400', [200,400]),
    ('HT-400To600', [400,600]),
    ('HT-600To800', [600,800]),
    ('HT-800To1200', [800, 1200]),
    ('HT-1200To2500', [1200,2500]),
    ('HT-2500ToInf', [2500]),
])
njet  = OD([
        ('0Jet', 0),
        ('1Jet', 1),
        ('2Jet', 2),
        ('3Jet', 3),
        ('4Jet', 4)
])

f = r.TFile('hadd_stitch.root')
nevents = OD()
for samplename in [d.GetName() for d in f.Get('evt').GetListOfKeys()]:
    hist = f.Get(f'evt/{samplename}/hCutFlow_central')
    assert(samplename not in nevents.keys())
    nevents[samplename] = {}
    for pt_idx, pt_key in enumerate(pt.keys()):
            nevents[samplename][pt_key] = {}
            for njet_idx, njet_key in enumerate(njet.keys()):
                    nevents[samplename][pt_key][njet_key] = hist.GetBinContent(pt_idx+1, njet_idx+1)
    sample_exist = check_sample(samplename, samples_sum_events)
    nevents[samplename]['nevents'] = hist.Integral()
    nevents[samplename]['xs'] = get_xs(samplename)
    #if '70To100' in samplename: print('nevents: ', samplename, '\t', hist.Integral())
    if len(sample_exist[0]):
            for sample in sample_exist[0]:
                    if samplename == sample: continue
                    hist = f.Get(f'evt/{sample}/hCutFlow_central')
                    for pt_idx, pt_key in enumerate(pt.keys()):
                            for njet_idx, njet_key in enumerate(njet.keys()):
                                    nevents[samplename][pt_key][njet_key] += hist.GetBinContent(pt_idx+1, njet_idx+1)
                    #if '70To100' in samplename: print('nevents: ', samplename)
                    nevents[samplename]['nevents'] += hist.Integral()

with open('stitching.json', "w") as fSampleInfo:
        json.dump(nevents, fSampleInfo, indent=4)

#sys.exit()
stitchinginfo = OD()
inclu_xs = 0
exclu_xs = 0
for htbin in pt.keys():
        #if htbin != 'HT-200To400': continue
        assert(htbin not in stitchinginfo.keys())
        stitchinginfo[htbin] = {}
        for nj in njet.keys():
                assert(nj not in stitchinginfo[htbin].keys())
                stitchinginfo[htbin][nj] = {}
                xs = 0
                nevent = 0
                inclusive = 0
                ht_exclusive = 0
                nj_exclusive = 0
                for sampleinfo in nevents.keys():
                        if ht_exclusive and htbin in sampleinfo:
                                continue
                        if nj_exclusive and nj in sampleinfo:
                                continue
                        if nevents[sampleinfo][htbin][nj] == 0:
                                continue
                        if htbin in sampleinfo: ht_exclusive = True
                        if nj in sampleinfo: nj_exclusive = True
                        if htbin not in sampleinfo and nj not in sampleinfo: inclusive = True
                        xs += nevents[sampleinfo]["xs"] * nevents[sampleinfo][htbin][nj] / nevents[sampleinfo]["nevents"]
                        #print('xs: ', sampleinfo, '\t', nevents[sampleinfo]["xs"], '\t', nevents[sampleinfo][htbin] , '\t', nevents[sampleinfo]["nevents"])
                        nevent += nevents[sampleinfo][htbin][nj]
                if ht_exclusive and nj_exclusive and inclusive:
                        xs /= 3.
                elif ht_exclusive and inclusive and not nj_exclusive:
                        xs *= 0.5
                elif nj_exclusive and inclusive and not ht_exclusive:
                        xs *= 0.5
                stitchinginfo[htbin][nj]['xs'] = xs
                stitchinginfo[htbin][nj]['nevent'] = nevent
                if ht_exclusive or nj_exclusive or inclusive:
                        exclu_xs += xs
print(exclu_xs, '\t', inclu_xs)
with open('stitchinginfo.json', "w") as fSampleInfo:
        json.dump(stitchinginfo, fSampleInfo, indent=4)
