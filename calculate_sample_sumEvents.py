def check_sample(proc, samples_sum):
        for sample_sum in samples_sum:
                if proc in sample_sum[0]:
                        print(f'{proc} found on {sample_sum[0]}')
                        return (True, sample_sum[1])
        return (False, 0)

def calculate_samples_sum_event(samplesInfo):
    samples_sum = []
    for idx1, sample_info in enumerate(samplesInfo):
        #samplesInfo[sample_info]['sumEvents'] = get_histvalue(samplesInfo[sample_info]['process_name'])
        proc = samplesInfo[sample_info]['process_name']
        if 'QCD' in proc: continue
        proc_1_mod = proc.lower()
        sample_exist = check_sample(proc, samples_sum)
        if sample_exist[0]: continue
        sample_sum = []
        gensum = 0
        for idx2, sample_info_2 in enumerate(samplesInfo):
                if not (idx2 > idx1): continue
                add_sample = False
                proc_2 = samplesInfo[sample_info_2]['process_name']
                if 'QCD' in proc_2: continue
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
                        if len(sample_sum) == 0:
                                sample_sum.append(proc)
                                gensum += samplesInfo[sample_info]['sumEvents']
                        sample_sum.append(proc_2)
                        gensum += samplesInfo[sample_info_2]['sumEvents']
        if len(sample_sum): samples_sum.append([sample_sum, gensum])
    return samples_sum
