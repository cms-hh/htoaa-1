import os
import sys
import json
import glob
import argparse
from collections import OrderedDict as OD
import numpy as np
import time
from datetime import datetime
import copy
import re
from condor_helper import *
from calculate_sample_sumEvents import *

from htoaa_Settings import *
from htoaa_Samples import (
    Samples2018,
    kData
)

sAnalysis         = "whtoaa_Analysis_wCoffea.py"  # "htoaa_Analysis.py"
sConfig           = "config_htoaa.json"
sRunCommandFile   = "1_RunCommand.txt"
sJobSubLogFile    = "1_JobSubmission.log"
sOpRootFile       = "analyze_htoaa_$SAMPLE_$STAGE_$IJOB.root"

printLevel = 3

if __name__ == '__main__':

    print("htoaa_Wrapper:: main: {}".format(sys.argv)); sys.stdout.flush()
    
    parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
    parser.add_argument('-era', dest='era',   type=str, default=Era_2018, choices=[Era_2016, Era_2017, Era_2018], required=False)
    parser.add_argument('-run_mode',          type=str, default='local', choices=['local', 'condor'])
    parser.add_argument('-v', '--version',    type=str, default=None, required=True)
    parser.add_argument('-samples',           type=str, default=None, help='samples to run seperated by comma')
    parser.add_argument('-nFilesPerJob',      type=int, default=10)
    parser.add_argument('-nResubMax',         type=int, default=80)
    parser.add_argument('-ResubWaitingTime',  type=int, default=15, help='Resubmit failed jobs after every xx minutes')
    parser.add_argument('-iJobSubmission',    type=int, default=0,  help='Job submission iteration. Specify previous last job submittion iteration if script terminated for some reason.')
    parser.add_argument('-dryRun',            action='store_true', default=False)
    parser.add_argument('-rf', '--run_file',    type=str, default='analysis', choices=['analysis', 'count_genweight', 'stitch'])
    parser.add_argument('-as', '--applystitching', action='store_false', default=True)

    args=parser.parse_args()
    print("args: {}".format(args))

    era              = args.era
    run_mode         = args.run_mode
    nFilesPerJob     = args.nFilesPerJob
    selSamplesToRun  = args.samples
    anaVersion       = args.version
    nResubmissionMax = args.nResubMax
    ResubWaitingTime = args.ResubWaitingTime
    iJobSubmission   = args.iJobSubmission
    dryRun           = args.dryRun
    run_file         = args.run_file
    applystitching   = args.applystitching

    if run_file == 'count_genweight':
            sAnalysis = 'countSumEventsInSample.py'
    elif run_file == 'stitch':
        sAnalysis = 'count_LHE_HT.py'

    pwd = os.getcwd()
    DestinationDir = "./%s" % (anaVersion)
    if "/afs/cern.ch/" in pwd: DestinationDir = "/afs/cern.ch/work/s/snandan/public/myforkhtoaa/htoaa-1/%s" % (anaVersion)
    sFileRunCommand = "%s/%s" % (DestinationDir, sRunCommandFile)
    sFileJobSubLog  = "%s/%s" % (DestinationDir, sJobSubLogFile)
    if not os.path.exists(DestinationDir): os.mkdir(DestinationDir)

    # save run command into a .txt tile
    with open(sFileRunCommand, 'a') as fRunCommand:
        datatime_now = datetime.now()
        sCommand = ' '.join(sys.argv)
        fRunCommand.write('%s %s \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), sCommand))
    do_hadd = False
    samplesList = None
    samplesInfo = None
    Luminosity  = None
    if era == Era_2018:
        samplesList = 'Samples_2018UL.json'
        with open(samplesList) as fSamplesInfo:
            samplesInfo = json.load(fSamplesInfo)
        Luminosity = Luminosities[era][0]

    selSamplesToRun_list = []
    if selSamplesToRun:
        selSamplesToRun_list = selSamplesToRun.split(',')

    condor_ana = condor(sAnalysis, sFileJobSubLog)
    while iJobSubmission <= nResubmissionMax:
        if dryRun and iJobSubmission >0: break
        condor_ana.initialize_list()
        print('\n\n%s \t Startiing iJobSubmission for analysis: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
        for leptonselection in ['Fake', 'Tight']:
            for met in ['met_GE20', 'met_L20']:
                #condor_ana.initialize_list()
                samples_sum_events = calculate_samples_sum_event(samplesInfo)
                for dbsname, sampleInfo in samplesInfo.items():
                    sample_category = sampleInfo['sample_category']
                    process_name = sampleInfo['process_name']
                    if len(selSamplesToRun_list) > 0:
                        skipThisSample = True
                        for selSample in selSamplesToRun_list:
                            if not applystitching:
                                if sample_category == 'WJets' and (re.findall('^W[1-4]Jets', process_name) or ('HT-' in process_name)):
                                    break
                            if leptonselection == 'Fake' and 'SUSY' in process_name: break
                            elif selSample in process_name or selSample in sample_category:
                                skipThisSample= False
                                break
                    if skipThisSample:
                        continue

                    fileList = sampleInfo[sampleFormat]
                    files = []
                    for iEntry in fileList:
                        if "*" in iEntry:  files.extend( glob.glob( iEntry ) )
                        else:  files.append( iEntry )
                    sample_cossSection = sampleInfo["cross_section"] if (sample_category != kData) else None
                    sample_nEvents     = sampleInfo["nEvents"]
                    sample_exist = check_sample(process_name, samples_sum_events)
                    if sample_exist[0]:
                        sample_sumEvents   = sample_exist[1]
                    else:
                        sample_sumEvents   = sampleInfo["sumEvents"] if (sample_category != kData) else None

                    nSplits = int(len(files) / nFilesPerJob) if nFilesPerJob > 0 else 1
                    nSplits = 1 if nSplits==0 else nSplits
                    files_splitted = np.array_split(files, nSplits)
                    sampledir = os.path.join(DestinationDir, leptonselection, met, process_name)
                    for iJob in range(len(files_splitted)):
                        #config = config_Template.deepcopy()
                        config = copy.deepcopy(config_Template)
                        # Job related files
                        outputdir = os.path.join(sampledir, 'output')
                        if not os.path.exists(outputdir): os.makedirs(outputdir)
                        sOpRootFile_to_use    = '%s/%s' % (outputdir, sOpRootFile)
                        sOpRootFile_to_use    = sOpRootFile_to_use.replace('$SAMPLE', process_name)
                        sOpRootFile_to_use    = sOpRootFile_to_use.replace('$STAGE', str(0))
                        sOpRootFile_to_use    = sOpRootFile_to_use.replace('$IJOB', str(iJob))
                        configdir = os.path.join(sampledir, 'config')
                        if not os.path.exists(configdir): os.makedirs(configdir)
                        sConfig_to_use = sOpRootFile_to_use.replace('.root', '_config.json').replace(outputdir, configdir)

                        condordir = os.path.join(sampledir, 'condor')
                        if not os.path.exists(condordir): os.makedirs(condordir)
                        sCondorExec_to_use    = sOpRootFile_to_use.replace('.root', '_condor_exec.sh').replace(outputdir, condordir)
                        sCondorSubmit_to_use  = sOpRootFile_to_use.replace('.root', '_condor_submit.sh').replace(outputdir, condordir)
                        sCondorLog_to_use     = sOpRootFile_to_use.replace('.root', '_condor.log').replace(outputdir, condordir)
                        sCondorOutput_to_use  = sOpRootFile_to_use.replace('.root', '_condor.out').replace(outputdir, condordir)
                        sCondorError_to_use   = sOpRootFile_to_use.replace('.root', '_condor.error').replace(outputdir, condordir)
                        condor_ana.initialize_files(sOpRootFile_to_use, configdir,
                                                    sConfig_to_use, condordir, sCondorExec_to_use,
                                                    sCondorSubmit_to_use, sCondorLog_to_use, sCondorOutput_to_use,
                                                    sCondorError_to_use
                                                )
                        jobStatusForJobSubmission  = [0, 3, 4, 5]
                        # Check if job related file exist or not
                        jobStatus = condor_ana.check_jobstatus(printLevel >= 3)
                        #if iJobSubmission == 0:
                        if jobStatus == 0:
                            config["era"] = era
                            config["inputFiles"] = list( files_splitted[iJob] )
                            config["outputFile"] = condor_ana.sOpRootFile_to_use
                            config["sampleCategory"] = sample_category
                            config['process_name'] = process_name
                            config["isMC"] = 1 if (sample_category != kData) else 0
                            config["leptonselection"] = leptonselection
                            config["met"] = met
                            config["nEvents"] = sample_nEvents
                            if (sample_category != kData):
                                config["crossSection"] = sample_cossSection
                                config["sumEvents"] = sample_sumEvents
                                if applystitching:
                                    if sample_category == 'WJets' and 'WJetsToQQ' not in config['process_name']:
                                        config['applystitching'] = True
                            else:
                                del config["crossSection"]
                                del config["sumEvents"]
                            if 'EGamma' in config['process_name']:
                                config['use_triggers_1mu'] = False
                                config['use_triggers_jet'] = False
                            elif 'SingleMuon' in config['process_name']:
                                config['use_triggers_1e'] = False
                                config['use_triggers_jet'] = False

                            with open(condor_ana.sConfig_to_use, "w") as fConfig:
                                json.dump( config,  fConfig, indent=4)

                            condor_ana.writeCondorExecFile()

                        if jobStatus in jobStatusForJobSubmission: #[0, 3, 4]:
                            if jobStatus == [3, 5]:
                                # save previos .out and .error files with another names
                                sCondorOutput_vPrevious = condor_ana.sCondorOutput_to_use.replace('.out', '_v%d.out' % (iJobSubmission-1))
                                sCondorError_vPrevious  = condor_ana.sCondorError_to_use.replace('.error', '_v%d.error' % (iJobSubmission-1))
                                os.rename(condor_ana.sCondorOutput_to_use, sCondorOutput_vPrevious)
                                os.rename(condor_ana.sCondorError_to_use,  sCondorError_vPrevious)

                            increaseJobFlavour = False
                            if jobStatus == 4:
                                increaseJobFlavour = True
                            condor_ana.writeCondorSumitFile(increaseJobFlavour)

                        if jobStatus in [1, 2]:
                            # job is either running or succeeded
                            continue

                        if run_mode == 'condor':
                            cmd1 = "condor_submit %s" % condor_ana.sCondorSubmit_to_use
                            print("Now:  %s " % cmd1)
                            if not dryRun:
                                os.system(cmd1)
                        else:
                            pass
        # write JobSubmission status report in JobSubLog file
        condor_ana.update_JobSubLog(iJobSubmission)
        if dryRun:
            print('%s \t druRun with iJobSubmission: %d  \nTerminating...\n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
            #continue#exit(0)
        elif len(condor_ana.OpRootFiles_Target) == len(condor_ana.OpRootFiles_Exist):
            break
        else:
            do_hadd = True
            time.sleep( ResubWaitingTime * 60 )
        iJobSubmission += 1

    with open(condor_ana.fJobSubLog, 'a') as fJobSubLog:
        fJobSubLog.write('%s \t Jobs are done. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
    print('%s \t Jobs are done. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
del condor_ana
if dryRun: sys.exit()
if 1:
    condor_hadd = condor('hadd.py', sJobSubLogFile)
    while iJobSubmission <= nResubmissionMax:
        print('\n\n%s \t Startiing iJobSubmission for hadd: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
        iJobSubmission = 0
        condor_hadd.initialize_list()
        for leptonselection in ['Fake', 'Tight']:
            for met in ['met_GE20', 'met_L20']:
                inputdirs = glob.glob(f'{DestinationDir}/{leptonselection}/{met}/*/output/')
                for inputdir in inputdirs:
                    process_name = inputdir.split('/')[-3]
                    if 'SUSY' in process_name and leptonselection == 'Fake': continue
                    if process_name == 'hadd': continue
                    sOpRootFile_to_use    = os.path.join(inputdir, f'hadd_{process_name}.root')
                    sConfig_to_use        = sOpRootFile_to_use.replace('.root', '_config.json').replace('output', 'config')
                    sCondorExec_to_use    = sOpRootFile_to_use.replace('.root', '_condor_exec.sh').replace('output', 'condor')
                    sCondorSubmit_to_use  = sOpRootFile_to_use.replace('.root', '_condor_submit.sh').replace('output', 'condor')
                    sCondorLog_to_use     = sOpRootFile_to_use.replace('.root', '_condor.log').replace('output', 'condor')
                    sCondorOutput_to_use  = sOpRootFile_to_use.replace('.root', '_condor.out').replace('output', 'condor')
                    sCondorError_to_use   = sOpRootFile_to_use.replace('.root', '_condor.error').replace('output', 'condor')
                    condor_hadd.initialize_files(sOpRootFile_to_use, configdir,
                                                 sConfig_to_use, condordir, sCondorExec_to_use,
                                                 sCondorSubmit_to_use, sCondorLog_to_use, sCondorOutput_to_use,
                                                 sCondorError_to_use
                                             )
                    jobStatus = condor_hadd.check_jobstatus(printLevel >= 3)
                    config = {}
                    if jobStatus == 0:
                        config["inputdir"] = inputdir
                        config["outputFile"] = sOpRootFile_to_use
                        with open(sConfig_to_use, "w") as fConfig:
                            json.dump( config,  fConfig, indent=4)
                        condor_hadd.writeCondorExecFile(f'-p {config["inputdir"]} -o {config["outputFile"]}')
                    if jobStatus in jobStatusForJobSubmission: #[0, 3, 4]:
                        if jobStatus == [3, 5]:
                            # save previos .out and .error files with another names
                            sCondorOutput_vPrevious = sCondorOutput_to_use.replace('.out', '_v%d.out' % (iJobSubmission-1))
                            sCondorError_vPrevious  = sCondorError_to_use.replace('.error', '_v%d.error' % (iJobSubmission-1))
                            os.rename(condor_hadd.sCondorOutput_to_use, sCondorOutput_vPrevious)
                            os.rename(condor_hadd.sCondorError_to_use,  sCondorError_vPrevious)
                        increaseJobFlavour = False
                        if jobStatus == 4:
                            increaseJobFlavour = True
                        condor_hadd.writeCondorSumitFile(increaseJobFlavour)
                        if jobStatus in [1, 2]:
                            # job is either running or succeeded
                            continue
                        if run_mode == 'condor':
                            cmd1 = "condor_submit %s" % condor_hadd.sCondorSubmit_to_use
                            print("Now:  %s " % cmd1)
                            if not dryRun:
                                os.system(cmd1)
                        else:
                            pass
        if len(condor_hadd.OpRootFiles_Target) == len(condor_hadd.OpRootFiles_Exist):
            break
        else:
            time.sleep( ResubWaitingTime * 60 )
            iJobSubmission += 1
    with open(condor_hadd.fJobSubLog, 'a') as fJobSubLog:
        fJobSubLog.write('%s \t Jobs are done for hadd. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
    print('%s \t Jobs are done for hadd. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))

iJobSubmission = 0
while iJobSubmission <= nResubmissionMax:
    print('\n\n%s \t Startiing iJobSubmission for hadd_stage1: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
    for leptonselection in ['Fake', 'Tight']:
        for met in ['met_GE20', 'met_L20']:
            condor_hadd.initialize_list()
            inputdir = os.path.join(DestinationDir, leptonselection, met)
            outputdir = os.path.join(inputdir, 'hadd', 'output')
            if not os.path.exists(outputdir): os.makedirs(outputdir)
            configdir = outputdir.replace('output', 'config')
            if not os.path.exists(configdir): os.makedirs(configdir)
            sOpRootFile_to_use    = os.path.join(outputdir, f'hadd_stage1.root')
            sConfig_to_use        = sOpRootFile_to_use.replace('.root', '_config.json').replace('output', 'config')
            condordir = outputdir.replace('output','condor')
            if not os.path.exists(condordir): os.makedirs(condordir)
            sCondorExec_to_use    = sOpRootFile_to_use.replace('.root', '_condor_exec.sh').replace('output', 'condor')
            sCondorSubmit_to_use  = sOpRootFile_to_use.replace('.root', '_condor_submit.sh').replace('output', 'condor')
            sCondorLog_to_use     = sOpRootFile_to_use.replace('.root', '_condor.log').replace('output', 'condor')
            sCondorOutput_to_use  = sOpRootFile_to_use.replace('.root', '_condor.out').replace('output', 'condor')
            sCondorError_to_use   = sOpRootFile_to_use.replace('.root', '_condor.error').replace('output', 'condor')
            condor_hadd.initialize_files(sOpRootFile_to_use, configdir,
                                         sConfig_to_use, condordir, sCondorExec_to_use,
                                         sCondorSubmit_to_use, sCondorLog_to_use, sCondorOutput_to_use,
                                         sCondorError_to_use
                                     )
            jobStatus = condor_hadd.check_jobstatus(printLevel >= 3)
            config = {}
            if jobStatus == 1: continue
            if jobStatus == 0:
                config["inputdir"] = inputdir
                config["outputFile"] = sOpRootFile_to_use
                config['wildcard'] = '*/output/hadd*.root'
                with open(condor_hadd.sConfig_to_use, "w") as fConfig:
                    json.dump( config,  fConfig, indent=4)
                condor_hadd.writeCondorExecFile(f'-p {config["inputdir"]} -o {config["outputFile"]} -w {config["wildcard"]}')
            if jobStatus in [0, 3, 4, 5]: #[0, 3, 4]:
                if jobStatus == [3, 5]:
                    # save previos .out and .error files with another names
                    sCondorOutput_vPrevious = sCondorOutput_to_use.replace('.out', '_v%d.out' % (iJobSubmission-1))
                    sCondorError_vPrevious  = sCondorError_to_use.replace('.error', '_v%d.error' % (iJobSubmission-1))
                    os.rename(condor_hadd.sCondorOutput_to_use, sCondorOutput_vPrevious)
                    os.rename(condor_hadd.sCondorError_to_use,  sCondorError_vPrevious)
                increaseJobFlavour = False
                if jobStatus == 4:
                    increaseJobFlavour = True
                condor_hadd.writeCondorSumitFile(increaseJobFlavour)
            else:
                assert(jobStatus == 2)
                # job is either running or succeeded
                continue
            if jobStatus in [0,3,4,5]:
                if run_mode == 'condor':
                    cmd1 = "condor_submit %s" % condor_hadd.sCondorSubmit_to_use
                    print("Now:  %s " % cmd1)
                    if not dryRun:
                        os.system(cmd1)
                    else:
                        pass
    if len(condor_hadd.OpRootFiles_Target) == len(condor_hadd.OpRootFiles_Exist):
        break
    else:
        time.sleep( ResubWaitingTime * 60 )
        iJobSubmission += 1

with open(condor_hadd.fJobSubLog, 'a') as fJobSubLog:
        fJobSubLog.write('%s \t Jobs are done for hadd_stage2. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
print('%s \t Jobs are done for hadd_stage2. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
sys.exit()
if 1:#do_hadd:
    if run_file == 'count_genweight':
        outputfile = 'hadd_countgenWeight.root'
    elif run_file == 'stitch':
        outputfile = 'hadd_stitch.root'
    else:
        outputfile = 'hadd_apply_stitch.root'
    cmd = f'python3 hadd.py -p {DestinationDir} -o {outputfile}'
    os.system(cmd)
