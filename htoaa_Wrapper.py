import os
import sys
import json
import glob
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

sConfig           = "config_htoaa.json"
sRunCommandFile   = "1_RunCommand.txt"
sJobSubLogFile    = "1_JobSubmission.log"
sOpRootFile       = "analyze_htoaa_$SAMPLE_$STAGE_$IJOB.root"

printLevel = 3

class analysis_wrapper():
    def __init__(self,
                 era, run_mode, nFilesPerJob,
                 selSamplesToRun,
                 anaVersion,
                 nResubmissionMax,
                 ResubWaitingTime,
                 iJobSubmission,
                 dryRun, run_file,
                 applystitching,
                 sAnalysis,
                 lepton_selection,
                 msoftdrop
        ):
        self.era = era
        self.run_mode = run_mode
        self.nFilesPerJob = nFilesPerJob
        self.selSamplesToRun = selSamplesToRun,
        self.anaVersion = anaVersion,
        self.nResubmissionMax = nResubmissionMax
        self.ResubWaitingTime = ResubWaitingTime
        self.iJobSubmission = iJobSubmission
        self.dryRun = dryRun
        self.run_file = run_file
        self.applystitching = applystitching
        self.sAnalysis = sAnalysis
        self.lepton_selection = lepton_selection
        self.msoftdrop = msoftdrop

        self.pwd = os.getcwd()
        self.DestinationDir = "./%s" % (anaVersion)
        if "/afs/cern.ch/" in self.pwd: self.DestinationDir = "/afs/cern.ch/work/s/snandan/public/myforkhtoaa/htoaa-1/%s" % (anaVersion)
        self.sFileRunCommand = "%s/%s" % (self.DestinationDir, sRunCommandFile)
        self.sFileJobSubLog  = "%s/%s" % (self.DestinationDir, sJobSubLogFile)
        if not os.path.exists(self.DestinationDir): os.mkdir(self.DestinationDir)

        # save run command into a .txt tile
        with open(self.sFileRunCommand, 'a') as fRunCommand:
            datatime_now = datetime.now()
            sCommand = ' '.join(sys.argv)
            fRunCommand.write('%s %s \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), sCommand))

        samplesList = None
        samplesInfo = None
        Luminosity  = None
        if self.era == Era_2018:
            samplesList = 'Samples_2018UL.json'
            with open(samplesList) as fSamplesInfo:
                self.samplesInfo = json.load(fSamplesInfo)
            Luminosity = Luminosities[era][0]

        selSamplesToRun_list = []
        if selSamplesToRun:
            self.selSamplesToRun_list = selSamplesToRun.split(',')
        self.run_analysisJobs()
        if self.dryRun: sys.exit()
        self.run_haddJobs(stage=1)
        self.run_haddJobs(stage=1.5)
        if 'data_obs' in self.selSamplesToRun_list:
            self.run_addBackgrounds()

    def run_analysisJobs(self):
        condor_ana = condor(self.sAnalysis, self.sFileJobSubLog)
        iJobSubmission = 0
        while iJobSubmission <= self.nResubmissionMax:
            if self.dryRun and iJobSubmission >0: break
            condor_ana.initialize_list()
            print('\n\n%s \t Startiing iJobSubmission for analysis: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
            #for leptonselection in ['Tight']:
            #for met in ['met_GE20', 'met_L20']:
            ##condor_ana.initialize_list()
            samples_sum_events = calculate_samples_sum_event(self.samplesInfo)
            for dbsname, sampleInfo in self.samplesInfo.items():
                sample_category = sampleInfo['sample_category']
                process_name = sampleInfo['process_name']
                if len(self.selSamplesToRun_list) > 0:
                    skipThisSample = True
                    for selSample in self.selSamplesToRun_list:
                        if 'TTG' in process_name: break
                        #if 'Enriched' not in process_name: break
                        if not self.applystitching:
                            if sample_category == 'WJets' and (re.findall('^W[1-4]Jets', process_name) or ('HT-' in process_name)):
                                break
                            #if leptonselection == 'Fake' and 'SUSY' in process_name:
                            #break
                        if selSample in process_name or selSample in sample_category:
                            skipThisSample= False
                            break
                if skipThisSample:
                    continue

                fileList = sampleInfo[sampleFormat]
                files = []
                for iEntry in fileList:
                    if "*" in iEntry:
                        files.extend( glob.glob( iEntry ) )
                    else:  files.append( iEntry )
                sample_cossSection = sampleInfo["cross_section"] if (sample_category != kData) else None
                sample_nEvents     = sampleInfo["nEvents"]
                sample_exist = check_sample(process_name, samples_sum_events)
                if sample_exist[0]:
                    sample_sumEvents   = sample_exist[1]
                else:
                    sample_sumEvents   = sampleInfo["sumEvents"] if (sample_category != kData) else None
                nSplits = int(len(files) / self.nFilesPerJob) if self.nFilesPerJob > 0 else 1
                nSplits = 1 if nSplits==0 else nSplits
                files_splitted = np.array_split(files, nSplits)
                sampledir = os.path.join(self.DestinationDir, process_name)
                for iJob in range(len(files_splitted)):
                    #config = config_Template.deepcopy()
                    config = copy.deepcopy(config_Template)
                    config['lepton_selection'] = self.lepton_selection
                    config['msoftdrop'] = self.msoftdrop
                    # Job related files
                    outputdir = os.path.join(sampledir, 'output')
                    sOpRootFile_to_use = sOpRootFile.replace('$SAMPLE', f'{process_name}')
                    sOpRootFile_to_use = sOpRootFile_to_use.replace('$STAGE', str(0))
                    sOpRootFile_to_use = sOpRootFile_to_use.replace('$IJOB', str(iJob))
                    self.init_condor_files(condor_ana, outputdir, sOpRootFile_to_use)
                    jobStatusForJobSubmission  = [0, 3, 4, 5]
                    jobStatus = condor_ana.check_jobstatus(printLevel >= 3)
                    if jobStatus == 0:
                        config["era"] = self.era
                        config["inputFiles"] = list( files_splitted[iJob] )
                        config["outputFile"] = condor_ana.sOpRootFile_to_use
                        config["sampleCategory"] = sample_category
                        config['process_name'] = process_name
                        config["isMC"] = 1 if (sample_category != kData) else 0
                        config["leptonselection"] = self.lepton_selection
                        config["met"] = self.msoftdrop
                        config["nEvents"] = sample_nEvents
                        if (sample_category != kData):
                            config["crossSection"] = sample_cossSection
                            config["sumEvents"] = sample_sumEvents
                            if self.applystitching:
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
                    elif jobStatus in [1, 2]:
                        # job is either running or succeeded
                        continue
                    increaseJobFlavour = False
                    self.submitJob(condor_ana, jobStatus)
            # write JobSubmission status report in JobSubLog file
            condor_ana.update_JobSubLog(iJobSubmission)
            if self.dryRun:
                print('%s \t druRun with iJobSubmission: %d  \nTerminating...\n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
            elif len(condor_ana.OpRootFiles_Target) == len(condor_ana.OpRootFiles_Exist):
                break
            else:
                time.sleep( self.ResubWaitingTime * 60 )
            iJobSubmission += 1

        with open(condor_ana.fJobSubLog, 'a') as fJobSubLog:
            fJobSubLog.write('%s \t Jobs are done. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
        print('%s \t Jobs are done. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
        del condor_ana

    def run_haddJobs(self, stage):
        condor_hadd = condor('hadd.py', sJobSubLogFile)
        iJobSubmission = 0
        while iJobSubmission <= self.nResubmissionMax:
            print('\n\n%s \t Startiing iJobSubmission for hadd: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission)) 
            condor_hadd.initialize_list()
            #for leptonselection in ['Tight']:
            #for met in ['met_GE20', 'met_L20']:
            if stage == 1:
                inputdirs = glob.glob(f'{self.DestinationDir}/*/output/')
            else:
                inputdirs = [os.path.join(self.DestinationDir)]
                #if not len(glob.glob(f'{inputdirs}/*')): continue
            for inputdir in inputdirs:
                if stage == 1:
                    process_name = inputdir.split('/')[-3]
                    #if 'SUSY' in process_name and leptonselection == 'Fake': continue
                    if process_name == 'hadd': continue
                    outputdir = inputdir
                    sOpRootFile_to_use = f'hadd_{process_name}.root'
                else:
                    outputdir = os.path.join(inputdir, 'hadd', 'output')
                    sOpRootFile_to_use = f'hadd_stage1.root'
                self.init_condor_files(condor_hadd, outputdir, sOpRootFile_to_use)
                jobStatus = condor_hadd.check_jobstatus(printLevel >= 3)
                config = {}
                if jobStatus == 0:
                    config["inputdir"] = inputdir
                    config["outputFile"] = os.path.join(outputdir, sOpRootFile_to_use)
                    if stage == 1:
                        config['wildcard'] = 'ana'
                    else:
                        config['wildcard'] = '*/output/hadd'
                    with open(condor_hadd.sConfig_to_use, "w") as fConfig:
                        json.dump( config,  fConfig, indent=4)
                        condor_hadd.writeCondorExecFile(f'-p {config["inputdir"]} -o {config["outputFile"]} -w {config["wildcard"]}')
                elif jobStatus in [1, 2]:
                    # job is either running or succeeded
                    continue
                self.submitJob(condor_hadd, jobStatus)
            # write JobSubmission status report in JobSubLog file
            condor_hadd.update_JobSubLog(iJobSubmission)
            if self.dryRun:
                print('%s \t druRun with iJobSubmission: %d  \nTerminating...\n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
            elif len(condor_hadd.OpRootFiles_Target) == len(condor_hadd.OpRootFiles_Exist):
                break
            else:
                time.sleep( self.ResubWaitingTime * 60 )
            iJobSubmission += 1
        with open(condor_hadd.fJobSubLog, 'a') as fJobSubLog:
            fJobSubLog.write(f'%s \t Jobs are done for hadd stage_{stage}. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
        print(f'%s \t Jobs are done for hadd stage_{stage} iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))

    def init_condor_files(self, con, outputdir, sOpRootFile_to_use):
        if not os.path.exists(outputdir): os.makedirs(outputdir)
        sOpRootFile_to_use    = '%s/%s' % (outputdir, sOpRootFile_to_use)
        configdir = outputdir.replace('output', 'config')
        if not os.path.exists(configdir): os.makedirs(configdir)
        sConfig_to_use = sOpRootFile_to_use.replace('.root', '_config.json').replace(outputdir, configdir)
        condordir = outputdir.replace('output', 'condor')
        if not os.path.exists(condordir): os.makedirs(condordir)
        sCondorExec_to_use = sOpRootFile_to_use.replace('.root', '_condor_exec.sh').replace(outputdir, condordir)
        sCondorSubmit_to_use = sOpRootFile_to_use.replace('.root', '_condor_submit.sh').replace(outputdir, condordir)
        sCondorLog_to_use = sOpRootFile_to_use.replace('.root', '_condor.log').replace(outputdir, condordir)
        sCondorOutput_to_use = sOpRootFile_to_use.replace('.root', '_condor.out').replace(outputdir, condordir)
        sCondorError_to_use = sOpRootFile_to_use.replace('.root', '_condor.error').replace(outputdir, condordir)
        con.initialize_files(sOpRootFile_to_use, configdir,
            sConfig_to_use, condordir, sCondorExec_to_use,
            sCondorSubmit_to_use, sCondorLog_to_use, sCondorOutput_to_use,
            sCondorError_to_use
        )
    def submitJob(self, con, jobStatus):
        increaseJobFlavour = False
        if jobStatus == [3, 5]:
            # save previos .out and .error files with another names
            sCondorOutput_vPrevious = con.sCondorOutput_to_use.replace('.out', '_v%d.out' % (iJobSubmission-1))
            sCondorError_vPrevious  = con.sCondorError_to_use.replace('.error', '_v%d.error' % (iJobSubmission-1))
            os.rename(con.sCondorOutput_to_use, sCondorOutput_vPrevious)
            os.rename(con.sCondorError_to_use,  sCondorError_vPrevious)
        elif jobStatus == 4:
            increaseJobFlavour = True
        con.writeCondorSumitFile(increaseJobFlavour)
        if self.run_mode == 'condor':
            cmd1 = "condor_submit %s" % con.sCondorSubmit_to_use
            print("Now:  %s " % cmd1)
            if not self.dryRun:
                os.system(cmd1)
        else:
            pass
    
    def run_addBackgrounds(self):
        
        samples = [s for s in self.selSamplesToRun_list if ('SUSY' not in s and 'data_obs' not in s)]
        samples = ' '.join(samples)
        for leptonselection in ['Tight']:
            for met in ['met_L20', 'met_GE20']:#, 'met_L20']:
                ip = os.path.join(self.DestinationDir, leptonselection, met, 'hadd', 'output')
                op = os.path.join(self.DestinationDir, leptonselection, met, 'addBackgrounds')
                os.makedirs(op, exist_ok=True)
                ifile = os.path.join(ip, f'hadd_stage1_{leptonselection}_{met}.root')
                ofile = os.path.join(op, f'addBackground_{leptonselection}_{met}.root')
                cmd = f'python3 addBackgrounds.py -p {ifile} -o {ofile} -s {samples} -syst topt_reweigingUp topt_reweigingDown'
                print('cmd: ', cmd)
                os.system(cmd)
'''
if 1:#do_hadd:
    if run_file == 'count_genweight':
        outputfile = 'hadd_countgenWeight.root'
    elif run_file == 'stitch':
        outputfile = 'hadd_stitch.root'
    else:
        outputfile = 'hadd_apply_stitch.root'
    cmd = f'python3 hadd.py -p {DestinationDir} -o {outputfile}'
    os.system(cmd)
'''
