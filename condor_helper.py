import os
import sys
from collections import OrderedDict as OD
from datetime import datetime

class condor():
    def __init__(self, sAnalysis, sFileJobSubLog):
        self.sAnalysis = sAnalysis
        self.fJobSubLog = sFileJobSubLog

    def initialize_list(self):
        self.OpRootFiles_Exist = []
        self.OpRootFiles_Target = []
        self.OpRootFiles_iJobSubmission = []
        self.jobStatus_dict = {}

    def initialize_files(self, sOpRootFile_to_use, configdir,
                         sConfig_to_use, condordir, sCondorExec_to_use,
                         sCondorSubmit_to_use, sCondorLog_to_use, sCondorOutput_to_use,
                         sCondorError_to_use
    ):
        self.sOpRootFile_to_use = sOpRootFile_to_use
        self.configdir = configdir
        self.sConfig_to_use = sConfig_to_use
        self.condordir = condordir
        self.sCondorExec_to_use = sCondorExec_to_use
        self.sCondorSubmit_to_use = sCondorSubmit_to_use
        self.sCondorLog_to_use = sCondorLog_to_use
        self.sCondorOutput_to_use = sCondorOutput_to_use
        self.sCondorError_to_use = sCondorError_to_use

    def writeCondorExecFile(self, sConfig_to_use=''):
        pwd = os.getcwd()
        if not os.path.isfile(self.sCondorExec_to_use):    
            with open(self.sCondorExec_to_use, 'w') as f:
                f.write("#!/bin/bash  \n\n")
                f.write("cd %s \n" % pwd)
                f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
                f.write("export SCRAM_ARCH=slc6_amd64_gcc700  \n")
                f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
                #f.write("cd ")
                #f.write("export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  \n")
                # Using x509 proxy without shipping it with the job  
                #https://batchdocs.web.cern.ch/tutorial/exercise2e_proxy.html
                f.write("export X509_USER_PROXY=$1 \n")
                f.write("voms-proxy-info -all \n")
                f.write("voms-proxy-info -all -file $1 \n")
                #f.write("eval \n")
                f.write("cd %s \n" % (pwd))
                f.write("source /afs/cern.ch/user/s/snandan/.bashrc \n")
                f.write("which conda \n")
                f.write("time conda env list \n")
                f.write("conda activate myEnv \n")
                f.write("time conda env list \n")

                f.write("time conda list \n")
                f.write("which python3 \n")
                f.write("python3 -V \n")
                #f.write(" \n")
                f.write("time /afs/cern.ch/work/s/snandan/anaconda3/envs/myEnv/bin/python3 %s/%s  %s \n" % (pwd, self.sAnalysis, sConfig_to_use if sConfig_to_use!='' else self.sConfig_to_use))

            os.system("chmod a+x %s" % self.sCondorExec_to_use)
        return;

    def writeCondorSumitFile(self, increaseJobFlavour=False):
        '''
        Job Flavours::
        espresso     = 20 minutes
        microcentury = 1 hour
        longlunch    = 2 hours
        workday      = 8 hours
        tomorrow     = 1 day
        testmatch    = 3 days
        nextweek     = 1 week
        '''

        jobFlavours = OD([
            (0, 'espresso'),
            (1, 'microcentury'),
            (2, 'longlunch'),
            (3, 'workday'),
            (4, 'tomorrow'),
            (5, 'testmatch'),
            (6, 'nextweek'),
        ])
        #iJobFlavour = 2 # 2, 'longlunch' 2 hours
        iJobFlavour = 1 # 1, 'microcentury'
        if increaseJobFlavour: iJobFlavour += 1

        #if not os.path.isfile(condor_submit_file):
        with open(self.sCondorSubmit_to_use, 'w') as f:
            f.write("universe = vanilla \n")
            f.write("executable = %s \n" % self.sCondorExec_to_use)
            f.write("getenv = TRUE \n")
            f.write("log = %s \n" % self.sCondorLog_to_use)
            f.write("output = %s \n" % self.sCondorOutput_to_use)
            f.write("error = %s \n" % self.sCondorError_to_use)
            f.write("notification = never \n")
            f.write("should_transfer_files = YES \n")
            f.write("when_to_transfer_output = ON_EXIT \n")

            #f.write("x509userproxy = /afs/cern.ch/user/s/ssawant/x509up_u108989 \n")
            #f.write("use_x509userproxy = true \n")
            f.write("x509userproxy=proxy \n")
            f.write("arguments = $(x509userproxy)  \n")

            #f.write("+JobFlavour = \"longlunch\" \n")
            f.write("+JobFlavour = \"%s\" \n" % (jobFlavours[iJobFlavour]))
            f.write("queue \n")

        os.system("chmod a+x %s" % self.sCondorSubmit_to_use)
        return

    def searchStringInFile(self, sFileName, searchString, nLinesToSearch, SearchFromEnd=True):
        lines = None
        with open(sFileName, 'r') as f:
            if SearchFromEnd:
                # https://stackoverflow.com/questions/260273/\
                    #most-efficient-way-to-search-the-last-x-lines-of-a-file
                # f.seek(offset, whence)
                # offset − This is the position of the read/write pointer within the file.
                # whence − This is optional and defaults to 0 which means absolute file positioning,\
                #other values are 1 which means seek relative to the current position and
                #2 means seek relative to the file's end.
                f.seek(0, 2)  # Seek @ EOF
                fsize = f.tell() # Get size of file
                f.seek(max(fsize-1024, 0), 0) # Set pos @ last n chars
                lines = f.readlines()  # Read to end

        lines = lines[-1*nLinesToSearch : ] # Get last x lines
        searchStringFound = False
        for line in lines:
            if searchString in line:
                searchStringFound = True
                break
        return searchStringFound

    def check_jobstatus(self, verbose):
        jobStatus = -1
        jobStatusForJobSubmission = [0, 3, 4, 5]
        isConfigExist        = os.path.isfile(self.sConfig_to_use)
        isOpRootFileExist    = os.path.isfile(self.sOpRootFile_to_use)
        isCondorExecExist    = os.path.isfile(self.sCondorExec_to_use)
        isCondorSubmitExist  = os.path.isfile(self.sCondorSubmit_to_use)
        isCondorLogExist     = os.path.isfile(self.sCondorLog_to_use)
        isCondorOutputExist  = os.path.isfile(self.sCondorOutput_to_use)
        isCondorErrorExist   = os.path.isfile(self.sCondorError_to_use)
        execfile = self.sCondorExec_to_use.split('/')[-1]
        if not isConfigExist:
            jobStatus = 0 # job not yet submitted
            if verbose:
                print(f"  jobStatus = 0 for {self.sConfig_to_use}")
        elif isOpRootFileExist:
            jobStatus = 1 # job ran successfully
            self.OpRootFiles_Exist.append(self.sOpRootFile_to_use)
            if verbose:
                print(f"  jobStatus = 1")
        elif not os.system(f'condor_q -nobatch | grep {execfile}'):
            jobStatus = 2 # job is running
            if verbose:
                print(f"  jobStatus = 2 for {self.sOpRootFile_to_use}")
        else:
            if isCondorLogExist:
                if (self.searchStringInFile(
                        sFileName       = self.sCondorLog_to_use,
                        searchString    = 'Job terminated',
                        nLinesToSearch  = 3,
                        SearchFromEnd   = True)):
                    # check wheter the job was terminated or not
                    jobStatus = 3 # job failed due to some other error
                    if verbose:
                        print(f"  jobStatus = 3")
                elif (self.searchStringInFile(
                        sFileName       = self.sCondorLog_to_use,
                        searchString    = 'Job was aborted',
                        nLinesToSearch  = 3,
                        SearchFromEnd   = True)):
                    # check wheter sCondorError does not exist due to Job was aborted
                    jobStatus = 4 # job aborted
                    if verbose:
                        print(f"  jobStatus = 4 for {self.sOpRootFile_to_use}")
                # check if job failed due to XRootD error
            elif isCondorErrorExist:
                if (self.searchStringInFile(
                        sFileName       = self.sCondorError_to_use,
                        searchString    = 'OSError: XRootD error: [ERROR]', 
                        nLinesToSearch  = 150,
                        SearchFromEnd   = True) or \
                    self.searchStringInFile(
                        sFileName       = self.sCondorError_to_use,
                        searchString    = '[ERROR] Invalid redirect URL', 
                        nLinesToSearch  = 150,
                        SearchFromEnd   = True) ):
                    jobStatus = 5 # job failed due to XRootD error
                    if verbose:
                        print(f"  jobStatus = 5")
        self.OpRootFiles_Target.append(self.sOpRootFile_to_use)
        if jobStatus in jobStatusForJobSubmission : # [0, 3, 4]:
            self.OpRootFiles_iJobSubmission.append(self.sOpRootFile_to_use)
        if jobStatus not in self.jobStatus_dict.keys():
            self.jobStatus_dict[jobStatus] = [self.sOpRootFile_to_use]
        else:
            self.jobStatus_dict[jobStatus].append(self.sOpRootFile_to_use)
        return jobStatus

    def update_JobSubLog(self, iJobSubmission):
        with open(self.fJobSubLog, 'a') as fJobSubLog:
            fJobSubLog.write('%s \t iJobSubmission %d \t OpRootFiles_Target (%d):  \n'\
                             % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), 
                                iJobSubmission, len(self.OpRootFiles_Target)))
            if iJobSubmission == 0:
                for f in self.OpRootFiles_Target:
                    fJobSubLog.write('\t %s \n' % (f))
            else:
                fJobSubLog.write('OpRootFiles_Exist %d out of %d. \n'\
                                 % (len(self.OpRootFiles_Exist), len(self.OpRootFiles_Target)))
                fJobSubLog.write('OpRootFiles_iJobSubmission (%d): '\
                                 % (len(self.OpRootFiles_iJobSubmission)))
                for f in self.OpRootFiles_iJobSubmission:
                    fJobSubLog.write('\t %s \n' % (f))
                fJobSubLog.write('\n\nJob status wise output files: \n')
                for jobStatus in self.jobStatus_dict.keys():
                    fJobSubLog.write('\t jobStatus %d (%d) \n' \
                                     % (jobStatus, len(self.jobStatus_dict[jobStatus])))
                    if jobStatus in [0, 1]: continue
                    for f in self.jobStatus_dict[jobStatus]:
                        fJobSubLog.write('\t\t %s \n' % (f))
                    
                fJobSubLog.write('%s\n\n\n' % ('-'*10))

            jobStatus_list = [ (jobStatus, len(self.jobStatus_dict[jobStatus])) for jobStatus in self.jobStatus_dict.keys() ]
            print('\n\n\n%s \t iJobSubmission %d \t OpRootFiles_Exist %d out of %d.\
            No. of jobs submitted in this resubmission: %d:  '\
                  % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), 
                     iJobSubmission, len(self.OpRootFiles_Exist), len(self.OpRootFiles_Target), len(self.OpRootFiles_iJobSubmission)))
            print(f"jobStatus_list: {jobStatus_list} \n"); sys.stdout.flush()
            #return
