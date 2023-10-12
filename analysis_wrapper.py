from htoaa_Wrapper import analysis_wrapper
import sys
import argparse

Era_2016 = '2016'
Era_2017 = '2017'
Era_2018 = '2018'

if __name__ == '__main__':

    print("htoaa_Wrapper:: main: {}".format(sys.argv)); sys.stdout.flush()

    parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
    parser.add_argument('-era', dest='era',   type=str, default=Era_2018, choices=[Era_2016, Era_2017, Era_2018], required=False)
    parser.add_argument('-run_mode',          type=str, default='local', choices=['local', 'condor'])
    parser.add_argument('-v', '--version',    type=str, default=None, required=True)
    parser.add_argument('-samples',           type=str, default=None, help='samples to run seperated by comma')
    parser.add_argument('-nFilesPerJob',      type=int, default=8)
    parser.add_argument('-nResubMax',         type=int, default=80)
    parser.add_argument('-ResubWaitingTime',  type=int, default=5, help='Resubmit failed jobs after every xx minutes')
    parser.add_argument('-iJobSubmission',    type=int, default=0,  help='Job submission iteration. Specify previous last job submittion iteration if script terminated for some reason.')
    parser.add_argument('-dryRun',            action='store_true', default=False)
    parser.add_argument('-rf', '--run_file',    type=str, default='analysis', choices=['analysis', 'count_genweight', 'stitch'])
    parser.add_argument('-as', '--applystitching', action='store_false', default=True)
    parser.add_argument('-ls', dest='lepton_selection', type=str, required=True)
    parser.add_argument('-ms', dest='msoftdrop', type=str, nargs='+', required=True)

    args=parser.parse_args()
    print("args: {}".format(args))

    sAnalysis         = "whtoaa_Analysis_wCoffea.py"

    if args.run_file == 'count_genweight':
        sAnalysis = 'countSumEventsInSample.py'
    elif args.run_file == 'stitch':
        sAnalysis = 'count_LHE_HT.py'

    analysis = analysis_wrapper(
        era              = args.era,
        run_mode         = args.run_mode,
        nFilesPerJob     = args.nFilesPerJob,
        selSamplesToRun  = args.samples,
        anaVersion       = args.version,
        nResubmissionMax = args.nResubMax,
        ResubWaitingTime = args.ResubWaitingTime,
        iJobSubmission   = args.iJobSubmission,
        dryRun           = args.dryRun,
        run_file         = args.run_file,
        applystitching   = args.applystitching,
        sAnalysis        = sAnalysis,
        lepton_selection = args.lepton_selection,
        msoftdrop        = args.msoftdrop
    )
