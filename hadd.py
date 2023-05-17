import glob
import os
import argparse

parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
parser.add_argument('-p', dest='input_path', type=str, required=True)
parser.add_argument('-o', dest='output_file', type=str, required=True)
args=parser.parse_args()

p = args.input_path#'/afs/cern.ch/work/s/snandan/public/myforkhtoaa/htoaa-1/test_v2'
output_file = args.output_file
files = glob.glob(f'{p}/*root')
nfiles=5
dohadd = True
count = 0
final_hadd = ''
while dohadd:
    outputs = []
    for idx, f in enumerate(range(0, len(files), nfiles)):
        files_tohadd = ' '.join(files[f:f+nfiles])
        output = f'hadd_{count}_{idx}.root' if len(files) > nfiles else output_file
        outputs.append(output)
        cmd = f'hadd -f {output} {files_tohadd}'
        print('cmd: ', cmd)
        os.system(cmd)
    if len(outputs) == 1:
        dohadd = False
        final_hadd = outputs[0]
    else:
        files = outputs
        count += 1
    
