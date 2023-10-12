from glob import glob
import os
import argparse

parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
parser.add_argument('-p', dest='input_path', type=str, required=True)
parser.add_argument('-o', dest='output_file', type=str, required=True)
parser.add_argument('-w', dest='wildcard', type=str, default='ana')
args=parser.parse_args()

p = args.input_path
output_file = args.output_file
wc = args.wildcard
files = glob(f'{p}/{wc}*root')
nfiles=5
dohadd = True
assert(len(files))
count = 0
final_hadd = ''
outputdir = output_file[:output_file.rfind('/')]
rmfiles = []
while dohadd:
    outputs = []
    for idx, f in enumerate(range(0, len(files), nfiles)):
        files_tohadd = ' '.join(files[f:f+nfiles])
        output = os.path.join(outputdir, f'hadd_{count}_{idx}.root') if len(files) > nfiles else os.path.join(outputdir,output_file)
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
        rmfiles.extend(outputs)

for file_ in rmfiles:
    os.system(f'rm {file_}')
