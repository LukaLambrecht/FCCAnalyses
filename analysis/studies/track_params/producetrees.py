import os
import sys
import glob
import time
import argparse
import subprocess
import multiprocessing as mp
from subprocess import Popen,PIPE
from datetime import date


def read_samplelist(samplelist):
    '''
    Read a sample list.
    Input arguments:
    - samplelist: path to txt file listing samples.
      every line in the sample list is assumed to be a path to a root file
      (potentially containing unix-style wildcards).
    '''
    # read each line from the sample list
    with open(samplelist, 'r') as f:
        lines = f.readlines()
    # cleaning and filtering
    lines = [l.strip(' \t\n') for l in lines]
    lines = [l for l in lines if (l.endswith('.root') and not l.startswith('#'))]
    # expand wildcards
    input_files = sum([glob.glob(l) for l in lines], [])
    return input_files


if __name__ == '__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, nargs='+',
      help='Input .root files, OR path to a .txt file listing input .root files (one per line)')
    parser.add_argument('-o', '--outputfile', required=True,
      help='Output .root file')
    parser.add_argument('-n', '--nevents', type=int, default=-1,
      help='Number of events to process (default: all available events)')
      # note: when multiple files are provided, this number seems to be the total number of events,
      #       not the number of events per file!
      #       so some files may not be read at all if nevents is reached in an earlier file.
    args = parser.parse_args()

    # find input files
    input_files = []
    for el in args.input:
        if el.endswith('.root'): input_files.append(el)
        elif el.endswith('.txt'):
            input_files += read_samplelist(el)
    print(f'Found following input files ({len(input_files)}):')
    for f in input_files: print(f'  - {f}')

    # check output file
    # (naming depends on output file ending in '.root', maybe make more robust/flexible later)
    if not args.outputfile.endswith('.root'):
        msg = 'Output file must end with ".root"'
        raise Exception(msg)

    # make output directory
    outputdir = os.path.dirname(args.outputfile)
    if not os.path.exists(outputdir): os.makedirs(outputdir)

    # make command for stage 1
    tempfile = args.outputfile
    cmd_stage1 = 'fccanalysis run analysis.py'
    if args.nevents > 0: cmd_stage1 += ' --nevents {}'.format(args.nevents)
    cmd_stage1 += ' --output {}'.format(tempfile)
    cmd_stage1 += ' --files-list {}'.format(' '.join(input_files))

    # run stage 1
    print(f'Now running stage 1...')
    print(cmd_stage1)
    start_time = time.time()
    subprocess.check_call(cmd_stage1, shell=True)
    end_time = time.time()
    msg = "Stage 1 time: {:.3f}s.\n".format(end_time - start_time)
    print(msg)
