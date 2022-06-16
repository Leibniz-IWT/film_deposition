#!/usr/bin/python3

import sys,argparse,os

help = """%s [options]\n\n """
SCRIPT_NAME = os.path.basename(sys.argv[0])
parser = argparse.ArgumentParser(help % SCRIPT_NAME)
parser.add_argument('file', help='Enter the file you want to simulate (and multiply)')
parser.add_argument('files',type=int,help='Enter the amount of files that will be created')
parser.add_argument('-np','--ncpu',type=int,default = 2,help="The amount of cpu per simulation; default: 2")
parser.add_argument('-em','--email',default=False,help='Send email after completion;',action='store_true')
args = parser.parse_args()

with open(args.file) as fpr:
    foo = fpr.readlines()
    fpr.close()

name = foo[13].split('=')[1][2:-2]

for fil in range(args.files):
    if fil < 9:
        file_name = name + '_00'+str(fil+1)
        ff = 'DLA_batch_00'+str(fil+1)+'.py'
    elif fil < 99:
        file_name = name + '_0'+str(fil+1)
        ff = 'DLA_batch_0'+str(fil+1)+'.py'
    else:
        file_name = name + '_'+str(fil+1)
        ff = 'DLA_batch_'+str(fil+1)+'.py'

    for lin in range(len(foo)):
        if foo[lin].startswith('experiment ='):
            foo[lin] = ("experiment = '%s'\n" %file_name)

    with open(ff,'w') as fpw:
        for lin in range(len(foo)):
             fpw.write(foo[lin])
        fpw.close()

        folder = './' + file_name +'/'
        if not os.path.exists(folder):
           os.makedirs(folder)
        if args.email:
            os.system('dla_mpi %s -np %i -em' %(ff,args.ncpu))
        else:
            os.system('dla_mpi %s -np %i' %(ff,args.ncpu))
        os.system('mv %s %s' %(ff,folder))

