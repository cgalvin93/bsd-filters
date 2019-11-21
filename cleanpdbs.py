#ipython ~/desktop/bsd_filters/cleanpdbs.py *.pdb

import sys
for filename in sys.argv[2:]:
    f=open(filename,'r')
    lines=[line for line in f.readlines() if line[0:4]=='ATOM' or line[0:3]=='TER' or line[0:6]=='HETNAM']
    ofile=open('clean_'+filename,'w')
    for line in lines:
        ofile.write(line)
    f.close()
    ofile.close()
