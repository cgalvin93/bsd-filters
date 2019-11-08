#ipython merge_scorefiles.py *.sc
import sys

ofile=open('all_scores.sc','w')

startfile=sys.argv[1]
f=open(startfile,'r')
lines=[line for line in f.readlines()]
for line in lines:
    ofile.write(line)
f.close()


for file in sys.argv[2:]:
    f=open(file,'r')
    flines=[line for line in f.readlines()]
    for line in flines[2:]:
        ofile.write(line)
    f.close()

ofile.close()
