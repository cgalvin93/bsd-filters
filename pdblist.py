#ipython pdblist.py pdblist.txt *.pdb
import sys

ofilename=sys.argv[1]
ofile=open(ofilename,'w')


for strc in sys.argv[2:]:
	s = str(strc).strip()
	ofile.write(s+'\n')

ofile.close()
