#writes file names of all pdb files in working directory to a text file
#python ~/desktop/bsd_filters/pdblist.py pdblist.txt *.pdb
import sys

ofilename=sys.argv[1]
ofile=open(ofilename,'w')


for strc in sys.argv[2:]:
	s = str(strc).strip()
	ofile.write(s+'\n')

ofile.close()
