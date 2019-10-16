#writes a resfile for an input pdb file, requires text file specifying which
#residues to allow packing or design, fixes current sidechain of other residues
#ipython resfile.py protein.pdb options.txt out.resfile
#options file format is desiredcode,num1chain,num2chain,etc. with diff code on newline
#ex.
#ALLAA,1A,2B
#NATAA,3A,4A,1X

import sys

ptnfile = sys.argv[1]
pdbfile = open(ptnfile, 'r')

optfile = sys.argv[2]
options = open(optfile, 'r')

outfilename = sys.argv[3]
ofile = open(outfilename,'w')

#parse pdb file to get list of str w resnum, icode, chain, for each res
rawres=[]
for line in pdbfile.readlines():
	if line[0:4]=='ATOM':
		resnum = line[22:27]
		chain = line[21:22]
		rawres.append((resnum,chain))
residues=[]
for i,x in rawres:
	if (i,x) not in residues:
		residues.append((i,x))

#parse options file to get lists of res w certain options
allaa=[]
nataa=[]
for line in options.readlines():
	if line[0:5]=='ALLAA':
		l = line.split(',')
		l2=[i.strip() for i in l[1:]]
		for i in l2:
			allaa.append(i)
	elif line[0:5]=='NATAA':
		l = line.split(',')
		l2=[i.strip() for i in l[2:]]
		for i in l2:
			nataa.append(i)

#write header
ofile.write('USE_INPUT_SC')
ofile.write('\nEX1 EX2')
ofile.write('\nstart'+'\n')

#write line for each residue + it's corresponding code
for i,x in residues:
	res_id = str(i.strip())+str(x)
	if res_id in allaa:
		s = str(i) + ' ' + str(x) + ' ALLAA'	#resnum of diff length puts char in diff place on line!!!!
		ofile.write(s+'\n')
	elif res_id in nataa:
		s = str(i) + ' ' + str(x) + ' NATAA'
		ofile.write(s+'\n')
	else:
		s = str(i) + ' ' + str(x) + ' NATRO'
		ofile.write(s+'\n')

ofile.close()

'''
ALLAA           # allow all amino acids
EX 1 EX 2       # allow extra chi rotamers at chi-id 1 and 2 (note: multiple commands can be on the same line.)
USE_INPUT_SC    # allow the use of the input side chain conformation   ( see below for more detailed description of commands)
start
<PDBNUM>[<ICODE>] <CHAIN>  <COMMANDS> 	#<PDBNUM>[<ICODE>] corresponds to columns 22-26
40A Q ALLAA   # Residue 40, insertion code A, on chain Q, use any residue type
'''
