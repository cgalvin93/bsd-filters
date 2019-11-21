# MUTATION PROFILES
# take starting strc that designs were generated off
# for each position, look at all the mutations (ala to glu; ala to ile; etc.)
# generate matrix where y is all residues for native protein,
# x is all residue for design proteins
# entries are number of times a y->x mutation occurred (5 ala to glus? 10 ala to ile?)
# get a mutation profile for each structure, sum for all strc for each method
# compare with mutation profile of native -> all other natives
# 1. mut profiles for each strc native->des
# 2. ^ native->all natives
# 3. sum each for all strc to give des profiles for each method to compare with nat
# (nat profile same for all methods for a given strc), summed profile only 1
# 4. y all amino acids
#   x native, fd, cm, bre
#   values are number of mutations to that residue
# ^^^^can do one of these for each structure, at the least, which should def be valid
# ^^^^can maybe do another summed over all strc for each method, even though kinda sketch

#load str of starting strc bs seq
#load lists of strings for design bs seqs, native bs seqs
#for each bs position, record 2 letter sequence of mutations (eg AG,KR) for all
#seqs in des, native lists
#i dont think the pos rly matters,so i can just have one very long list for all
#muts (2 lists actually, 1 for des 1 for nat)
#---------------------------------------------------

#actually the mutations to x table should be really easy because I can just use
#frequency matrix for des and nats, cus thats all it really is
#doing that first then:
#time ipython aa_freq_table.py ptn_name native_seqs.txt design_pdbs.txt
#to be run from 'strc' directory, which contains 3 subdir, 1 for each method
#makes it easier to navigate to each of these to get design strc names
import sys
import pandas as pd
import Bio
from Bio import SeqUtils
from Bio.SeqUtils import IUPACData
#must put binding site res numbers here manually
#store pdb numbering of binding site residues for query protein
ptn_pdb_id=sys.argv[1]
if ptn_pdb_id == '2xbn':
    binding_site_res=[73, 133, 134, 135, 138, 159,161, 202, 231, 233,
                      234, 262,264, 265]
elif ptn_pdb_id == '1f4p':
    binding_site_res=[10, 11, 12, 13, 14, 15, 58, 59,
                      60, 61, 62, 68, 93, 94, 95, 98,
                      100, 101, 102]
elif ptn_pdb_id == '1zk4':
    binding_site_res=[13, 14, 15, 16, 17, 18, 19, 36,
                      37, 38, 62, 63, 89, 90, 91, 92,
                      112, 140, 142, 155, 159]
elif ptn_pdb_id == '3dk9':
    binding_site_res=[26, 27, 28, 29, 30, 31, 49, 50,
                      51, 52, 56, 57, 62,
                      66, 129, 130, 155, 156, 157,
                      177, 181, 197, 198, 201, 202,
                      291, 294, 298, 330, 331]
elif ptn_pdb_id == '3dlc':
    binding_site_res=[48, 49, 50,
                      51, 52, 53, 55, 72, 73, 74, 77,
                      100, 101, 102, 117, 118, 119,
                      122, 123]
elif ptn_pdb_id == '3r2q':
    binding_site_res=[9, 10, 11, 33, 34, 48, 49, 50,
                      62, 63, 64]
else:
    print 'the binding site residues of the query protein cannot be identified'

#convert 3 letter aa code to 1 letter
def threeto1(s):
    olseq=''
    for i in range(0,len(s),3):
        res = s[i:i+3]
        threelc=res[0]+res[1:3].lower()
        onelc = Bio.SeqUtils.IUPACData.protein_letters_3to1[threelc]
        olseq+=onelc
    return olseq

#grab binding site seqs from list of pdb files
def return_bs_seq(path,pdbfiles):
    seqs=[]
    for file in pdbfiles:
        fullseq=[]
        f=open(path+'/'+file,'r');lines=[line for line in f.readlines()]
        for line in lines:
            if line[0:4]=='ATOM':
                resname=line[17:21].strip()
                resnum = line[22:27].strip()
                fullseq.append((resnum, resname))
        fullseq=list(set(fullseq))
        s=''
        for i in binding_site_res:
            for num, name in fullseq:
                if i==int(num):
                    s+=name
        s2=threeto1(s)
        seqs.append(s2)
    return seqs

#get all the sequences from a fasta file
def capture_seqs(file):
    indices=[]
    f=open(file,'r')
    lines=[line for line in f.readlines()]
    f.close()
    for line in lines:
        if line[0]=='>':
            indices.append(lines.index(line))
    strs=[]
    for i in range(0,len(indices)-1):
        start=indices[i]
        end=indices[i+1]
        s=''
        for line in lines[start+1:end]:
            s+=str(line.strip('\n'))
        strs.append(s.upper())
    return list(set(strs))

#produce a list of (resname, res frequency)
def aa_freq_table(seqs):
    ala=[];asp=[];asn =[];arg =[];cys=[];phe =[];gly=[];glu =[];gln =[];his =[];
    ile =[];leu=[];lys =[];met =[];pro =[];ser =[];trp =[];tyr =[];thr=[];val=[]
    gap=[]
    n_seqs=float(len(seqs))
    for seq in seqs:
        for aa in seq:
            if aa=='A':
                ala.append(aa)
            elif aa=='C':
                cys.append(aa)
            elif aa=='D':
                asp.append(aa)
            elif aa=='E':
                glu.append(aa)
            elif aa=='F':
                phe.append(aa)
            elif aa=='G':
                gly.append(aa)
            elif aa=='H':
                his.append(aa)
            elif aa=='I':
                ile.append(aa)
            elif aa=='K':
                lys.append(aa)
            elif aa=='L':
                leu.append(aa)
            elif aa=='M':
                met.append(aa)
            elif aa=='N':
                asn.append(aa)
            elif aa=='P':
                pro.append(aa)
            elif aa=='Q':
                gln.append(aa)
            elif aa=='R':
                arg.append(aa)
            elif aa=='S':
                ser.append(aa)
            elif aa=='T':
                thr.append(aa)
            elif aa=='V':
                val.append(aa)
            elif aa=='W':
                trp.append(aa)
            elif aa=='Y':
                tyr.append(aa)
            else:
                gap.append(aa)
    nala=len(ala)/n_seqs;nasp=len(asp)/n_seqs;nasn=len(asn)/n_seqs;
    narg=len(arg)/n_seqs;ncys=len(cys)/n_seqs;nphe=len(phe)/n_seqs;
    ngly=len(gly)/n_seqs;nglu=len(glu)/n_seqs;ngln=len(gln)/n_seqs;
    nhis=len(his)/n_seqs;nile=len(ile)/n_seqs;nleu=len(leu)/n_seqs;
    nlys=len(lys)/n_seqs;nmet=len(met)/n_seqs;npro=len(pro)/n_seqs;
    nser=len(ser)/n_seqs;ntrp=len(trp)/n_seqs;ntyr=len(tyr)/n_seqs;
    nthr=len(thr)/n_seqs;nval=len(val)/n_seqs;
    freqs=[('ala',nala),('asp',nasp),('asn',nasn),('arg',narg),('cys',ncys),('phe',nphe),
          ('gly',ngly),('glu',nglu),('gln',ngln),('his',nhis),('ile',nile),('leu',nleu),
          ('lys',nlys),('met',nmet),('pro',npro),('ser',nser),('trp',ntrp),('tyr',ntyr),
          ('thr',nthr),('val',nval)]
    return freqs

#load native seqs
native_seqs=capture_seqs(sys.argv[2])
#load bre seqs
bre_path='bre/'+ptn_pdb_id+'/pdblist.txt'
bre_list=[]
bre_file=open(bre_path,'r')
for line in bre_file.readlines():
    bre_list.append(line.strip('\n'))
bre_seqs=return_bs_seq(bre_path[:-12],bre_list)
#load cm seqs
cm_path='cm/'+ptn_pdb_id+'/pdblist.txt'
cm_list=[]
cm_file=open(cm_path,'r')
for line in cm_file.readlines():
    cm_list.append(line.strip('\n'))
cm_seqs=return_bs_seq(cm_path[:-12],cm_list)
#load fd seqs
fd_path='fd/'+ptn_pdb_id+'/pdblist.txt'
fd_list=[]
fd_file=open(fd_path,'r')
for line in fd_file.readlines():
    fd_list.append(line.strip('\n'))
fd_seqs=return_bs_seq(fd_path[:-12],fd_list)

#get frequency table data
native_data=aa_freq_table(native_seqs)
bre_data=aa_freq_table(bre_seqs)
cm_data=aa_freq_table(cm_seqs)
fd_data=aa_freq_table(fd_seqs)

#put everything into a table
df=pd.DataFrame()
res_names=[]
nat_freqs=[]
for i,x in native_data:
    res_names.append(i)
    nat_freqs.append(x)
df['res']=pd.Series(res_names)
df['nat']=pd.Series(nat_freqs)
bre_freqs=[]
for i,x in bre_data:
    bre_freqs.append(x)
df['bre']=pd.Series(bre_freqs)
cm_freqs=[]
for i,x in cm_data:
    cm_freqs.append(x)
df['cm']=pd.Series(cm_freqs)
fd_freqs=[]
for i,x in fd_data:
    fd_freqs.append(x)
df['fd']=pd.Series(fd_freqs)
print df

#print it all to a csv file
df.to_csv('freq_table.csv')

#time ipython ~/desktop/bsd_filters/aa_freq_table.py 1f4p ~/desktop/prjk/analysis/native_seqs/1f4p_binding_site.txt pdblist.txt
#next, to do for all strc, + one table for all strc together