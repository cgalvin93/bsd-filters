#time ipython ~/desktop/bsd_filters/pareto.py bre-1f4p-all-score.sc ~/desktop/prjk/analysis/strc/native_seqs/1f4p_binding_site.txt 1f4p

import numpy as np
import sys
import pandas as pd
import Bio
from Bio import SeqUtils
from Bio.SeqUtils import IUPACData
import scipy
from scipy import spatial

ptn_pdb_id=sys.argv[3]
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


#store lines of scorefile
scorefile=open(sys.argv[1],'r')
lines=[line for line in scorefile.readlines()[1:]]
scorefile.close()
#acquire the names of the metrics calculated in the scorefile from first line
terms=[]
for i in lines[0].split()[1:]:
    terms.append(i)
del lines[0]

#parse the score file to create a list where each element is (structure name,property name, property value)
#for all structures and properties in the score file
strc_names=[]
datas=[]
for line in lines:
    properties=line.split()[1:]
    name=properties[-1]
    strc_names.append(name)
    for index,term in enumerate(terms[:-1]):
        val=properties[index]
        datas.append((name,term,val))

df=pd.DataFrame()
rel_terms=['hbond_sc','ref','exphyd']
rel_terms.sort()
df['terms']=pd.Series(rel_terms)
names=[a for a,b,c in datas]
names=list(set(names))
for name in names:
    vals = [c for a,b,c in datas if a==name and b in rel_terms]
    df[name]=pd.Series(vals)

scores=[]
for name in names:
    x,y,z = df[name]
    scores.append((x,y,z))
scores=np.array(scores)

def pareto(points):
    size=points.shape[0]
    ids=np.arange(size)
    pareto_front=np.ones(size,dtype=bool)
    for i in range(size):
        for j in range(size):
            if all(points[j]<=points[i])and any(points[j]<points[i]):
                pareto_front[i]=0
                break
    return ids[pareto_front]

pareto_indices=pareto(scores)
pareto_opt=[i[:-5]+'.pdb' for i in np.array(names)[pareto_indices]]

def threeto1(s):
    olseq=''
    for i in range(0,len(s),3):
        res = s[i:i+3]
        threelc=res[0]+res[1:3].lower()
        onelc = Bio.SeqUtils.IUPACData.protein_letters_3to1[threelc]
        olseq+=onelc
    return olseq

#from a list of pdb files, return binding site sequences from those files
def return_bs_seq(pdbfiles):
    seqs=[]
    for file in pdbfiles:
        fullseq=[]
        f=open(file,'r');lines=[line for line in f.readlines()]
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

def p_dist(seqs):
    n_positions=len(seqs[0])
    n_seqs=len(seqs)
    a=np.zeros(shape=(20,n_positions))
    ala=[];asp=[];asn =[];arg =[];cys=[];phe =[];gly=[];glu =[];gln =[];his =[];
    ile =[];leu=[];lys =[];met =[];pro =[];ser =[];trp =[];tyr =[];thr=[];val=[]
    for i in range(0,n_positions):
        nala=1;nasp=1;nasn=1;narg=1;ncys=1;nphe=1;ngly=1;nglu=1;ngln=1;nhis=1;
        nile=1;nleu=1;nlys=1;nmet=1;npro=1;nser=1;ntrp=1;ntyr=1;nthr=1;nval=1;
        ngap=0
        for seq in seqs:
            aa=seq[i]
            if aa=='A':
                nala+=1
            elif aa=='C':
                ncys+=1
            elif aa=='D':
                nasp+=1
            elif aa=='E':
                nglu+=1
            elif aa=='F':
                nphe+=1
            elif aa=='G':
                ngly+=1
            elif aa=='H':
                nhis+=1
            elif aa=='I':
                nile+=1
            elif aa=='K':
                nlys+=1
            elif aa=='L':
                nleu+=1
            elif aa=='M':
                nmet+=1
            elif aa=='N':
                nasn+=1
            elif aa=='P':
                npro+=1
            elif aa=='Q':
                ngln+=1
            elif aa=='R':
                narg+=1
            elif aa=='S':
                nser+=1
            elif aa=='T':
                nthr+=1
            elif aa=='V':
                nval+=1
            elif aa=='W':
                ntrp+=1
            elif aa=='Y':
                ntyr+=1
            else:
                ngap+=1
        ntot=float((n_seqs-ngap)+20);fala=nala/ntot
        fasp=nasp/ntot;fasn=nasn/ntot;farg=narg/ntot;
        fcys=ncys/ntot;fphe=nphe/ntot;fgly=ngly/ntot;
        fglu =nglu/ntot;fgln=ngln/ntot;fhis=nhis/ntot;
        ffile=nile/ntot;fleu=nleu/ntot;flys=nlys/ntot;
        fmet=nmet/ntot;fpro=npro/ntot;fser=nser/ntot;
        ftrp=ntrp/ntot;ftyr=ntyr/ntot;fthr=nthr/ntot;
        fval=nval/ntot
        ala.append(fala);asp.append(fasp);asn.append(fasn);arg.append(farg);
        cys.append(fcys);phe.append(fphe);gly.append(fgly);glu.append(fglu);
        gln.append(fgln);his.append(fhis);ile.append(ffile);leu.append(fleu);
        lys.append(flys);met.append(fmet);pro.append(fpro);ser.append(fser);
        trp.append(ftrp);tyr.append(ftyr);thr.append(fthr);val.append(fval)
    a[0]=ala;a[1]=asp;a[2]=asn;a[3]=arg;a[4]=cys;a[5]=phe;a[6]=gly
    a[7]=glu;a[8]=gln;a[9]=his;a[10]=ile;a[11]=leu;a[12]=lys;a[13]=met
    a[14]=pro;a[15]=ser;a[16]=trp;a[17]=tyr;a[18]=thr;a[19]=val
    return a

def calc_pps(p,q):
    n_pos = p.shape[1]
    pps_scores=[]
    for i in range(0,n_pos):
        jsd = (scipy.spatial.distance.jensenshannon(p[:,i],q[:,i]))**2
        pps=1-jsd
        pps_scores.append(pps)
    return pps_scores

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

native_seqs=capture_seqs(sys.argv[2])
p_nat=p_dist(native_seqs)
pareto_seqs=return_bs_seq(pareto_opt)
p_pareto=p_dist(pareto_seqs)
pps_scores=calc_pps(p_nat,p_pareto)
print pps_scores
print np.median(pps_scores)
