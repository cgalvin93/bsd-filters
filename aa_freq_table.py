#time ipython ~/desktop/bsd_filters/aa_freq_table.py



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
import numpy as np

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

ptn_ids=['2xbn','1f4p','1zk4','3dk9','3r2q','3dlc','3s6f']

all_native=[]
all_fd=[]
all_bre=[]
all_cm=[]
res_names=[]

for ptn_pdb_id in ptn_ids:
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
    elif ptn_pdb_id == '3s6f':
        binding_site_res=[78, 79, 80, 85, 86, 87, 88, 90,
                          91, 114, 115, 118, 119, 121]
    else:
        print 'the binding site residues of the query protein cannot be identified'
    native_path='native_seqs/'+ptn_pdb_id.upper()+'_binding_site.txt'
    native_seqs=capture_seqs(native_path)
    bre_path='bre/'+ptn_pdb_id+'/pdblist.txt'
    bre_list=[]
    bre_file=open(bre_path,'r')
    for line in bre_file.readlines():
        bre_list.append(line.strip('\n'))
    bre_seqs=return_bs_seq(bre_path[:-12],bre_list)
    cm_path='cm/'+ptn_pdb_id+'/pdblist.txt'
    cm_list=[]
    cm_file=open(cm_path,'r')
    for line in cm_file.readlines():
        cm_list.append(line.strip('\n'))
    cm_seqs=return_bs_seq(cm_path[:-12],cm_list)
    fd_path='fd/'+ptn_pdb_id+'/pdblist.txt'
    fd_list=[]
    fd_file=open(fd_path,'r')
    for line in fd_file.readlines():
        fd_list.append(line.strip('\n'))
    fd_seqs=return_bs_seq(fd_path[:-12],fd_list)
    native_data=aa_freq_table(native_seqs)
    bre_data=aa_freq_table(bre_seqs)
    cm_data=aa_freq_table(cm_seqs)
    fd_data=aa_freq_table(fd_seqs)
#    df=pd.DataFrame()
#    res_names=[]
    nat_freqs=[]
    for i,x in native_data:
        res_names.append(i)
        nat_freqs.append(x)
#    df['res']=pd.Series(res_names)
#    df['nat']=pd.Series(nat_freqs)
    bre_freqs=[]
    for i,x in bre_data:
        bre_freqs.append(x)
#    df['bre']=pd.Series(bre_freqs)
    cm_freqs=[]
    for i,x in cm_data:
        cm_freqs.append(x)
#    df['cm']=pd.Series(cm_freqs)
    fd_freqs=[]
    for i,x in fd_data:
        fd_freqs.append(x)
#    df['fd']=pd.Series(fd_freqs)
#    ofilename=ptn_pdb_id+'_aa_freqs.csv'
#    df.to_csv(ofilename)
    all_native.append(nat_freqs)
    all_bre.append(bre_freqs)
    all_cm.append(cm_freqs)
    all_fd.append(fd_freqs)

#for getting avged results for all strc
df=pd.DataFrame()
df['res']=pd.Series(res_names)
nat_avg_freqs=[]
bre_avg_freqs=[]
cm_avg_freqs=[]
fd_avg_freqs=[]
for x in range(0,20):
    aa_nat_freq_vals=[]
    for i in all_native:
        aa_nat_freq_vals.append(i[x])
    avg_nat_aa_freq=np.mean(aa_nat_freq_vals)
    nat_avg_freqs.append(avg_nat_aa_freq)
    aa_bre_freq_vals=[]
    for i in all_bre:
        aa_bre_freq_vals.append(i[x])
    avg_bre_aa_freq=np.mean(aa_bre_freq_vals)
    bre_avg_freqs.append(avg_bre_aa_freq)
    aa_cm_freq_vals=[]
    for i in all_cm:
        aa_cm_freq_vals.append(i[x])
    avg_cm_aa_freq=np.mean(aa_cm_freq_vals)
    cm_avg_freqs.append(avg_cm_aa_freq)
    aa_fd_freq_vals=[]
    for i in all_fd:
        aa_fd_freq_vals.append(i[x])
    avg_fd_aa_freq=np.mean(aa_fd_freq_vals)
    fd_avg_freqs.append(avg_fd_aa_freq)
df['nat']=pd.Series(nat_avg_freqs)
df['bre']=pd.Series(bre_avg_freqs)
df['fd']=pd.Series(fd_avg_freqs)
df['cm']=pd.Series(cm_avg_freqs)
df.to_csv('allstrc_aa_freqs.csv')


'''
import Bio
from Bio import SeqUtils
from Bio.SeqUtils import IUPACData
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def ratios(df):
    allnewrows=[]
    for i in range(0,df.shape[0]):
        rowvals=[x for x in df.iloc[i][1:]]
        nat=rowvals[1];res=rowvals[0]
        newrow=[];newrow.append(res)
        for k in rowvals[1:]:
            ratio=k/nat
            newrow.append(ratio)
        allnewrows.append(newrow)
    newdf=pd.DataFrame(allnewrows,columns=['res','nat','bre','cm','fd'])
    return newdf

def lists(df):
    labels=[]
    nat=[]
    bre=[]
    cm=[]
    fd=[]
    for i in df['res']:
        labels.append(i)
    for i in df['nat']:
        nat.append(i)
    for i in df['bre']:
        bre.append(i)
    for i in df['cm']:
        cm.append(i)
    for i in df['fd']:
        fd.append(i)
    return labels, nat, bre, cm, fd

def threeto1(s):
    olseq=''
    for i in range(0,len(s),3):
        res = s[i:i+3].upper()
        threelc=res[0]+res[1:3].lower()
        onelc = Bio.SeqUtils.IUPACData.protein_letters_3to1[threelc]
        olseq+=onelc
    return olseq

def plot(filename):
    df=pd.read_csv(filename)
    dfr=ratios(df)
    a,b,c,d,e=lists(dfr)
    a2=[]
    for i in a:
        a2.append(threeto1(i))
    yticks=range(5)
    x=np.arange(len(a2))
    width = 0.3
    fig, ax = plt.subplots()
    rects2 = ax.bar(x- width, c, width, label='BRE',color='green')
    rects3 = ax.bar(x, d, width, label='CM',color='red')
    rects4 = ax.bar(x + width, e, width, label='FD',color='blue')
    ax.set_ylabel('Relative Abundance')
    ax.set_title('Amino Acid Frequencies Amongst Designs, Relative to Native Frequencies')
    ax.set_xticks(x)
    ax.set_yticks(yticks)
    ax.set_xticklabels(a2)
    ax.legend()
    plt.grid()
    plt.savefig(filename[:-3]+'pdf')

a,b,c,d,e=
labs,nat,bre,fd,cm

x=np.arange(4)
y=np.arange(4)
dat=[]
for i in range(len(a)):
    dat.append((a[i],b[i],e[i]))
fig,ax=plt.subplots()
for elm in dat:
    ax.scatter(elm[1],elm[2],marker=r"$ {} $".format(elm[0]), edgecolors='none')
ax.set_ylabel('CM Designs')
ax.set_xlabel('Native')
ax.set_title('Amino Acid Frequencies')
ax.plot(x,y)
plt.savefig('CM.pdf')
plt.close()
'''
