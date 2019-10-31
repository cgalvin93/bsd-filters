#python filter_analysis.py scorefile native_seqs outfile
#must be run in directory with design pdb files so it can get sequences
#to calculate pps with

#1. pdblist.py
#2. feed list to scoring xml
#3. feed scorefile here

#time python ~/desktop/bsd_filters/filter_analysis.py 3r2q-iam-score.sc ~/desktop/prjk/analysis/native_seqs/3r2q_binding_site.txt fd-3r2q-iam.pdf
#time python ~/desktop/bsd_filters/filter_analysis.py 1f4p-iam-score.sc ~/desktop/prjk/analysis/native_seqs/1f4p_binding_site.txt fd-1f4p-iam.pdf
#time python ~/desktop/bsd_filters/filter_analysis.py 1zk4-iam-score.sc ~/desktop/prjk/analysis/native_seqs/1zk4_binding_site.txt fd-1zk4-iam.pdf



#must put binding site res numbers here manually
#pdb numbering
#2xbn:
# binding_site_res=[73, 133, 134, 135, 138, 159,161, 202, 231, 233,
#                   234, 262,264, 265]
#1f4p:
# binding_site_res=[10, 11, 12, 13, 14, 15, 58, 59,
#                   60, 61, 62, 68, 93, 94, 95, 98,
#                   100, 101, 102]
#1zk4:
binding_site_res=[13, 14, 15, 16, 17, 18, 19, 36,
                  37, 38, 62, 63, 89, 90, 91, 92,
                  112, 140, 142, 155, 159]
#3dk9:
# binding_site_res=[26, 27, 28, 29, 30, 31, 49, 50,
#                   51, 52, 56, 57, 62,
#                   66, 129, 130, 155, 156, 157,
#                   177, 181, 197, 198, 201, 202,
#                   291, 294, 298, 330, 331]
#3dlc:
# binding_site_res=[48, 49, 50,
#                   51, 52, 53, 55, 72, 73, 74, 77,
#                   100, 101, 102, 117, 118, 119,
#                   122, 123]
#3r2q:
# binding_site_res=[9, 10, 11, 33, 34, 48, 49, 50,
#                   62, 63, 64]

#terms in the score file that are more favorable when value higher must be
#specified:
higher_better=['packstat','dSASA_hphobic','dSASA_int','dSASA_polar',
               'hbonds_int','nres_int']

import sys
import Bio
from Bio import SeqUtils
from Bio.SeqUtils import IUPACData
import sys
import numpy as np
import scipy
from scipy import spatial
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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


#find the lowest and highest values for each term,
#store in a list (term, low, high) that I can then feed to
#'scan feature' function that will do five intervals between low and hi
to_scan=[]
for term in terms:
    if term != 'description':
        vals=[]
        for entry in datas:
            if entry[1]==term:
                vals.append(float(entry[2]))
        high=max(vals); low=min(vals)
        to_scan.append((term,low,high))

#find four equally spaced values between high and low value
def foursteps(min, max):
    d = max - min
    int_size = d/5
    p1=min+int_size;p2=min+2*int_size;p3=min+3*int_size
    p4=min+4*int_size
    threshvals=[p1,p2,p3,p4]
    return threshvals

#based on a condition and eq (< or >) return strc from score data that
#satisfy the condition
def return_filtered_strc(datas, condition, eq):
    strc=[]
    for x in datas:
        if x[1]==condition[0]:
            observed_val=x[2]
            if eq=='>=':
                if float(observed_val) >= condition[1]:
                    strc.append(x[0][:-5]+'.pdb')  #formatting pdbfile names
            elif eq=='<=':
                if float(observed_val) <= condition[1]:
                    strc.append(x[0][:-5]+'.pdb')
    return strc


#convert 3 letter aa code to 1 letter
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


#return matrix of probability distribution for each amino acid at each position
#uses pseudo counts to avoid any P=0
#rows = amino acids
#columns = frequencies at each position
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
        strs.append(s)
    return list(set(strs))

#get sequences and frequency matrix of native binding site sequences
native_seqs=capture_seqs(sys.argv[2])
p_nat=p_dist(native_seqs)

#calculate the pps score for each position in 2 distributions
def calc_pps(p,q):
    n_pos = p.shape[1]
    pps_scores=[]
    for i in range(0,n_pos):
        jsd = (scipy.spatial.distance.jensenshannon(p[:,i],q[:,i]))**2
        pps=1-jsd
        pps_scores.append(pps)
    return pps_scores

#open the output pdf file to put all the graphs on
outfilename=sys.argv[3]
pdf = PdfPages(outfilename)


#for every term, create 5 conditions within low/high range
#for each condition, iterate over datas and store strc that meet
#for strc that meet, get their bs seq from names
#for set of bs seq , eval pps and store
#plot the 5 pps distributions for the term, along w n strc, name, thresh vals
for term, low, high in to_scan:
    vals = foursteps(low,high)
    term_scores=[]  #will store pps scores for strc meeting 5 thresholds
    nstrcs=[]   #number of structures meeting each threshold
    if term in higher_better: #gotta do terms that obs > thresh and obs < thresh sep
        conditions=[]
        conditions.append((term,low)) #low for > gives all strc
        for i in vals:
            conditions.append((term,i))
        for condition in conditions:
            strc = return_filtered_strc(datas, condition, '>=')
            nstrcs.append(str(len(strc)))
            seqs = return_bs_seq(strc)
            p_mut=p_dist(seqs)
            pps_scores=calc_pps(p_nat,p_mut)
            term_scores.append(pps_scores)
    else:
        conditions=[]
        for i in vals:
            conditions.append((term,i))
        conditions.append((term, high)) #< high gives all strc
        conditions.reverse()
        for condition in conditions:
            strc = return_filtered_strc(datas, condition, '<=')
            nstrcs.append(str(len(strc)))
            seqs = return_bs_seq(strc)
            p_mut=p_dist(seqs)
            pps_scores=calc_pps(p_nat,p_mut)
            term_scores.append(pps_scores)
    fig,ax=plt.subplots()       #plotting the five pps distributions
    bp_dict = ax.boxplot(term_scores)
    ax.set_title(term)
    ax.set_ylabel('PPS')
    xlocs=[x+1 for x in range(len(term_scores))]
    xlabs=[str(b)[0:7] for a,b in conditions]
    plt.xticks(xlocs, xlabs)
    for line in bp_dict['medians']:
        x,y = line.get_xydata()[1]
        plt.text(x,y, '%.3f' % y,verticalalignment='center')
    for i, line in enumerate(bp_dict['boxes']):
        x, y = line.get_xydata()[0]
        plt.text(x,y, nstrcs[i],
             horizontalalignment='center',
             verticalalignment='top')
    pdf.savefig()
    plt.clf()

pdf.close()
