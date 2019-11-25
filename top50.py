#find the top quintile for each of the relevant metrics, minimum 40 or add 2nd q
#amongst the total set of structures generated (those belonging to the top
#50 of at least one metric), report the number of sets that each protein
#belongs to, along with what those sets are.
#Take the top 50 strc by most sets and find pps.

#python filter_analysis.py scorefile native_seqs outfile pdbID
#must be run in directory with design pdb files so it can get sequences
#to calculate pps with

#time ipython ~/desktop/bsd_filters/multi_filter.py bre-1f4p-all-score.sc ~/desktop/prjk/analysis/native_seqs/1f4p_binding_site.txt 1f4p_top_fifty.pdf 1f4p


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
import pandas as pd

#must put binding site res numbers here manually
#store pdb numbering of binding site residues for query protein
ptn_pdb_id=sys.argv[4]
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

#terms in the score file that are more favorable when value higher must be
#specified:
higher_better=['packstat','dSASA_hphobic','dSASA_int','dSASA_polar',
               'hbonds_int','nres_int','acc']

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


#find the lowest and highest values for relevant terms,
#store in a list (term, low, high) that will be used for calculations
#need to look at plots and edit this list
# acc
# buns_bb_heavy
# cav
# dSASA_int
# exphyd
#dgsep/dsasa
# fa_intra_sol_xover4
# lk_ball_wtd
#omega
# p_aa_pp
# packstat
# rama_prepro
# ref
# fa_atr
#using 14 terms
to_scan=[]
for term in terms:
    if term == 'acc' or term == 'buns_bb_heavy' or term=='cav' or term == 'dG_separated/dSASAx100' or term == 'dSASA_int' or term == 'exphyd' or term=='fa_intra_sol_xover4' or term == 'lk_ball_wtd' or term == 'omega' or term=='p_aa_pp' or term=='packstat' or term=='rama_prepro' or term=='ref' or term=='fa_atr':
        vals=[]
        for entry in datas:
            if entry[1]==term:
                vals.append(float(entry[2]))
        high=max(vals); low=min(vals)
        to_scan.append((term,low,high))
    else:
        pass

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
        strs.append(s.upper())
    return list(set(strs))


#calculate the pps score for each position in 2 distributions
def calc_pps(p,q):
    n_pos = p.shape[1]
    pps_scores=[]
    for i in range(0,n_pos):
        jsd = (scipy.spatial.distance.jensenshannon(p[:,i],q[:,i]))**2
        pps=1-jsd
        pps_scores.append(pps)
    return pps_scores

#return the structures in top quintile of specified term
#if less than 35 it includes second quintile as well
def return_topq_strc(term,low,high):
    vals = foursteps(low,high)
    strc=[]
    if term in higher_better: #gotta do terms that obs > thresh and obs < thresh sep
        conditions=[]
        conditions.append((term,low)) #low for > gives all strc
        for i in vals:
            conditions.append((term,i))
        firstqstrc = return_filtered_strc(datas, conditions[4], '>=')
        for i in firstqstrc:
            strc.append(i)
        if len(firstqstrc) < 35:
            secondqstrc = return_filtered_strc(datas, conditions[3], '>=')
            for i in secondqstrc:
                strc.append(i)
    else:
        conditions=[]
        for i in vals:
            conditions.append((term,i))
        conditions.append((term, high)) #< high gives all strc
        firstqstrc = return_filtered_strc(datas, conditions[0], '<=')
        for i in firstqstrc:
            strc.append(i)
        if len(firstqstrc) < 35:
            secondqstrc = return_filtered_strc(datas, conditions[1], '<=')
            for i in secondqstrc:
                strc.append(i)
    return strc

#get a list of all unique strc in top quintiles and
#pandas table with strc in topq of each term
all_strc = []
df = pd.DataFrame()
for term,low,high in to_scan:
    top_strc = return_topq_strc(term,low,high)
    for i in top_strc:
        all_strc.append(i)
    df[term]=pd.Series(top_strc)
all_strc = list(set(all_strc))

#get the number of sets that each strc is in
n_sets_each_strc=[]
for structure in all_strc:
    count=0
    for i in to_scan:
        if structure in list(df[i[0]]):
            count+=1
    n_sets_each_strc.append((structure,count))

#sort the structures by number of sets theyre in and get the top 50
n_sets_each_strc.sort(key=lambda entry:entry[1])
n_sets_each_strc.reverse()
top_fifty=[i[0] for i in n_sets_each_strc[0:50]]

#open the output pdf file to put the graphs on
outfilename=sys.argv[3]
pdf = PdfPages(outfilename)
#ofile=open(outfilename[:-3]+'txt')
#find the pps scores for the top 50 and top 25 structures
native_seqs=capture_seqs(sys.argv[2])
p_nat=p_dist(native_seqs)
mut_seqs = return_bs_seq(top_fifty)
p_mut=p_dist(mut_seqs)
top_twofive_seqs=return_bs_seq(top_fifty[0:25])
p_twofive=p_dist(top_twofive_seqs)
pps_scores=calc_pps(p_nat,p_mut)
twofive_pps=calc_pps(p_nat,p_twofive)
#get the pps scores for all designs to compare against
all_designs=[]
for i in strc_names:
    all_designs.append(i[:-5]+'.pdb')
whole_set_seqs=return_bs_seq(all_designs)
p_all_mut=p_dist(whole_set_seqs)
ref_pps=calc_pps(p_nat,p_all_mut)
#get the pps scores for all strc in at least one top quintile group
top_q_strc_filenames=[]
for i in all_strc:
    top_q_strc_filenames.append(i[:-4]+'.pdb') #why -4 here idk ???????
topq_seqs=return_bs_seq(top_q_strc_filenames)
p_topq=p_dist(topq_seqs)
topq_pps=calc_pps(p_nat,p_topq)

#plot the pps distributions for the top 50 and topq structures
fig,ax=plt.subplots()
both_pps=[]
both_pps.append(ref_pps)
both_pps.append(topq_pps)
both_pps.append(pps_scores)
both_pps.append(twofive_pps)
bp_dict = ax.boxplot(both_pps)
ax.set_ylabel('PPS')
xlocs=[x+1 for x in range(len(both_pps))]
xlabs=['all','topq_strc','top 50','top 25']
plt.xticks(xlocs, xlabs)
for line in bp_dict['medians']:
    x,y = line.get_xydata()[1]
    plt.text(x,y, '%.3f' % y,verticalalignment='center')
pdf.savefig()
plt.clf()

#plot showing the pps distributions for the topq structures of each term
pps_each_term=[]
ordered_terms=[]
for i in to_scan:
    term=i[0]
    ordered_terms.append(term)
    term_strc=[]
    for strc in list(df[term]):
        if type(strc)!=float:
            term_strc.append(strc)
    term_seqs = return_bs_seq(term_strc)
    p_term=p_dist(term_seqs)
    term_pps=calc_pps(p_nat,p_term)
    pps_each_term.append(term_pps)
fig,ax=plt.subplots()
bp_dict = ax.boxplot(pps_each_term)
ax.set_ylabel('PPS')
xlocs=[x+1 for x in range(len(ordered_terms))]
plt.xticks(xlocs, ordered_terms, rotation='vertical')
for line in bp_dict['medians']:
    x,y = line.get_xydata()[1]
    plt.text(x,y, '%.3f' % y,verticalalignment='center')
pdf.savefig()
plt.clf()


#ofile.close()
pdf.close()
