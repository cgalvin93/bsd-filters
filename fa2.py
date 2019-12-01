#time ipython ~/desktop/bsd_filters/fa2.py

import Bio
from Bio import SeqUtils
from Bio.SeqUtils import IUPACData
import sys
import numpy as np
import scipy
from scipy import spatial
import pandas as pd

#terms in the score file that are more favorable when value higher must be
#specified:
higher_better=['packstat','dSASA_hphobic','dSASA_int','dSASA_polar',
               'hbonds_int','nres_int','acc']

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
def return_bs_seq(pdbfiles,method,pdbid):
    seqs=[]
    for file in pdbfiles:
        fullseq=[]
        path=method+'/'+pdbid+'/'+file
        f=open(path,'r');lines=[line for line in f.readlines()]
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
def return_topq_strc(term,low,high,datas):
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


def return_toscan(file):
    scorefile=open(file,'r')
    lines=[line for line in scorefile.readlines()[1:]]
    scorefile.close()
    terms=[]
    for i in lines[0].split()[1:]:
        terms.append(i)
    del lines[0]
    strc_names=[]
    datas=[]
    for line in lines:
        properties=line.split()[1:]
        name=properties[-1]
        strc_names.append(name)
        for index,term in enumerate(terms[:-1]):
            val=properties[index]
            datas.append((name,term,val))
    to_scan=[]
    for term in terms:
        if term != 'description':
            vals=[]
            for entry in datas:
                if entry[1]==term:
                    vals.append(float(entry[2]))
            high=max(vals); low=min(vals)
            to_scan.append((term,low,high))
    return datas,to_scan


def return_top_ts(data,scan_list):
    ts_top=[]
    for term,low,high in scan_list:
        if term=='total_score':
            ts_top=return_topq_strc(term,low,high,data)
        else:
            pass
    return ts_top


#must put binding site res numbers here manually
#store pdb numbering of binding site residues for query protein
ptn_ids=['2xbn','1f4p','1zk4','3dk9','3r2q'] #3dlc,3s6f didnt work to get all scores
bre_all=[]
cm_all=[]
fd_all=[]
#ptn_ids=['2xbn','1f4p'] #for testing
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
    p_nat=p_dist(native_seqs)
    bre_sc='bre/'+ptn_pdb_id+'/bre-'+ptn_pdb_id+'-all-score.sc'
    cm_sc='cm/'+ptn_pdb_id+'/cm-'+ptn_pdb_id+'-all-score.sc'
    fd_sc='fd/'+ptn_pdb_id+'/fd-'+ptn_pdb_id+'-all-score.sc'
    bredatas,brescan=return_toscan(bre_sc);bre_ts=return_top_ts(bredatas,brescan)
    cmdatas,cmscan=return_toscan(cm_sc);cm_ts=return_top_ts(cmdatas,cmscan)
    fddatas,fdscan=return_toscan(fd_sc);fd_ts=return_top_ts(fddatas,fdscan)
    bre_ts_seqs=return_bs_seq(bre_ts,'bre',ptn_pdb_id);bre_ts_dist=p_dist(bre_ts_seqs)
    cm_ts_seqs=return_bs_seq(cm_ts,'cm',ptn_pdb_id);cm_ts_dist=p_dist(cm_ts_seqs)
    fd_ts_seqs=return_bs_seq(fd_ts,'fd',ptn_pdb_id);fd_ts_dist=p_dist(fd_ts_seqs)
    bre_ts_pps=calc_pps(p_nat,bre_ts_dist);bre_ts_med_pps=np.median(bre_ts_pps)
    cm_ts_pps=calc_pps(p_nat,cm_ts_dist);cm_ts_med_pps=np.median(cm_ts_pps)
    fd_ts_pps=calc_pps(p_nat,fd_ts_dist);fd_ts_med_pps=np.median(fd_ts_pps)
    bre_results=[] #strc,term,medpps,percd
    cm_results=[]
    fd_results=[]
    for term,low,high in brescan:
        if term!='total_score'and term!='description':
            top_strc=return_topq_strc(term,low,high,bredatas)
            top_strc_seqs=return_bs_seq(top_strc,'bre',ptn_pdb_id)
            top_strc_dist=p_dist(top_strc_seqs)
            top_strc_pps=calc_pps(p_nat,top_strc_dist)
            top_strc_med_pps=np.median(top_strc_pps)
            percd=((top_strc_med_pps-bre_ts_med_pps)/bre_ts_med_pps)*100
            bre_results.append((ptn_pdb_id,term,top_strc_med_pps,percd))
    for term,low,high in cmscan:
        if term!='total_score'and term!='description':
            top_strc=return_topq_strc(term,low,high,cmdatas)
            top_strc_seqs=return_bs_seq(top_strc,'cm',ptn_pdb_id)
            top_strc_dist=p_dist(top_strc_seqs)
            top_strc_pps=calc_pps(p_nat,top_strc_dist)
            top_strc_med_pps=np.median(top_strc_pps)
            percd=((top_strc_med_pps-cm_ts_med_pps)/cm_ts_med_pps)*100
            cm_results.append((ptn_pdb_id,term,top_strc_med_pps,percd))
    for term,low,high in fdscan:
        if term!='total_score'and term!='description':
            top_strc=return_topq_strc(term,low,high,fddatas)
            top_strc_seqs=return_bs_seq(top_strc,'fd',ptn_pdb_id)
            top_strc_dist=p_dist(top_strc_seqs)
            top_strc_pps=calc_pps(p_nat,top_strc_dist)
            top_strc_med_pps=np.median(top_strc_pps)
            percd=((top_strc_med_pps-fd_ts_med_pps)/fd_ts_med_pps)*100
            fd_results.append((ptn_pdb_id,term,top_strc_med_pps,percd))
    bre_all.append(bre_results)
    cm_all.append(cm_results)
    fd_all.append(fd_results)


def result_table(datalist):
    df=pd.DataFrame()
    termnames=[]
    for i in datalist[0]:
        termnames.append(i[1])
    df['metric']=pd.Series(termnames)
    for i in datalist:
        strc_name=i[0][0]
        vals_each_strc=[]
        for strc,term,medpps,percd in i:
            vals_each_strc.append(percd)
        df[strc_name]=pd.Series(vals_each_strc)
    return df

bre_df=result_table(bre_all)
bre_df.to_csv('bre_percd.csv')
cm_df=result_table(cm_all)
cm_df.to_csv('cm_percd.csv')
fd_df=result_table(fd_all)
fd_df.to_csv('fd_percd.csv')

'''
not part of script, rather i am running this in terminal using csv files
generated by script

bre_df=pd.read_csv('bre_percd.csv')
bre_percd_means=[]
for i in range(0,bre_df.shape[0]): #n_rows
    rowvals=[]; term=bre_df.iloc[i][1]
    for x in bre_df.iloc[i][2:7]: #2:7 skips row nums, metric names
        rowvals.append(x)
    mean=np.mean(rowvals);std=np.std(rowvals)
    bre_percd_means.append((term,mean,std))

bre_percd_means.sort(key=lambda means : means[1]) #sort by mean
bre_percd_means.reverse() #reverse so highest means first

bre_percd_std=[]
for i in range(0,bre_df.shape[0]): #n_rows
    rowvals=[]; term=bre_df.iloc[i][1]
    for x in bre_df.iloc[i][2:7]: #2:7 skips row nums, metric names
        rowvals.append(x)
    mean=np.mean(rowvals);std=np.std(rowvals)
    bre_percd_std.append((term,mean,std))
bre_percd_std.sort(key=lambda std : std[2]) #sort by lowest std

df_bre_means=pd.DataFrame()
def make_df(list):
    df=pd.DataFrame()
    names=[]
    means=[]
    stds=[]
    for a,b,c in list:
        names.append(a)
        means.append(b)
        stds.append(c)
    df['metric']=pd.Series(names)
    df['mean %diff']=pd.Series(means)
    df['std']=pd.Series(stds)
    return df
bremeandf=make_df(bre_percd_means)
bremeandf.to_csv('bremeans.csv')

cm_df=pd.read_csv('cm_percd.csv')
cm_percd_means=[]
for i in range(0,cm_df.shape[0]): #n_rows
    rowvals=[]; term=cm_df.iloc[i][1]
    for x in cm_df.iloc[i][2:7]: #2:7 skips row nums, metric names
        rowvals.append(x)
    mean=np.mean(rowvals);std=np.std(rowvals)
    cm_percd_means.append((term,mean,std))

cm_percd_means.sort(key=lambda means : means[1]) #sort by mean
cm_percd_means.reverse() #reverse so highest means first

cmmeandf=make_df(cm_percd_means)
cmmeandf.to_csv('cmmeans.csv')

fd_df=pd.read_csv('fd_percd.csv')
fd_percd_means=[]
for i in range(0,fd_df.shape[0]): #n_rows
    rowvals=[]; term=fd_df.iloc[i][1]
    for x in fd_df.iloc[i][2:7]: #2:7 skips row nums, metric names
        rowvals.append(x)
    mean=np.mean(rowvals);std=np.std(rowvals)
    fd_percd_means.append((term,mean,std))

fd_percd_means.sort(key=lambda means : means[1]) #sort by mean
fd_percd_means.reverse() #reverse so highest means first

fdmeandf=make_df(fd_percd_means)
fdmeandf.to_csv('fdmeans.csv')

#finally, to take the mean across all methods
fd_df=pd.read_csv('fd_percd.csv') #..... for bre and cm also
all_means=[]
for i in range(0,fd_df.shape[0]): #n_rows
    rowvals=[]; term=fd_df.iloc[i][1]
    for x in fd_df.iloc[i][2:7]: #2:7 skips row nums, metric names
        rowvals.append(x)
    for x in bre_df.iloc[i][2:7]: #2:7 skips row nums, metric names
        rowvals.append(x)
    for x in cm_df.iloc[i][2:7]: #2:7 skips row nums, metric names
        rowvals.append(x)
    mean=np.mean(rowvals);std=np.std(rowvals)
    all_means.append((term,mean,std))

all_means.sort(key=lambda means : means[1])
all_means.reverse()

all_df=make_df(all_means)
all_df.to_csv('allmeans.csv')

'''

#take product of mean and ratio mean/std ??? rank by that ??? (to balance both)
#next will be to take these tables and get the mean and var
#across each row (ie for each term)
#to get table thats
#x mean, var
#y term names
#finally, I will just need to repeat this process using all values for the
#3 methods

#Dataframe.[ ] ; str of column title; + can add int of row (['column'][0])
# Dataframe.loc[ ] : This function is used for row labels. .loc[0] = 1st row
#columns come after ie Dataframe.loc[["row1", "row2"], ["column1", "column2", "column3"]]
#to select all rows and some columns:
# Dataframe.loc[:, ["column1", "column2", "column3"]]
# Dataframe.iloc[ ] : This function is used for positions or integer based
#The df.iloc indexer is very similar to df.loc but only uses integer locations to make its selections.
#uses row, column
#data.iloc [[3, 4], [1, 2]]
# Dataframe.ix[] : This function is used for both label and integer based

#pareto front of relevant term scores
#pps of non dominated structures
