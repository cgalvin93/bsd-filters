'''
Natural sequences of these binding sites were obtained using the protein family
 alignment from Pfam and filtered to remove all redundant sequences.


 Once the
multiple alignment is de®ned, the pro®le is constructed by counting the numbers of each amino
acid at each position along the multiple alignment.
These counts are transformed into probabilities by
normalizing the counts by the total number of
amino acids and gaps observed at that position.
These empirical probabilities re¯ect the likelihood
of observing any amino acid k at position i. Since
the counts are based on a ®nite set of sequences it
can happen that not all 20 amino acids are
observed at each position. Therefore, pseudo
counts are introduced so that no amino acid has a
zero probability to occur at position i. For more
information on pro®le generating techniques, see
Gribskov & Veretnik.26

1. get the sequences
2. function to turn sequences into probability distributions
3. function to calculate js divergence for each position in 2 profiles
pps = 1-js
is it just for positions that are designable? or does it not matter since those
only differences, so doing for whole thing captures anyway

pseudo counts you can just add one to every event
need to add to this code, change lists to counts and start at 1
then string together into script that takes 2 alignment fastas and returns pps
for each position
then im gonna need a way to generate seq alignments from list of pdb files
output by filter script. at that point i have
everything but the actual strc and native sequences to run
score->filter->subset_alignments->pps_subset/natural_binding_pos

actual strc
native sequences
position of binding site residues in native alignment
position of bionding site residues in designs alignment 
'''
#return a list of the unique sequences in a fasta file seq alignment
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


#return matrix of probability distribution for each amino acid at each position
#rows = amino acids
#columns = frequencies at each position
import numpy as np
def p_dist(seqs):
    n_positions=len(seqs[0])
    n_seqs=len(seqs)
    a=np.zeros(shape=(20,n_positions))
    ala=[];asp=[];asn =[];arg =[];cys=[];phe =[];gly=[];glu =[];gln =[];his =[];
    ile =[];leu=[];lys =[];met =[];pro =[];ser =[];trp =[];tyr =[];thr=[];val=[]
    for i in range(0,n_positions):
        nala=[];nasp=[];nasn =[];narg =[];ncys=[];nphe =[];ngly=[];nglu =[];ngln =[];nhis =[];
        nile =[];nleu=[];nlys =[];nmet =[];npro =[];nser =[];ntrp =[];ntyr =[];nthr=[];nval=[];
        ngap=[]
        for seq in seqs:
            aa=seq[i]
            if aa=='A':
                nala.append(aa)
            elif aa=='C':
                ncys.append(aa)
            elif aa=='D':
                nasp.append(aa)
            elif aa=='E':
                nglu.append(aa)
            elif aa=='F':
                nphe.append(aa)
            elif aa=='G':
                ngly.append(aa)
            elif aa=='H':
                nhis.append(aa)
            elif aa=='I':
                nile.append(aa)
            elif aa=='K':
                nlys.append(aa)
            elif aa=='L':
                nleu.append(aa)
            elif aa=='M':
                nmet.append(aa)
            elif aa=='N':
                nasn.append(aa)
            elif aa=='P':
                npro.append(aa)
            elif aa=='Q':
                ngln.append(aa)
            elif aa=='R':
                narg.append(aa)
            elif aa=='S':
                nser.append(aa)
            elif aa=='T':
                nthr.append(aa)
            elif aa=='V':
                nval.append(aa)
            elif aa=='W':
                ntrp.append(aa)
            elif aa=='Y':
                ntyr.append(aa)
            else:
                ngap.append(aa)
        ntot=float(n_seqs-len(ngap));fala=len(nala)/ntot
        fasp=len(nasp)/ntot;fasn =len(nasn)/ntot;farg =len(narg)/ntot;
        fcys=len(ncys)/ntot;fphe =len(nphe)/ntot;fgly=len(ngly)/ntot;
        fglu =len(nglu)/ntot;fgln =len(ngln)/ntot;fhis =len(nhis)/ntot;
        ffile =len(nile)/ntot;fleu=len(nleu)/ntot;flys =len(nlys)/ntot;
        fmet =len(nmet)/ntot;fpro =len(npro)/ntot;fser =len(nser)/ntot;
        ftrp =len(ntrp)/ntot;ftyr =len(ntyr)/ntot;fthr=len(nthr)/ntot;
        fval=len(nval)/ntot
        ala.append(fala);asp.append(fasp);asn.append(fasn);arg.append(farg);
        cys.append(fcys);phe.append(fphe);gly.append(fgly);glu.append(fglu);
        gln.append(fgln);his.append(fhis);ile.append(ffile);leu.append(fleu);
        lys.append(flys);met.append(fmet);pro.append(fpro);ser.append(fser);
        trp.append(ftrp);tyr.append(ftyr);thr.append(fthr);val.append(fval)
    a[0]=ala;a[1]=asp;a[2]=asn;a[3]=arg;a[4]=cys;a[5]=phe;a[6]=gly
    a[7]=glu;a[8]=gln;a[9]=his;a[10]=ile;a[11]=leu;a[12]=lys;a[13]=met
    a[14]=pro;a[15]=ser;a[16]=trp;a[17]=tyr;a[18]=thr;a[19]=val
    return a

#calculate the pps score for each position in 2 distributions
import scipy
from scipy import spatial
def calc_pps(p,q):
    n_pos = p.shape[1]
    pps_scores=[]
    for i in range(0,n_pos):
        jsd = (scipy.spatial.distance.jensenshannon(p[:,i],q[:,i]))**2
        pps=1-jsd
        pps_scores.append(pps)
    return pps_scores
