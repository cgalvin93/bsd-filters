#time ipython scorefile thresholds outfile

import sys

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

#read threshold conditions from threshold file
#ex. threshold file:
#delta_unsatHbonds <= 1.0
# packstat >= 0.60
# dG_separated/dSASAx100 <= 0.0
# dSASA_int >= 300.0
# sc_value >= 0.7
thresh_infile=open(sys.argv[2],'r')
threshlines=[line.strip('\n') for line in thresh_infile.readlines()]
thresh_infile.close()
conditions=[]
for i in threshlines:
    thresh_name,equality_condition,thresh_value=i.split()
    conditions.append((thresh_name,equality_condition,float(thresh_value)))

#for each condition, go through the data and store the names of structures that
#satisfy the condition
strc_sat_counter=[]
for condition in conditions:
    for x in datas:
        if x[1]==condition[0]:
            observed_val=x[2]
            if condition[1]=='>=':
                if float(observed_val) >= condition[2]:
                    strc_sat_counter.append(x[0])
            elif condition[1]=='<=':
                if float(observed_val) <= condition[2]:
                    strc_sat_counter.append(x[0])

#identify the names of structures that satisfy all conditions,
#which will have as many entries in the strc_sat_counter list
#as there are conditions in the conditions list
final_structures=[]
for structure in set(strc_sat_counter):
    count=0
    for entry in strc_sat_counter:
        if structure==entry:
            count+=1
    if count==len(conditions):
        final_structures.append(structure)

#write the names of the structures that satisfy all conditions to text file
ofile=open(sys.argv[3],'w')
ofile.write('CONDITIONS:\n')
for condition in conditions:
    s=condition[0] + ' ' + condition[1] + ' ' + str(condition[2])
    ofile.write(s+'\n')
nfinstruct=str(len(final_structures))
ofile.write('\nNUMBER OF STRUCTURES SATISFYING: ' + nfinstruct +'\n')
ofile.write('\nSTRUCTURES:\n')
for s in final_structures:
    ofile.write(s+'\n')
ofile.close()
