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
