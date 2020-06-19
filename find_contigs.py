import numpy as np

# This function returns an array of contigs containing at least one portal gene protein. 
# Input: a file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
def get_portal(filename): 
    portals = open(filename, 'r')
    portal_res = []
    for line in portals:
        if line.startswith('#'): # ignore subordinate lines
            continue
        line = line.split()
        portal_res.append(line[0]) # memorize the name of contig with protein number
    portals.close()
    n = len(portal_res)
    for i in range(n): # to get the pure contig name we have to remove protein numbers
        elem = portal_res[i].split('_')
        elem = elem[:-1]
        elem = '_'.join(elem)
        portal_res[i] = elem
    portal_res = np.asarray(portal_res, dtype=str)
    portal_res_unique = np.unique(portal_res) # here we remove the possibly repeated contig names 
    return portal_res_unique
    
# This function returns the number of contigs containing at least one portal gene protein. 
# Input: a file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
def get_portal_number(filename):
    p = get_portal(filename)
    return len(p)
    
# This function returns an array of contigs containing at least one terminase gene protein. 
# Input: a file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
def get_terminase(filename):
    terminases = open(filename, 'r')
    term_res = []
    for line in terminases:
        if line.startswith('#'): # ignore subordinate lines
            continue
        line = line.split()
        term_res.append(line[0]) # memorize the name of contig with protein number
    terminases.close()
    n = len(term_res)
    for i in range(n): # to get the pure contig name we have to remove protein numbers
        elem = term_res[i].split('_')
        elem = elem[:-1]
        elem = '_'.join(elem)
        term_res[i] = elem
    term_res = np.asarray(term_res, dtype=str)
    term_res_unique = np.unique(term_res) # here we remove the possibly repeated contig names
    return term_res_unique
    
# This function returns the number of contigs containing at least one terminase gene protein. 
# Input: a file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
def get_terminase_number(filename):
    t = get_terminase(filename)
    return len(t)
    
# This function returns an array of contigs containing at least one capsid gene protein. 
# Input: a file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def get_capsid(filename):
    capsids = open(filename, 'r')
    cap_res = []
    for line in capsids:
        if line.startswith('#'): # ignore subordinate lines
            continue
        line = line.split()
        cap_res.append(line[0]) # memorize the name of contig with protein number
    capsids.close()
    n = len(cap_res)
    for i in range(n): # to get the pure contig name we have to remove protein numbers
        elem = cap_res[i].split('_')
        elem = elem[:-1]
        elem = '_'.join(elem)
        cap_res[i] = elem
    cap_res = np.asarray(cap_res, dtype=str)
    cap_res_unique = np.unique(cap_res) # here we remove the possibly repeated contig names
    return cap_res_unique
    
# This function returns the number of contigs containing at least one capsid gene protein. 
# Input: a file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def get_capsid_number(filename):
    c = get_capsid(filename)
    return len(c)
