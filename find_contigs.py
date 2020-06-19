import numpy as np

# This function returns an array of contigs containing at least one portal gene. 
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
    
# This function returns the number of contigs containing at least one portal gene. 
# Input: a file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
def get_portal_number(filename):
    p = get_portal(filename)
    return len(p)
    
# This function returns an array of contigs containing at least one terminase gene. 
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
    
# This function returns the number of contigs containing at least one terminase gene. 
# Input: a file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
def get_terminase_number(filename):
    t = get_terminase(filename)
    return len(t)
    
# This function returns an array of contigs containing at least one capsid gene. 
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
    
# This function returns the number of contigs containing at least one capsid gene. 
# Input: a file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def get_capsid_number(filename):
    c = get_capsid(filename)
    return len(c)

# This function returns an array of contigs containing only terminase genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def terminase_only(port_filename, term_filename, cap_filename):
    t = get_terminase(term_filename)
    p = get_portal(port_filename)
    c = get_capsid(cap_filename)
    t_only = np.asarray([name for name in t if name not in c and name not in p])
    return t_only

# This function returns the number of contigs containing only terminase genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def terminase_only_number(port_filename, term_filename, cap_filename):
    t_only = terminase_only(port_filename, term_filename, cap_filename)
    return len(t_only)

# This function returns an array of contigs containing only capsid genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def capsid_only(port_filename, term_filename, cap_filename):
    t = get_terminase(term_filename)
    p = get_portal(port_filename)
    c = get_capsid(cap_filename)
    c_only = np.asarray([name for name in c if name not in t and name not in p])
    return c_only

# This function returns the number of contigs containing only capsid genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def capsid_only_number(port_filename, term_filename, cap_filename):
    c_only = capsid_only(port_filename, term_filename, cap_filename)
    return len(c_only)

# This function returns an array of contigs containing only portal genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def portal_only(port_filename, term_filename, cap_filename):
    t = get_terminase(term_filename)
    p = get_portal(port_filename)
    c = get_capsids(cap_filename)
    p_only = np.asarray([name for name in p if name not in t and name not in c])
    return p_only

# This function returns the number of contigs containing only portal genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def portal_only_number(port_filename, term_filename, cap_filename):
    p_only = portal_only(port_filename, term_filename, cap_filename)
    return len(p_only)

# This function returns an array of contigs containing only terminase and capsid genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def terminase_and_capsid_only(port_filename, term_filename, cap_filename):
    t = get_terminase(term_filename)
    p = get_portal(port_filename)
    c = get_capsid(cap_filename)
    t_and_c_only = np.asarray([name for name in t if name in c and name not in p])
    return t_and_c_only

# This function returns the number of contigs containing only terminase and capsid genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def terminase_and_capsid_only_number(port_filename, term_filename, cap_filename):
    t_and_c_only = terminase_and_capsid_only(port_filename, term_filename, cap_filename)
    return len(t_and_c_only)

# This function returns an array of contigs containing only terminase and portal genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def terminase_and_portal_only(port_filename, term_filename, cap_filename):
    t = get_terminase(term_filename)
    p = get_portal(port_filename)
    c = get_capsid(cap_filename)
    t_and_p_only = np.asarray([name for name in t if name in p and name not in c])
    return t_and_p_only

# This function returns the number of contigs containing only terminase and portal genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def terminase_and_portal_only_number(port_filename, term_filename, cap_filename):
    t_and_p_only = terminase_and_portal_only(port_filename, term_filename, cap_filename)
    return len(t_and_p_only)

# This function returns an array of contigs containing only capsid and portal genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def portal_and_capsid_only(port_filename, term_filename, cap_filename):
    t = get_terminase(term_filename)
    p = get_portal(port_filename)
    c = get_capsid(cap_filename)
    p_and_c_only = np.asarray([name for name in c if name in p and name not in t])
    return p_and_c_only

# This function returns the number of contigs containing only capsid and portal genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def portal_and_capsid_only_number(port_filename, term_filename, cap_filename):
    p_and_c_only = portal_and_capsid_only(port_filename, term_filename, cap_filename)
    return len(p_and_c_only)

# This function returns an array of contigs containing capsid, terminase and portal genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def portal_and_terminase_and_capsid(port_filename, term_filename, cap_filename):
    t = get_terminase(term_filename)
    p = get_portal(port_filename)
    c = get_capsid(cap_filename)
    p_and_t_and_c = np.asarray([name for name in t if name in c and name in p])
    return p_and_t_and_c

# This function returns the number of contigs containing capsid, terminase and portal genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def portal_and_terminase_and_capsid_number(port_filename, term_filename, cap_filename):
    p_and_t_and_c = portal_and_terminase_and_capsid(port_filename, term_filename, cap_filename)
    return len(p_and_t_and_c)

# This function returns the average number of contigs containing capsid, terminase or portal genes. 
# Input (in this order): 
# - file target.out extracted from hmmsearch for hidden Markov models of portal genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of terminase genes (see Notebook)
# - file target.out extracted from hmmsearch for hidden Markov models of capsid genes (see Notebook)
def average_capsid_terminase_portal_number(port_filename, term_filename, cap_filename):
    t = get_terminase(term_filename)
    p = get_portal(port_filename)
    c = get_capsid(cap_filename) 
    return (len(p) + len(t) + len(c)) / 3
