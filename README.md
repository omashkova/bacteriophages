The project was aimed at pipeline development for estimation of alpha-diversity of bacteriophages in metagenomes via hidden Markov models of single virus genes. 

Files All_pVOG_capsids.hmm, All_terminase.hmm, All_portal_hmms.hmm contain hidden Markov models of capside, terminase and portal proteins extracted from the pVOGs database (http://dmk-brain.ecn.uiowa.edu/pVOGs/) which were used within the framework of this research. Protein prediction from metagenomic assemblies were carried out using Prodigal 2.6.3, then the retrieved proteins were compared with the above mentioned hidden Markov models via hmmsearch 3.3. 

In Notebook you can find a detailed description of the analysis methods used.

File find_contigs.py contains a Python3 script calculating the number of contigs with capside, terminase and portal genes and some other statistics. 
