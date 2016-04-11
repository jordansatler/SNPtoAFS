#!/usr/bin/env python
"""This script will convert a SNP output file from pyRAD to an input 
SNP file ready for use in AFS_FSC_total.py for conversion to 
an observed AFS. 

author: J. Satler
date: 18 Jan 2016

python SNPtoAFSready.py file(s)
"""

import sys

#dictionary of ambiguity codes
ambig = {"Y":["C", "T"], "R":["A", "G"], "W":["A", "T"], "S":["C", "G"],
         "K":["G","T"], "M":["A", "C"]}

def duplicate(file):
    """duplicate all individuals so they are diploid"""
    with open(file, 'r') as Infile:
        d = []
        
        for line in Infile:
            if not line.startswith("#"):
                d.append(line)
                d.append(line)
        return d
        
def header(SNP_list):
    """create a header for each SNP to record linkage patterns"""
    getL = SNP_list[0]
    getL = getL.strip().split()
    
    locus = 1
    L = [len(i) for i in getL[1:] if i != "_"]
    
    head = ['']
    for j in L:
        name = "LocusNumber" + str(locus)
        #write N number of times to list
        for reps in range(j):
             head.append(name)
        locus += 1
    return head

def phase(SNP_dup):
    """resolve SNPs that contain ambiguity codes"""
    count = 2
    
    phased = []
    for line in SNP_dup:
        line = line.strip().split()

        if count % 2 == 0:
            seq = [line[0] + "_a"]
        else:
            seq = [line[0] + "_b"]
        
        for bp in line[1:]:

            for nc in bp:
                
                if nc in ambig:
            
                    if count % 2 == 0:
                        seq.append(ambig[nc][0])
                        
                    else:
                        seq.append(ambig[nc][1])
                        
                #Skip loci that are invariant
                elif nc == "_":
                    pass
                
                else:
                    seq.append(nc)
        phased.append(seq)
        count += 1
    return phased

            
def write_out(name, header, SNPs):
    """write the new SNP file"""
    
    with open(name[:-5] + "_AFSready.txt", "w") as out:
        for i in range(len(header)):
            if i != len(header) - 1:
                out.write(header[i] + "\t")
            else:
                out.write(header[i] + "\n")
        
        for j in range(len(SNPs)):
            for k in range(len(SNPs[j])):
                if k != len(SNPs[j]) - 1:
                    out.write(SNPs[j][k] + "\t")
                else:
                    out.write(SNPs[j][k] + "\n")

def Main():
    dip = duplicate(sys.argv[1])
    headL = header(dip)
    ph = phase(dip)
    write_out(sys.argv[1], headL, ph)

if __name__ == "__main__":
    Main()