#!/usr/bin/env python

#THIS SCRIPT NEEDS TO BE DOCUMENTED AND TESTED, HEAVILY!!!

"""
This script will build an allele frequency spectrum from an input 
matrix of SNPs. For use in a two population model. SNP matrix is 
filtered to all bi-allelic SNPs, that equal or surpass the population 
thresholds. Options allow for the user to use linked SNPs or to 
subsample a single SNP per locus, and if a threshold that requires
subsampling is used, the user can replicate the observed AFS N times 
due to the subsampling of alleles per SNP. Output is an observed allele 
frequency spectrum for use in fastsimcoal2.

For use with both pyRAD and AftrRAD output.
                                                     
python AFS_FSC_modular_Replicate_linkage.py SNP_infile.txt traits.txt 
       Monomorphics Threshold Replicate Linked/Unlinked

Author: Jordan Satler
Date: 25 Jan 2016
"""

import sys
import os
import shutil
import csv
import random


def pop_association(Traits):
    """Sets the individuals to their respective populations. 
       Also returns the sample counts per population."""
    
    with open(Traits, 'r') as traits:
        
        Pops = {}
        Pop_counts = {}
        next(traits)
        for line in traits:
            line = line.strip().split()
            Pops[line[0]] = line[1]
            
            if line[1] in Pop_counts:
                Pop_counts[line[1]] += 1
            else:
                Pop_counts[line[1]] = 1
        return Pops, Pop_counts

def build_AFS_LoL(file, Pop_counts, Threshold):
    """Build an empty allele frequency spectrum. Return population
       names and sample numbers. Uses lists of lists."""
    
    Populations = Pop_counts.items()
    
    pop1 = Populations[0]
    pop2 = Populations[1]
    
    t = float(Threshold) / 100
    
    pop1_thr = int(pop1[1] * t)
    pop2_thr = int(pop2[1] * t)
    
    AFS = [[0.0] * (pop1_thr + 1) for _ in range(pop2_thr + 1)]
    return AFS, pop1, pop2, pop1_thr, pop2_thr 


def Biallelic_SNPs(Infile, Pops, pop1, pop2, pop1_thr, pop2_thr):
    """Filter SNP matrix and retain only bi-allelic SNPs that are
       equal to or above both population thresholds."""
       
    RawSNPs = csv.reader(open(Infile), delimiter = '\t')
    Data = zip(*RawSNPs)

    indivs = Data[0]
    columns = Data[1:]
    
    #List of the unique polymorphic loci.
    PolyLoci = []

    #Alleles allowed.
    allowed = ['A', 'C', 'G', 'T']

    #List of bi-allelic SNPs that meet both population thresholds.
    Bi_Thr = []
    
    for i in columns:
    
        #Create list of polymorphic loci.
        if i[0] not in PolyLoci:
            PolyLoci.append(i[0])
    
        Alleles_set = []
        Pop1_thresh = 0
        Pop2_thresh = 0

        for j in range(len(i)):
            if j > 0 and i[j] in allowed:
                if i[j] not in Alleles_set:
                    Alleles_set.append(i[j])
                
                if Pops[indivs[j]] == pop1[0]:
                    Pop1_thresh += 1 
                else:
                    Pop2_thresh += 1

        if len(Alleles_set) == 2:
            if Pop1_thresh >= pop1_thr and Pop2_thresh >= pop2_thr:
                Bi_Thr.append(i)
    return Bi_Thr, indivs, len(columns), len(PolyLoci)

def subsample(Bi_Thr):
    """Randomly subsample a single linked SNP."""

    prevName = ""
    Unlink = []
    TempLink = []
    Num = 1

    for name in Bi_Thr:
    
        GeneName = name[0].strip()
    
        if Num == 1:
            prevName = GeneName
            TempLink.append(name)
            Num += 1
        
        elif Num > 1 and Num < len(Bi_Thr):
        
            #print "Middle SNPs: %s\tNum: %d" % (GeneName, Num)
            if GeneName != prevName:
            
                #Randomly sample single SNP from TempLink, 
                #and append to Unlink.
                Single = random.choice(TempLink)
                Unlink.append(Single)
            
                #Empty list.
                TempLink = []
                prevName = GeneName
            
                #Start new Temp list with this SNP.
                TempLink.append(name)
                Num += 1
            
            else:
        
                TempLink.append(name)
                prevName = GeneName
                Num += 1
   
        elif Num == len(Bi_Thr):
        
            if GeneName != prevName:
                Single = random.choice(TempLink)
                Unlink.append(Single)
                Unlink.append(name)
            
            else:
                TempLink.append(name)
                Single = random.choice(TempLink)
                Unlink.append(Single)
    return Unlink

def create_AFS(Unlink_Matrix, Pops, individuals, pop1, pop2, 
               pop1_thr, pop2_thr, AFS):
    """Count the minor allele in each population 
       and populate the observed AFS."""

    for snp in Unlink_Matrix:

        Allowed = ['A', 'C', 'G', 'T']

        alleles = set([x for x in snp[1:] if x in Allowed])
        Total = ''.join(alleles)

        #Try the Count alleles as list comprehensions.
        Pop1 = [snp[i] for i in range(len(snp)) if i > 0 and snp[i] 
                in Allowed and Pops[individuals[i]] == pop1]
        Pop2 = [snp[i] for i in range(len(snp)) if i > 0 and snp[i] 
                in Allowed and Pops[individuals[i]] == pop2]
        
        
        #Determine the minor allele.
        for c in Total:
            Allele = 100 * snp[1:].count(c) / \
                     float(len([x for x in snp[1:] if x in Allowed]))
            
            # If there is a minor allele
            if Allele < 50:
                
                MinP1 = 0
                MinP2 = 0
                
                #print "SNP # %s" % snp[0]
                if len(Pop1) == pop1_thr:
                    MinP1 = Pop1.count(c)
                    
                elif len(Pop1) > pop1_thr:
                    SubP1 = []
                    for sub in range(pop1_thr):
                        ch = random.choice(Pop1)
                        SubP1.append(ch)
                    MinP1 = SubP1.count(c)
                        
                if len(Pop2) == pop2_thr:
                    MinP2 = Pop2.count(c)
                    
                elif len(Pop2) > pop2_thr:
                    SubP2 = []
                    for sub in range(pop2_thr):
                        ch = random.choice(Pop2)
                        SubP2.append(ch)
                    MinP2 = SubP2.count(c)

                AFS[MinP2][MinP1] += 1.0
                
            elif Allele == 50:
                
                #All1 and All2 are the two alleles.
                All1 = Total[0]
                All2 = Total[1]
    
                MinP1_T0 = 0
                MinP2_T0 = 0
            
                MinP1_T1 = 0
                MinP2_T1 = 0
            
                if len(Pop1) == pop1_thr:
                    MinP1_T0 = Pop1.count(All1)
                    MinP1_T1 = Pop1.count(All2)
                    
                elif len(Pop1) > pop1_thr:
                    SubP1 = []
                    for sub in range(pop1_thr):
                        ch = random.choice(Pop1)
                        SubP1.append(ch)
                    MinP1_T0 = SubP1.count(All1)
                    MinP1_T1 = SubP1.count(All2)
                        
                if len(Pop2) == pop2_thr:
                    MinP2_T0 = Pop2.count(All1)
                    MinP2_T1 = Pop2.count(All2)
                    
                elif len(Pop2) > pop2_thr:
                    SubP2 = []
                    for sub in range(pop2_thr):
                        ch = random.choice(Pop2)
                        SubP2.append(ch)
                    MinP2_T0 = SubP2.count(All1)
                    MinP2_T1 = SubP2.count(All2)

                AFS[MinP2_T0][MinP1_T0] += 0.5
                AFS[MinP2_T1][MinP1_T1] += 0.5
                break 
                
    TotalBi_SNPs = 0
    
    for i in AFS:
        for j in i:
            TotalBi_SNPs += j
    return AFS, TotalBi_SNPs, TotalBi_SNPs - AFS[0][0]

def Count_loci(Monomorphics, threshold):
    """Count number of individuals sampled per locus. Retain those \
       loci that are >= the threshold."""
    
    with open(Monomorphics, 'r') as Infile:
    
        LineNumber = 0
        TotalLoci = 0
        TotalBP = 0
        
        MonoLociLengths = []

        for i in Infile:
            line = i.strip().split('\t')
            MonoLociLengths.append(len(line[0]))
            
            count = 0    
            for j in line[1:]:
        
                #If locus was called in an individual, count that.
                if int(j) > 0:
                    count += 1
            
            #If enough individuals were sequenced to retain a locus, 
            #count that.    
            if count >= int(float(threshold) / 100 * len(line[1:])):
                TotalLoci += 1
                
                #For each locus retained, count the total number of bps
                TotalBP += len(line[0])
            
            LineNumber += 1
        return TotalLoci, TotalBP + 6 * TotalLoci, MonoLociLengths

def write_out_AFS(AFS, Mono_Cell, Replicate):
    """Write AFS to an output file for analysis."""
    
    with open("Rep" + str(Replicate) + "_jointMAFpop1_0.obs", 'w') as Out:
        header = "1 observations"
        Out.write(header + '\n' + '\t')
        
        #Add the monomorphic sites to cell[0][0]
        AFS_comp = AFS
        AFS_comp[0][0] = Mono_Cell

        #Write the header for population 1
        for i in range(len(AFS_comp[0])):
            if i < len(AFS_comp[0]) - 1:
                head = "d0_" + str(i) + '\t'
            else:
                head = "d0_" + str(i) + '\n'
            Out.write(head)
            
        #Write the AFS table
        for i in range(len(AFS_comp)):
            for j in range(len(AFS_comp[i])):
                if j == 0:
                    Out.write('d1_' + str(i) + '\t' + 
                              str(AFS_comp[i][j]) + '\t')
                elif j == len(AFS_comp[i]) - 1:
                    Out.write(str(AFS_comp[i][j]) + '\n')
                else:
                    Out.write(str(AFS_comp[i][j]) + '\t')

def totalbp(file):
    """Count total number of sequenced base pairs"""
    
    with open(file, 'r') as Infile:
        loci = ''
        Count = 0
        Loci_count = 0
        
        start = False
        for line in Infile:
            
            #count first locus
            if Loci_count == 0:
                line = line.strip().split()
                Count += len(line[1]) + 6
                Loci_count += 1
            
            #find line breaks
            if '//' in line:
                start = True
            else:
                if start == True and line.startswith(">"):
                    line = line.strip().split()
                    #6 takes into account RE site
                    Count += len(line[1]) + 6
                    Loci_count += 1
                    start = False                
    return Count, Loci_count

def get_mono_cell(locus_file, TotalSNPs, TotalBi_SNPs_used):
    """Determine value to add to [0,0] cell"""
    TotalBP, Loci_count = totalbp(locus_file)
    return int((TotalBi_SNPs_used * TotalBP) / TotalSNPs) - TotalBi_SNPs_used, \
           TotalBP, Loci_count


def get_monomorphic_cell(TotalPolyLoci, TotalSNPs, TotalBi_SNPs_used, 
                         TotalMonoLoci, MonoLociLengths):
    """Determine value to add to [0,0] cell to standardize matrix
       and allow calculation of real parameter values with fsc. Used
       for unlinked SNP data set"""
    RE_length = 6
    TotalLoci = TotalPolyLoci + TotalMonoLoci
    
    TotalBP = 0
    for i in range(TotalLoci):
        bp = random.choice(MonoLociLengths)
        TotalBP += bp + RE_length
    #Get ratio of mono/poly sites. Multiply by BI_SNPs used to populate
    #the AFS to scale the invariant sites in the genome
    return int((TotalBi_SNPs_used * TotalBP) / TotalSNPs) - TotalBi_SNPs_used, \
           TotalBP


def Main_pyRAD():
    """Commands to run"""

    os.mkdir("./AFS/")
    os.mkdir("./log_files/")
    
    for rep in range(int(sys.argv[5])):
        Pops, Pop_counts = pop_association(sys.argv[2])
        AFS, pop1, pop2, pop1_thr, pop2_thr = build_AFS_LoL(file, Pop_counts, 
                                                            int(sys.argv[4]))
        Bi_Thr, individuals, TotalSNPs, TotalPolyLoci = Biallelic_SNPs(
            sys.argv[1], Pops, pop1, pop2, pop1_thr, pop2_thr)
        
        unlinked_data_set = True
        
        if sys.argv[6].upper().startswith('U'):
            Unlinked = subsample(Bi_Thr)
            AFS_fill, TotalBi_SNPs, TotalBi_SNPs_used = create_AFS(
                Unlinked, Pops, individuals, pop1[0], pop2[0], 
                pop1_thr, pop2_thr, AFS)
        else:
            unlinked_data_set = False
            AFS_fill, TotalBi_SNPs, TotalBi_SNPs_used = create_AFS(
                Bi_Thr, Pops, individuals, pop1[0], pop2[0], 
                pop1_thr, pop2_thr, AFS)
        
        mono_sites, totalBP, Loci_count = get_mono_cell(sys.argv[3], TotalSNPs,
                                                        TotalBi_SNPs_used)
            
        write_out_AFS(AFS, mono_sites, rep)
        shutil.move("./Rep" + str(rep) + "_jointMAFpop1_0.obs", 
                    "./AFS/Rep" + str(rep) + "_jointMAFpop1_0.obs")
            
        #Logging variables
        Pop1_samples = [i for i in individuals[1:] if Pops[i] == pop1[0]]
        Pop2_samples = [i for i in individuals[1:] if Pops[i] == pop2[0]]

        with open("Rep" + str(rep) + "_log.txt", 'w') as log:
            log.write("pyRAD input file\n\n")
            log.write("Threshold: %s%%\n" % sys.argv[4])
            if unlinked_data_set == True:
                log.write("Data set: Unlinked SNPs\n\n\n")
            else:
                log.write("Data set: Full SNPs\n\n\n")
            log.write("Population 1: %s\n" % pop1[0])
            log.write("Samples: %s\n" % Pop1_samples)
            log.write("Total Alleles: %s\tTotal Alleles Subsampled: %d\n\n\n" %
                       (pop1[1], pop1_thr))
            log.write("Population 2: %s\n" % pop2[0])
            log.write("Samples: %s\n" % Pop2_samples)
            log.write("Total Alleles: %s\tTotal Alleles Subsampled: %d\n\n\n" %
                       (pop2[1], pop2_thr))
            #log.write("Poly loci: %d\nMono Loci: %d\nTotal Loci: %d\n\n\n" %
            #         (TotalPolyLoci, TotalMonoLoci, TotalMonoLoci + TotalPolyLoci))
            log.write("Total Loci: %d\n" % Loci_count)
            log.write("Total SNPs: %d\t\t\tFinal SNPs: %d\nTotal BPs: %d\t\t\tAdjusted BPs: %d" %
                     (TotalSNPs, TotalBi_SNPs_used, totalBP, mono_sites + TotalBi_SNPs_used))
        shutil.move("./Rep" + str(rep) + "_log.txt", "./log_files/Rep" + str(rep) + "_log.txt")

def Main_AftrRAD():
    """Commands to run"""

    os.mkdir("./AFS/")
    os.mkdir("./log_files/")
    
    for rep in range(int(sys.argv[5])):
        Pops, Pop_counts = pop_association(sys.argv[2])
        AFS, pop1, pop2, pop1_thr, pop2_thr = build_AFS_LoL(file, Pop_counts, 
                                                            int(sys.argv[4]))
        Bi_Thr, individuals, TotalSNPs, TotalPolyLoci = Biallelic_SNPs(
            sys.argv[1], Pops, pop1, pop2, pop1_thr, pop2_thr)
        
        unlinked_data_set = True
        
        if sys.argv[6].upper().startswith('U'):
            Unlinked = subsample(Bi_Thr)
            AFS_fill, TotalBi_SNPs, TotalBi_SNPs_used = create_AFS(
                Unlinked, Pops, individuals, pop1[0], pop2[0], 
                pop1_thr, pop2_thr, AFS)
        else:
            unlinked_data_set = False
            AFS_fill, TotalBi_SNPs, TotalBi_SNPs_used = create_AFS(
                Bi_Thr, Pops, individuals, pop1[0], pop2[0], 
                pop1_thr, pop2_thr, AFS)
        
        TotalMonoLoci, TotalMonoBP, MonoLociLengths = Count_loci(sys.argv[3], 
                                                                 int(sys.argv[4]))
        
        Mono_Cell, TotalBP = get_monomorphic_cell(
            TotalPolyLoci, TotalSNPs, TotalBi_SNPs_used,
            TotalMonoLoci, MonoLociLengths)

        write_out_AFS(AFS, Mono_Cell, rep)
        shutil.move("./Rep" + str(rep) + "_jointMAFpop1_0.obs", 
                    "./AFS/Rep" + str(rep) + "_jointMAFpop1_0.obs")
            
        #Logging variables
        Pop1_samples = [i for i in individuals[1:] if Pops[i] == pop1[0]]
        Pop2_samples = [i for i in individuals[1:] if Pops[i] == pop2[0]]

        with open("Rep" + str(rep) + "_log.txt", 'w') as log:
            log.write("AftrRAD input file\n\n")
            log.write("Threshold: %s%%\n" % sys.argv[4])
            if unlinked_data_set == True:
                log.write("Data set: Unlinked SNPs\n\n\n")
            else:
                log.write("Data set: Full SNPs\n\n\n")
            log.write("Population 1: %s\n" % pop1[0])
            log.write("Samples: %s\n" % Pop1_samples)
            log.write("Total Alleles: %s\tTotal Alleles Subsampled: %d\n\n\n" %
                       (pop1[1], pop1_thr))
            log.write("Population 2: %s\n" % pop2[0])
            log.write("Samples: %s\n" % Pop2_samples)
            log.write("Total Alleles: %s\tTotal Alleles Subsampled: %d\n\n\n" %
                       (pop2[1], pop2_thr))
            log.write("Poly loci: %d\nMono Loci: %d\nTotal Loci: %d\n\n\n" %
                     (TotalPolyLoci, TotalMonoLoci, TotalMonoLoci + TotalPolyLoci))
            log.write("Total SNPs: %d\t\t\tFinal SNPs: %d\nTotal BPs: %d\t\t\tAdjusted BPs: %d" %
                     (TotalSNPs, TotalBi_SNPs_used, TotalBP, Mono_Cell + TotalBi_SNPs_used))
        shutil.move("./Rep" + str(rep) + "_log.txt", "./log_files/Rep" + str(rep) + "_log.txt")

if __name__ == '__main__':
    if len(sys.argv) != 7:
        print "python AFS_FSC_total.py SNP_infile.txt traits.txt", \
               "Monomorphics.txt/species.loci Threshold Replicate Un/Linked"
        sys.exit()

    if sys.argv[3].endswith(".loci"):
        Main_pyRAD()
    else:
        Main_AftrRAD()
