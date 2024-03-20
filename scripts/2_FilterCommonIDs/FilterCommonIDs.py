#!/usr/bin/ python3

'''
Author Julia Lienard
Date 2024-03-01

Description: This program evaluates the common SNP ID between a .map file and a pre-filtered
ClinVar database, runs a second filtering step on the ClinVar file and output the list of SNP 
IDs to extract.

List of function:
- common_member

Procedure:
1) the program extracts the SNPs IDs from a .map file (column 1) and stores in a set()
2) It then finds SNPs IDs from ClinVar_SNPpatho.txt (filtered ClinVar database file) 
and stores in a second set()
3) The common SNP Ids are evaluated using the function "common_member"
4) The ClinVar_SNPpatho.txt is filtered based on these common found IDs: output file is called
"ClinVar_SNPpatho" + basename + "filtered.txt", where basename refers to the .map file provided.
3) A list of the common SNP IDs is output

Usage: python FilterCommonIDs.py file.map ClinVarfileToFilter.txt 

'''
# 1. Import required modules
import sys
import os


# 2. User Defined Functions
def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
    return (a_set & b_set)

## check correct number of arguments provided according to usage
if len(sys.argv) != 3:
    print("Wrong number of arguments provided!\nUsage: python FilterCommonIDs.py file.map ClinVarfileToFilter.txt")
    quit()
MapFile = sys.argv[1]
ClinVar = sys.argv[2]

# extracting the basename of input .map file to name the output files
FileName = os.path.basename(MapFile).split(".")
basename = FileName[0]

List1ID = set()
List2ID = set()
common_IDs = set()


# 3. Check for errors:

## check if output files do not exist already
if os.path.isfile(basename + "SNP_to_extract.txt"):
        print("The output list of SNP_to_extract already exists. Change input .map file or rename existing output file")
        quit()
if os.path.isfile("ClinVar_SNPpatho" + basename + "filtered.txt"):
        print("The output filtered ClinVar file already exists. Change input .map file or rename existing output file")
        quit()


with open(MapFile, "r") as MapList:
    for line in MapList:
        info = line.split("\t")
        SNP_ID = info[1].strip()
        formatID = SNP_ID[:2]
        if formatID != "rs": ## check the .map is in the proper format
            print("The input .map file is not the correct file or not in the right format: rs ID should be in the first column, with no header")
            quit()
        ID = SNP_ID[2:]
        List1ID.add(ID.strip())


with open(ClinVar, "r") as ClinIDs, open(basename + "SNP_to_extract.txt", "w") as SNPList:
    for Line in ClinIDs:
        if Line.startswith("#"):
            header = Line.strip().split("\t")
            RScolumn = header[3]
            if RScolumn != "RS#(dbSNP)":
                print("The input ClinVar file is not a pre-filtered format with RS#(dbSNP) in the 4th columns")
                quit()
        if not Line.startswith("#"):
            Infoline = Line.strip().split("\t")
            Rs = Infoline[3]
            List2ID.add(Rs.strip())

    common_IDs = common_member(List1ID, List2ID)
    for element in common_IDs:
        SNPList.write(f'rs{element}\n')

with open(ClinVar, "r") as ClinIDs, open("ClinVar_SNPpatho" + basename + "filtered.txt", "w") as ClinVar_UserReduced:
    ClinVar_UserReduced.write(f'#AlleleID\tGeneID\tClinicalSignificance\tRS#(dbSNP)\tPhenotypeList\tChromosome\tReferenceAlleleVCF\tAlternateAlleleVCF\n')
    for Line in ClinIDs:
        if not Line.startswith("#"):
            Infoline = Line.strip().split("\t")
            AlleleID = Infoline[0]
            GeneID = Infoline[1]
            ClinicalSignificance = Infoline[2]
            Rs = Infoline[3]
            PhenotypeList = Infoline[4]
            Chromosome = Infoline[5]
            RefAllele = Infoline[6]
            AltAllele = Infoline[7]
            for element in common_IDs:
                if element == Rs:
                    ClinVar_UserReduced.write(f'{AlleleID}\t{GeneID}\t{ClinicalSignificance}\t{Rs}\t{PhenotypeList}\t{Chromosome}\t{RefAllele}\t{AltAllele}\n')

    



