#!/usr/bin/ python3

'''
2024-03-08
Julia Lienard

Description:
The programs fiters the variant_summary.txt file containing the full database of ClinVar, into a pre-fitlered
table based on "Single nucleotide variant", the more accurate "Assembly GRCh38" and ClinicalSignificance as 
being "Pathogenic".
Lines where RS# (dbSNP) is missing (marked at "-1) or the phenotypelist is missing (marked as "not provided")
are excluded.
The output tab-delimited table contains 8 columns:
AlleleID\tGeneID\tClinicalSignificance\tRS#(dbSNP)\tPhenotypeList\tChromosome\tReferenceAlleleVCF\tAlternateAlleleVCF

Usage: FilterClinVar.py variant_summary.txt

'''
# 1. Import required modules
import sys
import os

# 2. check correct number of arguments provided according to usage
if len(sys.argv) != 2:
    print("Wrong number of arguments provided!\nUsage: python FilterCommonIDs.py file.map ClinVarfileToFilter.txt")
    quit()
inputClinVar = sys.argv[1]

# check if output files do not exist already
if os.path.isfile("Clinvar_SNPpatho.txt"):
        print("The output Clinvar_SNPpatho.txt already exists. Remove or Rename existing output file")
        quit()
if os.path.isfile("Clinvar_SNPpatho_IDs.txt"):
        print("The output Clinvar_SNPpatho_IDs.txt already exists. Remove or Rename existing output file")
        quit()

with open(inputClinVar, "r") as inputfile, \
open("Clinvar_SNPpatho.txt", "w")as output, open("Clinvar_SNPpatho_IDs.txt", "w") as IDs:
	output.write(f'#AlleleID\tGeneID\tClinicalSignificance\tRS#(dbSNP)\tPhenotypeList\tChromosome\tReferenceAlleleVCF\tAlternateAlleleVCF\n')
	header = inputfile.readline()
	if header.startswith("#"):
		pass
	else:
		print("No header line detected in input ClinVar file!")

	for line in inputfile:
		if not line.startswith("#"):
			infoline = line.strip().split("\t")
			if len(infoline) != 40:
				print("Warning: Input file has not the expected number of columns (=40).", 
					"If only columns 35 and above are modified, no error in the output is expected")
			AlleleID = infoline[0]
			TypeVariant = infoline[1]
			GeneID = infoline[3]
			Assembly = infoline[16]
			ClinicalSignificance = infoline[6]
			RS = infoline[9]
			PhenotypeList = infoline[13]
			Chromosome = infoline[18]
			RefAllele = infoline[32]
			AltAllele = infoline[33]
			if TypeVariant == "single nucleotide variant" and Assembly == "GRCh38" and ClinicalSignificance == "Pathogenic" and RS != "-1" and PhenotypeList != "not provided":
				output.write(f'{AlleleID}\t{GeneID}\t{ClinicalSignificance}\t{RS}\t{PhenotypeList}\t{Chromosome}\t{RefAllele}\t{AltAllele}\n')
				IDs.write(f'{RS}\n')
