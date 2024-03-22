#!/usr/bin/ python3

'''
Author Julia Lienard
Date : 2024-03-08

Description: the program takes the genotype file generated with PlinkFormat2GenotypeFiltered.py and 
compares for each SNP ID the user has, the alternative allele referenced in the pre-filtered ClinVar db, 
associated with a pathogenic phenotypes. If that is the case, the program outputs a tab-delimited table with header:
#GeneticID rsID[chromosome] Clinical_Significance Disease/Disorder AlleleID Nb_alleles_involved
and corresponding information.

Usage: python Ancient_ClinVarmarker.py genotype.txt ClinvarSNP_patho.txt

'''

import sys
import os

## check correct number of arguments provided according to usage
if len(sys.argv) != 3:
    print("Wrong number of arguments provided!\nUsage: python Ancient_ClinVarmarker.py genotype.txt ClinvarSNP_patho.txt")
    quit()

GenotypeFile = sys.argv[1]
ClinVardb = sys.argv[2]


path="./ClinVar_markers/" 
if os.path.exists(path) == False:
	os.mkdir(path) # if no output directory exists, one is created. The ouput files will be directed to ClinVar_markers/ dir.

## check if output files do not exist already
if os.path.isfile("./ClinVar_markers/Ancient_markers.txt"):
	print("The output file called Ancient_markers.txt already exists. Remove or rename existing output file")
	quit()

with open(GenotypeFile) as Genotype, open("./ClinVar_markers/Ancient_markers.txt", "w") as marker:
	marker.write(f'#ancientGeneticID\trsID[chromosome]\tClinical_Significance\tDisease/Disorder\tAlleleID\tNb_alleles_involved\n')
	for line in Genotype:
		if not line.startswith("#"):
			infoline = line.strip().split("\t")
			if len(infoline) != 5:
				print("the genotype file provided is incorrect. Number of columns not equals to 5, with following structure:"
					"#rsiD chromosome position genotype UserID")
				quit()
			IDrs = infoline[0]
			ID = IDrs[2:]
			chromosomeUser = infoline[1]
			position = infoline[2]
			genotype = infoline[3]
			GeneticID = infoline[4]
			allele1 = genotype[0:1]
			allele2 = genotype[1:]
			AlleleNb = 0
			with open(ClinVardb, "r") as Clinvar:
				for Line in Clinvar:
					if Line.startswith("#"):
						header = Line.strip().split("\t")
						RScolumn = header[3]
						if RScolumn != "RS#(dbSNP)":
							print("The input ClinVar file is not a pre-filtered format with RS#(dbSNP) in the 4th columns " 
								"and following structure: #AlleleID	GeneID	ClinicalSignificance	RS#(dbSNP)	PhenotypeList	Chromosome	ReferenceAlleleVCF	AlternateAlleleVCF")
							quit()
					if not Line.startswith("#"):
						Infoline = Line.strip().split("\t")
						Rs = Infoline[3]
						chromosomeClinVar = Infoline[5]
						AltAllele = Infoline[7]
						AlleleID = Infoline[0]
						ClinSig = Infoline[2]
						Disease = Infoline[4]
						if ID == Rs and chromosomeUser == chromosomeClinVar:
							if allele1 == AltAllele:
								AlleleNb += 1
							if allele2 == AltAllele:
								AlleleNb += 1
							if AlleleNb == 1:
								marker.write(f'{GeneticID}\trs{ID}[chr{chromosomeUser}]\t{ClinSig}\t{Disease}\t{AlleleID}\t{AlleleNb}allele(s)\n')
								AlleleNb = 0
							if AlleleNb == 2:
								marker.write(f'{GeneticID}\trs{ID}[chr{chromosomeUser}]\t{ClinSig}\t{Disease}\t{AlleleID}\t{AlleleNb}allele(s)\n')
								AlleleNb = 0
						


