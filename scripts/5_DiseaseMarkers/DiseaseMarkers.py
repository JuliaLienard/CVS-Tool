#!/usr/bin/ python3

'''
2024-03-09
Julia Lienard

Description: The program searches in the provided user's genotype file disease markers in the pre-filtered ClinVar (8 columns file)
and output that in a table. Then, idenfied markers also found among Ancient people are output in a second table, as well as which Ancient
people data it concerns.

List of function:
- common_member

Procedure:
1) Extract rs# dbSNP IDs from the User's genotype file and from the pre-filtered ClinVar file provided in two sets.
2) common SNP IDs are collected in a 3rd set, using the function "common_member"
3) Using this list of common SNP Ids, the ClinVar file is further filtered in a ClinVartemp.txt temporary file
4) Using this list of common SNP Ids, the User's genotype file is also filtered in a UserGenotypeTemp.txt temporary file
5) Comparing UserGenotypeTemp.txt to the ClinVartemp.txt to identify markers, which are output in a first txt file, named using the 
basename of the User's genotype file followed by "markers", in a newly generated output/ directory.
6) The idenfied markers for the user are compared with markers identified among Ancient people. The name of the Ancient people are
written in a second txt file if shared markers are found, named using the basename of the User's genotype file followed by 
"MarkersharedAncient", in a newly generated output/ directory.

Usage: python DiseaseMarkers.py UserGenotypeFile Clinvar_SNPpathoStrict.txt AncientMarker.txt

'''

import sys
import os

## check correct number of arguments provided according to usage
if len(sys.argv) != 4:
    print("Wrong number of arguments provided!\nUsage: python DiseaseMarker.py UserGenotypeFile \
    	ClinVarFilt_8col.txt AncientMarker.txt")
    quit()

# check prior existnce of temporary files needed for the script
if os.path.isfile("ClinVartemp.txt"):
	print("A file called ClinVartemp.txt has been found in the working directory.\n\
	Remove or Rename existing file to avoid conflict with this program")
	quit()
if os.path.isfile("UserGenotypeTemp.txt"):
	print("A file called UserGenotypeTemp.txt has been found in the working directory.\n\
	Remove or Rename existing file to avoid conflict with this program")
	quit()

UsersGenotype = sys.argv[1]
ClinVardb = sys.argv[2]
AncientMarkerFile = sys.argv[3]

# splits the filename to get the extension
genotypeFileExtensionArr = UsersGenotype.split(".")
# extension stands in the last part
genotypeFileExtension = genotypeFileExtensionArr[-1]

if genotypeFileExtension == "csv" :
	print("CSV file obtained")
	fileDelimiter=","
elif genotypeFileExtension == "txt" :
	print("TXT file obtained")
	fileDelimiter="\t"
else :
	print("unsupported genotype file extension")
	quit()

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
    return (a_set & b_set)

print("Process starts")
ListID_User = set()
ListID_ClinVar = set()
common_IDs = set()

# STEP1: identifying the common SNP IDs between the user's genotype file and the pre-filtered ClinVar file
with open(UsersGenotype, "r") as User, open(ClinVardb) as ClinVar:
    for line in User:
    	if not line.startswith("#"):
    		infoline = line.strip().split(fileDelimiter)
    		IDrs = infoline[0]
    		ID = IDrs[2:]
    		ListID_User.add(ID.strip())
    for Line in ClinVar:
    	if not Line.startswith("#"):
    		Infoline = Line.strip().split("\t")
    		Rs = Infoline[3]
    		ListID_ClinVar.add(Rs.strip())
    common_IDs = common_member(ListID_User, ListID_ClinVar)
    print(f"{len(common_IDs)} genotyped SNPs in the user's file found in the ClinVar database")

# STEP2: Filterting the pre-filtered ClinVar file, based on the common SNP IDs found, in a temporary Clinvar file
print("Extracting these SNPs from the ClinVar database")
with open(ClinVardb) as ClinVar, open("ClinVartemp.txt", "w") as TempClinVar:
	for Line in ClinVar:
		if Line.startswith("#"):
			TempClinVar.write(Line)
		if not Line.startswith("#"):
			Infoline = Line.strip().split("\t")
			Rs = Infoline[3]
			for element in common_IDs:
				if element == Rs:
					TempClinVar.write(Line)

# STEP3: Filterting the user's genotype, based on the common SNP IDs found, in a temporary genotype file
print("Extracting these SNPs from the user's file")
with open(UsersGenotype) as User, open("UserGenotypeTemp.txt", "w+") as TempGenotype:
	next(User) # skip the header line 
	for line in User:
		if not line.startswith("#"):
			infoline = line.strip().split(fileDelimiter) # used the correct fileDelimiter based on the extension
			IDrs = infoline[0]
			ID = IDrs[2:]
			for element in common_IDs:
				if element == ID:
					if len(infoline) == 4: # last 2 columns are already merged 
						# we write and specify the delimiter (to use the same format)
						TempGenotype.write(f'{infoline[0]}\t{infoline[1]}\t{infoline[2]}\t{infoline[3]}\n')
						break # to shorten the computing time, once the SNPs ID has been found once in common_IDs, go to next line
					elif len(infoline) == 5: # The alleles are split in two columns that we then merge
						TempGenotype.write(f'{infoline[0]}\t{infoline[1]}\t{infoline[2]}\t{infoline[3]}{infoline[4]}\n')
						break # to shorten the computing time, once the SNPs ID has been found once in common_IDs, go to next line
					else :
						print("Input genotypeFile in incorrect format. Expecte format: rsID/chromosome/position and allele 1 and 2 in \
							the 4th column or split in the 4th and 5th columns")
						quit()


# STEP4: Creating the output directory
path="./output/" 
if os.path.exists(path) == False:
	os.mkdir(path) # if no output directory exists, one is created. The ouput files will be directed to output/ dir.

FileName = os.path.basename(UsersGenotype).split(".")
basename = FileName[0] # we extract the name of the file to be used for the output file containing the markers


# STEP5: Comparing the SNP IDs between User's temporary genotype file and the temporary Clinvar file
print("Extracting SNPs associated with disease according to the ClinVar database")
with open("UserGenotypeTemp.txt", "r") as FilteredGenotype, open("output/" + basename + "marker" + ".txt", "w") as marker:
	marker.write(f'#rsID[chromosome]\tClinical_Significance\tDisease/Disorder\tAlleleID\tNb_alleles_involved\n')
	for line in FilteredGenotype:
		if not line.startswith("#"):
			infoline = line.strip().split("\t")
			IDrs = infoline[0]
			ID = IDrs[2:]
			chromosomeUser = infoline[1]
			position = infoline[2]
			genotype = infoline[3]
			allele1 = genotype[0:1]
			allele2 = genotype[1:]
			AlleleNb = 0 # variable to evaluate if the genotype is homozygous or heterozygous for each SNP
			with open("ClinVartemp.txt", "r") as FilteredClinvar:
				for Line in FilteredClinvar:
					if not Line.startswith("#"):
						Infoline = Line.strip().split("\t")
						AlleleID = Infoline[0]
						Rs = Infoline[3]
						chromosomeClinVar = Infoline[5]
						AltAllele = Infoline[7]
						ClinSig = Infoline[2]
						Disease = Infoline[4]
						if ID == Rs and chromosomeUser == chromosomeClinVar:
							if allele1 == AltAllele:
								AlleleNb += 1
							if allele2 == AltAllele:
								AlleleNb += 1
							if AlleleNb == 1:
								marker.write(f'rs{ID}[chr{chromosomeUser}]\t{ClinSig}\t{Disease}\t{AlleleID}\t{AlleleNb}allele(s)\n')
								AlleleNb = 0
							elif AlleleNb == 2:
								marker.write(f'rs{ID}[chr{chromosomeUser}]\t{ClinSig}\t{Disease}\t{AlleleID}\t{AlleleNb}allele(s)\n')
								AlleleNb = 0

os.remove("ClinVartemp.txt")
os.remove("UserGenotypeTemp.txt")

print("Comparing with disease markers found among Ancient people")
with open("output/" + basename + "marker" + ".txt", "r") as Usermarker, open("output/" + basename + "MarkersharedAncient" + ".txt", "w") as sharedMarker:
	sharedMarker.write(f'#Marker\tAlleleID\tDisease/disorder\tNumber_of_Ancient_people_with_similar_marker\tAncient_name_IDs\n')
	ListAncientNames = []
	for line in Usermarker:
		if not line.startswith("#"):
			infoline = line.strip().split("\t")
			rsIDuser = infoline[0].strip()
			DiseaseUser = infoline[2].strip()
			AlleleIDuser = infoline[3].strip()
			NbAncient = 0
			with open(AncientMarkerFile, "r") as AncientMarker:
				for Line in AncientMarker:
					if not Line.startswith("#"):
						infoLine = Line.strip().split("\t")
						AncientID = infoLine[0].strip()
						RSAncient = infoLine[1].strip()
						DiseaseAncient = infoLine[3]
						AlleleIDAncient = infoLine[4]
						if rsIDuser == RSAncient and AlleleIDuser == AlleleIDAncient: 
							NbAncient += 1
							ListAncientNames.append(AncientID)
			if NbAncient != 0:
				NamesAncient = ','.join(ListAncientNames)
				sharedMarker.write(f'>{rsIDuser}\t{AlleleIDuser}\t{DiseaseUser}\t{NbAncient} Ancient people\t{NamesAncient}\n\n')

print("Done!")






