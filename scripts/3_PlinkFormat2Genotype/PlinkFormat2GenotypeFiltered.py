# !/ usr/bin/ python3

'''
2024-03-06
Julia Lienard

Description: the program generates from .map and .ped files a table with:
header: #rsiD\tchromosome\tposition\tgenotype\tUserID

Procedure:
1) From the .map file, information about the rsIDs, the corresponding chromosome location and genome position,
are extracted into 3 separate lists.
2) From the .ped file, the userID is collected in a variable. Each line of the file is then parsed after being split by occurrence
of a whitespace.
The genotype of the two alleles for the first SNP ID starts at index position 6 until parsing goes until the end of the line.
Each allele genotype is collected (as a number from the .ped file). Their value corresponds
to the index of the list "nucleotides", where values 0,1,2,3,4 corresponds to nucleotides 0,A,C,G,T. The conversion from
numbers to letters is done.
3) Then #rsiD,chromosome,position,genotype and tUserID are output in a line in a tab-delimited table.
The genotype is the combination of allele 1 and allele 2 (Ex: AC), in the same column. Only genotypes not equals to 00
are output.

Usage PlinkFormat2genotypeFiltered.py file.map file.ped outputfile

'''
# 1. Import required modules
import sys
import os

# 2. check correct number of arguments provided according to usage
if len(sys.argv) != 4:
    print("Wrong number of arguments provided!\nUsage PlinkFormat2genotypeFiltered.py file.map file.ped outputfile")
    quit()

# 3. attribute arguments
MapFile = sys.argv[1]
PedFile = sys.argv[2]
OutputFile = sys.argv[3]

# 4. check if output files do not exist already
if os.path.isfile(OutputFile):
        print("The output genotype file already exists. Remove or Rename existing output file")
        quit()

List_chromosome =[]
List_SNP = []
List_position = []
with open(MapFile, "r") as MAP:
	for line in MAP:
		info = line.split("\t")
		if len(info) != 4:
			print("the .map file provided is incorrect. Number of columns not equals to 4!")
			quit()
		chromosome = info[0].strip()
		SNP_ID = info[1].strip()
		formatID = SNP_ID[:2]
		if formatID != "rs":
			print("The input .map file is not the correct file or not in the right format: rsID should be in the first column, with no header")
			quit()
		position = info[3].strip()
		List_chromosome.append(chromosome)
		List_SNP.append(SNP_ID)
		List_position.append(position)


with open(PedFile, "r") as PED, \
open(OutputFile, "w") as output:
	output.write(f'#rsiD\tchromosome\tposition\tgenotype\tUserID\n')
	PositionOnLine = 6 #corresponds to first nucleotide of genotype data
	LastNucleotide = ""
	nucleotides = ["0", "A", "C", "G", "T" ]
	for line in PED:
		infoline = line.strip().split(" ")
		UserID = infoline[1].strip()
		PositionOnLine = 6
		LastNucleotide = int(len(infoline)) -1
		index = 0 # index of the lists generated from the .map file (index=0 is 1st line on .map)
		while PositionOnLine < int(LastNucleotide): 
			genotypeAllele1 = infoline[PositionOnLine].strip() 
			allele1 = nucleotides[int(genotypeAllele1)]
			PositionOnLine += 1
			genotypeAllele2 = infoline[PositionOnLine].strip()
			allele2 = nucleotides[int(genotypeAllele2)]
			Genotype = allele1 + allele2
			PositionOnLine += 1
			if Genotype != "00":
				output.write(f'{List_SNP[index]}\t{List_chromosome[index]}\t{List_position[index]}\t{Genotype}\t{UserID}\n')
				index += 1
			else:
				index += 1
			
