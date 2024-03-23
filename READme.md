# READme

CVS-Tool was developped to be able to identify disease markers among Ancient people's available genomic data, using the database ClinVar. Then, identify disease markers for any given user's genotype file and see if any of the identified markers were shared with Ancient people.

Available datasets are:
- variant_summary.txt.gz file (version February 2024) from ClinVar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
- Allen Ancient DNA Resource (AADR, v54.1) files in plink format: v54.1_1240K_public (.bed, .bim and .fam)
https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
In addition the AADR Annotation.xlsx provides information about all the samples available. Genetic ID for ancient people can be extracted
from this excel file (one ID per line) as it is provided in the file "Ancient_samples.txt"
- 5 Test_users Genotype files (Test1_DNA.txt, Test2.csv, Test3.csv, Test4_DNA.txt, Test5_DNA.txt) to be used to test CVS-Tool for identification of susceptibility disease markers.

# Setup

CVS-Tool is run on the terminal.
The following softwares are required:
- Python 3
- plink v1.90b6.12 https://anaconda.org/bioconda/plink
```sh
conda install bioconda::plink
```
If the necessary paths to plink does not exist, set it up following this example:
```sh
PATH=$PATH:path/To/plink/
```

See the file tree used here in the provided pdf file: "CVS_Tool_File_tree.pdf" 

-------------------------------------------------------------------------------------------------------------------------------------------
## DISEASE MARKERS AMONG ANCIENT PEOPLE
-------------------------------------------------------------------------------------------------------------------------------------------
## STEP 1 - Extracting the data subsets for the analysis

### A- Extract relevant lines of ClinVar for the subsequent analysis, from the variant_summary.txt.gz

The updated version of the variant_summary.txt can be used as this file is updated monthly by ClinVar.

#### Use of the FilterClinVar.py script to filter:
- "Single nucleotide variant" 
- the more accurate Assembly "GRCh38" 
- "Pathogenic" for the clinical significance
- RS# (dbSNP) as not missing (!= "-1") which represent 1,814,016 lines in the variant_summary.txt
- Disease description indicated only as "not provided" are filtered out

The output contains only certain relevant columns for the subsequent analysis (AlleleID, GeneID, ClinicalSignificance, RS# ID, PhenotypeList, Chromosome, RefAllele, AltAllele).
This will create a lighter version of the ClinVar database that can be checked for the additional information included if needed.

```sh
# in home directory of the project:
mkdir script
cd script
# All scripts such as FilterClinVar.py are located in this directory - Make the scripts executable is not already (chmod +x)

# in home directory of the project:
mkdir ClinVar
cd ClinVar
mkdir 1_rawdata # the inital ClinVar file "variant_summary.txt.gz" is located in this directory
cd 1_rawdata
gunzip variant_summary.txt.gz

# Create a filtered ClinVar database:
cd ../ClinVar
mkdir 2_Filtered
cd 2_Filtered
# Set soft link to FilterClinVar.py:
ln -s ~/script/FilterClinVar.py
python FilterClinVar.py ../1_rawdata/variant_summary.txt
```
output file is called Clinvar_SNPpatho.txt (version from version February 2024 gives 47759 lines without the header line)


### B - Extract info in the AADR database for individuals corresponding to "Ancient people"

#### Organize the data
```sh
# in home directory
mkdir AADR_v54.1
cd AADR_v54.1
mkdir 1_rawdata
# this directory contains the raw data from AADR: v54.1_1240K_public.bed, .bim and .fam

# In AADR_v54.1/:
mkdir 2_Ancient_filtered
# this directory contains the list of IDs corresponding to Ancient people in "Ancient_samples.txt"
```

#### uncompress files using plink

```sh
# uncompress using plink
cd AADR_v54.1/1_rawdata
plink --bfile v54.1_1240K_public --recode --out v54.1_1240K_public
```
2 ouput files called v54.1_1240K_public.ped and v54.1_1240K_public.map are generated


#### Extract data from the AADR database for individuals corresponding to "Ancient people" only

In the Ancient_samples.txt, there is a list of Genetic IDs, that can be used to filter the v54.1_1240K_public files to keep only ancient DNA data.

```sh
# in home directory:
cd AADR_v54.1/2_Ancient_filtered
plink --bfile ../1_rawdata/v54.1_1240K_public --keep Ancient_samples.txt --recode --out Ancient
```
Only Ancient people are kept in the generated Ancient.ped and Ancient.map files, and those files are directly uncompressed, using --recode


### C - Filtering the SNP IDs in common between Ancient.map and the filtered ClinVar database "Clinvar_SNPpatho.txt"

####  1) Filter "Clinvar_SNPpatho.txt"

```sh
# in home directory of the project:
cd ~/AADR_v54.1/2_Ancient_filtered
mkdir ClinVar_filtered
cd ClinVar_filtered
ClinVarpatho=~/ClinVar/2_filteredStrict/Clinvar_SNPpatho.txt 
python FilterCommonIDs.py ../Ancient.map $ClinVarpatho
```
output files are:
- ClinVar_SNPpathoAncientfiltered.txt 
- AncientSNP_to_extract.txt 
45 SNP IDs in common found

#### 2) Filter Ancient.map and Ancient.ped

```sh
# in ~/AADR_v54.1/2_Ancient_filtered/ClinVar_filtered:
Ancient=../Ancient
plink --file $Ancient --extract AncientSNP_to_extract.txt --recode --out Ancient_ClinVfilteredStrict
```
Output files are:
- Ancient_ClinVfilteredStrict.map
- Ancient_ClinVfilteredStrict.ped 
with the 45 common SNP Ids.
Total genotyping rate is 0.552615

#### 3) Create a genotype file for Ancient people for the common ID found, using the Ancient_ClinVfilteredStrict.map and.ped

The following information will be provided in the generated genotype file:
#rsiD	chromosome	position	genotype	ancientGeneticID

```sh
# From plink format files to a genotype file for Ancient_ClinVfilteredStrict.map and .ped:
cd ~/AADR_v54.1/2_Ancient_filtered/ClinVar_filtered:
# creating soft link for PlinkFormat2GenotypeFiltered.py in the working directory:
ln -s ~/script/PlinkFormat2GenotypeFiltered.py
python PlinkFormat2GenotypeFiltered.py Ancient_ClinVfilteredStrict.map Ancient_ClinVfilteredStrict.ped Ancient_ClinVfilteredStrict_Genotype.txt
```
Output file is Ancient_ClinVfiltered_Genotype.txt

```sh
cat Ancient_ClinVfilteredStrict_Genotype.txt | wc -l # here we get 250228 lines including header
```

## STEP 2 - Which Ancient people have ClinVar markers and output that to a table.

```sh
cd ~/AADR_v54.1/2_Ancient_filtered/
# Create soft link for required python script:
ln -s ../../script/Ancient_ClinVarmarker.py
python Ancient_ClinVarmarker.py ClinVar_filtered/Ancient_ClinVfiltered_Genotype.txt ClinVar_filtered/ClinVar_SNPpathoAncientfiltered.txt
```
An output directory is created, called ClinVar_markers/, in ~/AADR_v54.1/2_Ancient_filtered/
Ancient_markers.txt is created in this new directory.

7547 disease markers identified in Ancient people using the ClinVar database version February 2024

The output table looks like the following (first 5 lines):
```sh
#ancientGeneticID	rsID[chromosome]	Clinical_Significance	Disease/Disorder	AlleleID	Nb_alleles_involved
Ne30_genotyping_noUDG	rs5082[chr1]	Pathogenic	Hypercholesterolemia, familial, 1	32975	2allele(s)
Ne30_genotyping_noUDG	rs3755319[chr2]	Pathogenic	Lucey-Driscoll syndrome	27321	1allele(s)
Ne61_genotyping_noUDG	rs5082[chr1]	Pathogenic	Hypercholesterolemia, familial, 1	32975	2allele(s)
Ne35_genotyping_noUDG	rs5082[chr1]	Pathogenic	Hypercholesterolemia, familial, 1	32975	2allele(s)
```

The alleleIDs have been included because it is used later when comparing disease markers for any user with those found in Ancient people. Indeed, some disease markers can have the same rsID with different alleles (and thus different alleleIDs) and different associated Disease/disorder description. If the ClinVar database is filtered in the future by including also "likely pathogenic" for instance, CVS-Tool will still be able to accurately give shared markers between a user and Ancient people.

-------------------------------------------------------------------------------------------------------------------------------------------
## DISEASE MARKERS FOR ANY USER
-------------------------------------------------------------------------------------------------------------------------------------------
## Which markers the User's has? and which are shared with Ancient people

```sh
#in home directory:
mkdir Testusers
cd Testusers # Test users' files are in this directory
ln -s ../script/DiseaseMarkers.py # setting soft link for required python script
# Creating variable with paths to required input files:
AncientMarker=../AADR_v54.1/2_Ancient_filtered/ClinVar_markers/Ancient_markers.txt
ClinVar=../ClinVar/2_filteredStrict/Clinvar_SNPpatho.txt
# Finding disease markers for the user Test1
python DiseaseMarkers.py Test1_DNA.txt $ClinVar $AncientMarker
```
Outputs in a sub-directory, automatically created:
- Table with disease-associated SNPs for the user
(in the case example: Test1_DNAmarker.txt)
- Table with disease-associated SNPs for the user shared with Ancient people (each marker line starts with ">")
(in the case example: Test1_DNAMarkersharedAncient.txt)

Example stdout information during processing:
```sh
CSV file obtained
Process starts
1933 genotyped SNPs in the user's file found in the ClinVar database
Extracting these SNPs from the ClinVar database
Extracting these SNPs from the user's file
Extracting SNPs associated with disease according to the ClinVar database
Comparing with disease markers found among Ancient people
Done!
```

Example of output table for disease markers found in the user's genotype:
Here for Test1_DNAmarker.txt
```sh
#rsID[chromosome]	Clinical_Significance	Disease/Disorder	AlleleID	Nb_alleles_involved
rs5082[chr1]	Pathogenic	Hypercholesterolemia, familial, 1	32975	1allele(s)
rs3755319[chr2]	Pathogenic	Lucey-Driscoll syndrome	27321	2allele(s)
rs553668[chr10]	Pathogenic	Lipodystrophy, familial partial, type 8	33200	2allele(s)
```

Example of output table for the user's markers shared with Ancient people:
Here for Test1_DNAMarkersharedAncient.txt
```sh
#Marker	AlleleID	Disease/disorder	Number_of_Ancient_people_with_similar_marker	Ancient_name_IDs
>rs5082[chr1]	32975	Hypercholesterolemia, familial, 1	3296 Ancien people	Ne30_genotyping_noUDG,Ne61_genotyping_noUDG,Ne35_genotyping_noUDG,I13833,CAO009013,I13580,I6445,I6441,I1707,I13778,I0231,I1525,I14688,I14690,I14692,I14689,I13838,I13840,I13839,I14622,I14685,I15706,I14687,AltaiNeanderthal.DG,AltaiNeanderthal_snpAD.DG,Aconcagua_noUDG.SG,I0308,I2230,I12941,I12942,I12943,Yaghan890_noUDG.SG,Yaghan894_noUDG.SG,I12376,MA577_noUDG.SG,I12362,I12354,I12363,I12365,I1630,I1659,I1660,I1635,ARM001,ARM002,I13599,I14339,I14341,I14340,I16116,I16191,I16219,I19322,I16707,I19327,I15733,I15734,I15749,I18479,I19331,I19332,I19333,I19334,I19335,I19336,I19338,I19339,I19343,I19344,I19345,I19346,I19348,I19350,I19352,I19353,I20439,I20443,I14051,I14052,I14588,I14601,I14602,I14603,I14606,I14618,I14620,I15748,I15750,I18162,I18167,I18248,I18472,I18478,I18481,I18486,I19323,I19324,I19325,I19326,I19328,I14065,I16376,I14619,I14621.....
```



