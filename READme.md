# READme

In the context of an exercise for the course DNA sequences analyses of the Master program Bioinformatics at Lund University, this project was developped.

The overall aim was to identify disease markers among Ancient people's available genomic data, using the database ClinVar. Then, identify disease markers for any given user's genotype file and see if any of the identified markers were shared with Ancient people.

Available datasets are:
- variant_summary.txt.gz file (version February 2024) from ClinVar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
- Allen Ancient DNA Resource (AADR, v54.1) files in plink format: v54.1_1240K_public (.bed, .bim and .fam)
- 5 users' test files

---------------------------------- DISEASE MARKERS AMONG ANCIENT PEOPLE

## STEP 1 - Extracting the data subsets for the analysis

### A- Extract relavant lines of ClinVar for the subsequent analysis, from the variant_summary.txt 

The updated version of the variant_summary.txt can be used as this file is updated monthly by ClinVar.

#### Use of the FilterClinVarNotprovided.py script to filter 
- "Single nucleotide variant" 
- the more accurate Assembly "GRCh38" 
- "Pathogenic" for the clinical significance
- RS# (dbSNP) as not missing (!= "-1") which represent 1,814,016 lines in the variant_summary.txt
- Disease description indicated only as "not provided" are filtered out
and output only certain relevant columns for the subsequent analysis (AlleleID, GeneID, ClinicalSignificance, RS# ID, PhenotypeList, Chromosome, RefAllele, AltAllele).

This will create a lighter version of the ClinVar database that can be checked for the extra information included if needed.

```sh
# in home directory of the project:
mkdir script
cd script
# The FilterClinVar.py script is located in this directory - Make the script executable (chmod +x)
# in home directory of the project:
mkdir ClinVar
cd ClinVar
mkdir 1_rawdata # the inital ClinVar file "variant_summary.txt.gz" is located in this directory
cd 1_rawdata
gunzip variant_summary.txt.gz
cd ../ClinVar
mkdir 2_Filtered
cd 2_Filtered
# Set soft link to FilterClinVar.py:
ln -s ~/script/FilterClinVarNoprovided.py
python FilterClinVar.py ../variant_summary.txt

# outputs are Clinvar_SNPpatho.txt and Clinvar_SNPpatho_IDs.txt (both with 47759 lines without the header line)
```

### B - Extract info in the AADR database for individuals corresponding to "Ancient people"

#### Install plink v1.90b6.12 https://anaconda.org/bioconda/plink
```sh
conda install bioconda::plink
```

### Organize the data
```sh
# in home directory
mkdir AADR_v54.1
cd AADR_v54.1
mkdir 1_rawdata
# this directory contains the raw data from AADR: v54.1_1240K_public.bed, .bim and .fam
mkdir 2_Ancient_filtered
# this directory contains the list of IDs corresponding to Ancient people in "Ancient_samples.txt"
```

#### uncompress files using plink

```sh
# uncompress using plink
cd AADR_v54.1/1_rawdata
plink --bfile v54.1_1240K_public --recode --out v54.1_1240K_public
# 2 ouput files called v54.1_1240K_public.ped and v54.1_1240K_public.map are generated
```

#### Extract data from the AADR database for individuals corresponding to "Ancient people" only
In the Ancient_samples.txt, there is a list of Individual IDs, that can be used to filter the v54.1_1240K_public files
```sh
# in home directory:
cd AADR_v54.1/2_Ancient_filtered
plink --bfile ../1_rawdata/v54.1_1240K_public --keep Ancient_samples.txt --recode --out Ancient
# Only Ancient people are kept in the generated Ancient.ped and Ancient.map files, and those files are directly uncompressed, using --recode
```

### C - Filtering the SNP IDs in common between Ancient.map and the filtered ClinVar database "Clinvar_SNPpatho.txt"

####  1) Filter "Clinvar_SNPpatho.txt"

20230318 update: new common list ID with Clinvar_SNPpatho.txt in the 2_filteredStrict/
```sh
ClinVarpatho=~/ClinVar_search/ClinVar/2_filteredStrict/Clinvar_SNPpatho.txt 

```sh
# in home directory of the project:
cd ~/AADR_v54.1/2_Ancient_filtered
mkdir ClinVar_filtered
cd ClinVar_filtered
AncientMap=../Ancient.map
ClinVarpatho=../../../ClinVar/2_filteredStrict/Clinvar_SNPpatho.txt 
python FilterCommonIDs.py $AncientMap $ClinVarpatho
# output are the filtered ClinVar_SNPpathoAncientfiltered.txt and AncientSNP_to_extract.txt # 45 SNP IDs in common
```

#### 2) Filter Ancient.map and Ancient.ped

```sh
# in ~/AADR_v54.1/2_Ancient_filtered/ClinVar_filtered:
Ancient=../Ancient
plink --file $Ancient --extract AncientSNP_to_extract.txt --recode --out Ancient_ClinVfilteredStrict
# output files are Ancient_ClinVfilteredStrict.map and Ancient_ClinVfilteredStrict.ped with the 45 common SNP Ids.
# Total genotyping rate is 0.552615
```

#### 3) Extract Genotype, SNP ID, position, chromosome and user ID from .ped and .map files for the Ancient peple.

```sh
# From plink format files to a genotype file for Ancient_ClinVfilteredStrict
# in ~/AADR_v54.1/2_Ancient_filtered/ClinVar_filtered:
# creating soft link for PlinkFormat2GenotypeFiltered.py in the working directory:
ln -s ~/script/PlinkFormat2Genotype.py
python PlinkFormat2GenotypeFiltered.py Ancient_ClinVfilteredStrict.map Ancient_ClinVfilteredStrict.ped Ancient_ClinVfilteredStrict_Genotype.txt
# output file is Ancient_ClinVfiltered_Genotype.txt
cat Ancient_ClinVfilteredStrict_Genotype.txt | wc -l # 250228 lines
```

## STEP 2 - Which Ancient people have ClinVar markers and output that to a table.

```sh
cd ~/AADR_v54.1/2_Ancient_filtered/
ln -s ../../script/Ancient_ClinVarmarker.py
python Ancient_ClinVarmarker.py ClinVar_filtered/Ancient_ClinVfiltered_Genotype.txt ClinVar_filtered/ClinVar_SNPpathoAncientfiltered.txt
# Ancient_markers.txt is output in the newly generated ClinVar_markers/ directory in /AADR_v54.1/2_Ancient_filtered/
# 7548 disease markers identified in Ancient people using the ClinVar database version
```

---------------------------------- DISEASE MARKERS FOR ANY USER

## STEP 1 - Which markers the User's has? and which are shared with Ancient people

Là le All_inOne_UserTest.py marche sauf que je pers des Ancient homozygotes dans mon compte.. apparu avec l'intro de la comparaison des
AlleleID ! il est 17h04 14 mars

All_inOne_User.py script 
with inputs:
- UsersFile => I NEED TO MODIFY THE SCRIPT SUCH THAT DIFFERENT USER'S INPUT WORK!!!!!!!!!!!!!!!!!!! + check that alleleID is the same User/Ancient for the marker
- Clinvar_SNPpatho.txt
- Ancient_markers.txt
Outputs in a sub-directory:
- Table with disease-associated SNPs for the user
- Table with disease-associated SNPs for the user shared with Ancient people

