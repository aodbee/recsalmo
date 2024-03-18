# RECSALMO (Rapid Typing and Characterization Tool for Whole Genome Sequencing Data of Salmonella)

**Introduction:**

Prevalent foodborne pathogens that is Salmonella bacteria can cause serious illness which is mainly gastroenteritis with symptoms ranging from diarrhea, vomiting, nausea, muscle ache, fever and stomach cramps. In addition, the bacteria can cause a wide-scale outbreak; and their sources can be either plants, animals or humans such as onions, tomatoes, pigs and chickens. As such, performing epidemiological surveillance and outbreak control is essential in containing those Salmonella outbreaks. 
This project is intended to provide a python package for rapidly analyzing genome assemblies of Salmonella bacteria. This fast analysis of Salmonella genomes includes sequence-based typing (serotyping, MLST, cgMLST), antimicrobial Resistance genes, antibiotic classes/subclasses, phylogenetic tree construction, determination of Salmonella Pathogenicity Islands (SPI), and CRISPR determination.
This tool is intended to work only in Linux-like operating system e.g. Ubuntu. 


**Instructions:**

1. All genome assembly files must be contained in a single folder; and that folder must not have any sub-folders.
2. The assembly files must be of the FASTA format.
3. For each name of a genome assembly file, the program assigns a string in front of the first period symbol to be the main output name. For instances, if the name of an assembly file is Assembly1.scaffolds.fasta, the main output name (and also the output sub-folder name) becomes "Assembly1". Thus, names should be carefully adjusted to avoid duplicates that lead to errors.
4. The user must assign a single folder (in the command line) as the output folder. All the analysis results will be kept in this folder.


**Database:**

The followings are the embedded databases necessary for the analysis run of Salmonella genomes.
1.	Reference genome: Salmonella enterica subsp. enterica serovar Typhimurium str. LT2
2.	Reference genome: Salmonella enterica subsp. enterica serovar Typhimurium str. CT18
3.	Reference genome: Salmonella enterica subsp. enterica serovar Gallinarum str. 287/91
4.	SPIs (Salmonella Pathogenicity Islands): SPI-1 to SPI-17
5.	CRISPR spacers with assigned nametags


**Outputs:**

The outputs comprise both raw analytical reports (from dependent/custom packages and external software) and one summary file in “xlsx” and "csv" format. All the columns in the summary file include    
1.	Genome Assembly Name
2.	Salmonella sub-species
3.	Salmonella serovar
4.	MLST (Multi-Locus Sequence Typing)
5.	cgMLST (core-genome MLST)
6.	AMR:Gene/Antibiotic
7.	SPI (Salmonella Pathogenicity Islands) 
8.	Spacer-C1 (Spacer list in CRISPR locus 1)
9.	Spacer-C2 (Spacer list in CRISPR locus 2)
10.	NumSP-C1 (Number of spacers in CRISPR locus 1)
11.	NumSP-C2 (Number of spacers in CRISPR locus 2)
12.	DR-C1 (Direct repeat of CRISPR locus 1)
13.	DR-C2 (Direct repeat of CRISPR locus 2)
14.	Pos-C1 (Start-Stop Positions of CRISPR locus 1)
15.	Pos-C2 (Start-Stop Positions of CRISPR locus 2)
16.	LenC1 (Length in base-pair of  CRISPR locus 1)
17.	LenC2 (Length in base-pair of  CRISPR locus 2) 
    
Other than the summary file, there are also phylogenetic trees and pie-charts of all the genomes as follows
1. SNP-based Phylogenetic tree: this dendrogram is created by ParSNP tool with the Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 as the reference genome (NCBI Genome ID = SO4698-09 and accession = LN999997.1). The phylogenetic tree file is in the *.ggr format which can be opened and visualized by the Gingr software.
2. CRISPR-based Phylogenetic tree: this dendrogram is created by the alignment of CRISPR spacers of both locus 1 and 2.
3. Pie chart of serovar
4. Pie chart of ST

For the raw analytical report files, each analysis of one genome is composed of essential files as follows.
1.	“mlst_result.txt”
2.	“amrfinder_result.txt”
3.	“SPI_finderResult.txt”
4.	“crispr_result.txt”
5.	"spProfile.txt"


**Installation:**

This project is written mainly in Python; there are several dependent packages and software. The recommended version of Python is 3.8 and above. To make the installation a smooth experience, follow the steps below:
1. Download and install miniconda (Linux version): https://docs.anaconda.com/free/miniconda/index.html
2. Install Java
   >sudo apt update
   
   >sudo apt install default-jdk
   
   >sudo apt install default-jre
3. Create a custom conda environment (replacing "myenv" with the name of your choice)
   >conda create --name myenv
   
   >conda activate myenv
4. Install dependent conda packages
   >conda install -c bioconda fastmlst,sistr_cmd,ncbi-amrfinderplus,parsnp

   >conda install -c conda-forge openpyxl,seaborn
5. Update databases for fastmlst and ncbi-amrfinderplus packages
   >fastmlst --update-mlst -t 1
   
   >amrfinder -u
6. Create a folder (any name is fine) and put all the files of the RECSALMO project inside that folder (recsalmo.py is inside the project, first level).


**Usage:**

Users need to supply the input folder containing the genome assembly files. And the main output folder should be supplied as well but it is optional. If the output folder is not given, the folder “out_file” will be created in the current working directory and used as the main output folder.
Users need to supply 
  1. Absolute path to recsalmo.py
  2. Absolute path to input folder
  3. Absolute path to output folder (optional)
The format of program call is below
>python /path/recsalmo.py –input /path/input_folder –output /path/output_folder

Example program call supposing that all paths are under /home 
>python /home/recsalmo/recsalmo.py –input /home/input_folder –output /home/output_folder


  
