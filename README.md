# RECSALMO (Rapid Typing and Characterization Tool for Whole Genome Sequencing Data of Salmonella)

**Introduction:**

Prevalent foodborne pathogens that is Salmonella bacteria can cause serious illness which is mainly gastroenteritis with symptoms ranging from diarrhea, vomiting, nausea, muscle ache, fever and stomach cramps. In addition, the bacteria can cause a wide-scale outbreak; and their sources can be either plants, animals or humans such as onions, tomatoes, pigs and chickens. As such, performing epidemiological surveillance and outbreak control is essential in containing those Salmonella outbreaks. 
This project is intended to provide a python package for rapidly analyzing genome assemblies of Salmonella bacteria. This fast analysis of Salmonella genomes includes sequence-based typing (serotyping, MLST, cgMLST), antimicrobial Resistance genes, antibiotic classes/subclasses, phylogenetic tree construction, determination of Salmonella Pathogenicity Islands (SPI), and CRISPR determination.
This tool is intended to work only in Linux-like operating system e.g. Ubuntu. 


**Instructions:**

1. All genome assembly files must be contained in a single folder; and that folder must not have any sub-folders.
2. The assembly files must be of the FASTA format.
3. For each name of a genome assembly file, the program assigns a string in front of the first period symbol to be the main output name. For instances, if the name of an assembly file is Assembly1.scaffolds.fasta, the main output name (and also the output folder name) becomes "Assembly1". Thus, names should be carefully adjusted to avoid duplicates.
4. The user must assign a single folder (in the command line) as the output folder. All the analysis results will be kept in this folder.


**Database:**

There are embedded databases necessary for the analysis of Salmonella genome. The followings are built-in databases.
1.	Reference genome (for phylogenetic tree construction) which is Salmonella enterica subsp. enterica serovar Typhimurium str. LT2
2.	SPIs (Salmonella Pathogenicity Islands)
3.	CRISPR spacers


**Outputs:**

The outputs comprise both raw analytical reports (from dependent/custom packages and external software) and one summary file in “xlsx” format.  The displayed outputs in the summary file are    
1.	Identification of Salmonella species/sub-species
2.	Serotype (sequence-based)
3.	MLST (Multi-Locus Sequence Typing)
4.	cgMLST (core-genome MLST)
5.	Serovar cgMLST
6.	Serogroup
7.	cgMLST subspecies
8.	AMR classes (in the format of “CLASS1:SUBCLASS1|CLASS2:SUBCLASS2|…”
9.	SPI (Salmonella Pathogenicity Islands) 
10.	Number of CRISPR Locus
Other than the summary file, there are also the phylogenetic tree of all the genomes created by ParSNP tool using Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 as the reference genome (NCBI Genome ID = SO4698-09 and accession = LN999997.1). The phylogenetic tree file is in the *.ggr format which can be opened and visualized by the Gingr software. The phylogenetic tree-related files are in the folder “treeoutput”. Moreover, there two pie charts presenting Serotype and MLST. The summary file and pie chart images are in the “out_file” folder. 
For the raw analytical report files, each analysis of one genome is composed of essential files as follows.
1.	“SeqSero_result.txt”
2.	“mlst_result.txt”
3.	“amrfinder_result.txt”
4.	“SPI_result.txt”
5.	“crispr_result.txt”



**Usage:**

Users need to supply the input folder containing the genome assembly files. And the main output folder should be supplied as well but it is optional. If the output folder is not given, the folder “out_file” will be created in the current working directory and used as the main output folder.
Example of the program call 
>python recsalmo.py –input /input_folder/ –output /output_folder/  

**Installation:**

This project is written mainly in Python; there are several dependent packages and software. The list of those are as follows
- Secsero2
  Instruction on installation: https://anaconda.org/bioconda/seqsero2 
- SISTR
  Instruction on installation: https://anaconda.org/bioconda/sistr_cmd
- FastMLST
  Instruction on installation: https://anaconda.org/bioconda/fastmlst
- AMRFinderPlus
  Instruction on installation: https://anaconda.org/bioconda/ncbi-amrfinderplus
- ParSNP
  Instruction on installation: https://anaconda.org/bioconda/parsnp
- NCBI Blast+
  Instruction on installation: https://anaconda.org/bioconda/blast

  Alternatively, you can also install this on Ubuntu terminal directly:
  >sudo apt-get -y install ncbi-blast+

Notice that for all the required tools above, you can conveniently use conda to install (via bioconda channel). So, you can the following on the terminal
  >conda install -c bioconda seqsero2,sistr_cmd,fastmlst,ncbi-amrfinderplus,parsnp,blast

Other required python packages are as follows,
- biopython
- numpy
- seaborn
- matplotlib
- openpyxl
You can use pip3 command (since python3 is used) to install all of those python packages mentioned, like this
  >pip3 install biopython,numpy,seaborn,matplotlib,openpyxl

In summary, create a folder (any name is fine) and put all the files/folders above inside. Then do the followings on the Ubuntu terminal,
and you are good to go.

>conda install -c bioconda seqsero2,sistr_cmd,fastmlst,ncbi-amrfinderplut,parsnp,blast

>pip3 install biopython,numpy,seaborn,matplotlib,openpyxl




  
