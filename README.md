# bacWGS
Bacterial WGS pipeline for AHDC Molecular Diagnostics. Includes both read qc and main pipeline scripts, as well as a helper script to produce the AMR table. Also includes yml and spec files to set up conda environment.

## Basic Workflow
**Read QC**: the bacWGS_readQC.py script takes all fastq files (forward and reverse) as input, runs FastQC, parses the output, and returns a table with number of reads, average read length, average read Phred score, annd proportion of reads with an average Phred score of at least 30. Also includes a blank column in which to calculate estimated coverage based on read number, read length, and expected genome size.

If reads meet QC standards, proceed using the forward read only as input for bacWGS_pipeline.py, which does the following:

**Assembly**: Reads are assembled using SKESA. Assemblies FASTAs are stored in the directory from which the script was launced and log files will be moved to the SKESA_output directory. Assembly statistics are calculated using Quast and concatenated into a simplified report named assembly_stats.tsv.

**Mapping**: By default, reads are mapped back to the assembly in order to calculate coverage depth and insert size using BBMap. These metrics and some assembly statistics from Quast are then merged into a file called AssemblyQC_<datetime>.tsv. This behavior can be disabled using the --no_qc (-q) flag. More detailed BBMap output is stored in the BBMap_output directory. 
 
 **OPTIONAL SCREENING STEPS (included in bacWGS_pipeline.py)**
 
**Antimicrobial Resistance Gene Screening**: Conduct screening for AMR genes using NCBI AMRFinderPlus and produce a nice summary output table that includes the gene name, identity, coverage, and relevant antibiotic class. Enabled by using the -a or --amr flag. Does not include screening for the "plus" genes or point mutations. Full output is stored in the AMRFinder_output folder and a summary table called ARG_table_<datetime>.tsv is produced in the working directory. 
 
***Salmonella* Serotype Prediction**: Runs *Salmonella* serotyping using SISTR. Enabled using the --salm_sero (-s) flag. Full output is stored in the SISTR_output directory and a summary table called SALM_sero_<datetime> is produced in the working directory. 
 
***E. coli* Serotype Prediction**: Runs *E. coli* serotyping using ectyper. Enabled using the --ecol_sero (-e) flag. Full output is stored in the ectyper_out_<datetime> directory and a summary table called ECOL_sero_<datetime> is produced in the working directory. 
 
 **OTHER OPTIONS**
 
 **Trimming**: Read trimming using the ILLUMINACLIP (NexteraPE-PE.fa:2:30:10) and SLIDINGWINDOW (window size 4, Q 20) programs in    Trimmomatic can be run before assembly by using the --trim (-t) flag.
 
 **SPAdes**: Assemblies can be produced with SPAdes rather than SKESA by using the --spades (-d) flag. 
 
 ## Dependencies
 **NOTE**: Version numbers are those in the original conda environment. Other versions may or may not work.

 FastQC 0.11.8\
 SKESA 2.3.0\
 Quast 5.0.2\
 BBMap 38.58 (needed if you want to produce assembly QC report)\
 NCBI AMRFinderPlus 3.0.12 (only needed for AMR gene screening)\
 sistr_cmd 1.0.2 (only needed for _Salmonella_ serotyping)\
 ectyper 0.8.1 (only needed for _E. coli_ serotyping)\
 Trimmomatic 0.39 (only needed for trimming option)\
 SPAdes 3.13.0 (only needed if you want to use it instead of SKESA)\
