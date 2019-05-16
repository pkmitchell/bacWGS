# bacWGS
Bacterial WGS pipeline for AHDC Molecular Diagnostics. Includes both read qc and main pipeline scripts, as well as the script and table needed to translate abricate output into AMR table. Also includes yml and spec files to set up conda environment.

## Basic Workflow
**Read QC**: the bacWGS_readQC.py script takes all fastq files (forward and reverse) as input, runs FastQC, parses the output, and returns a table with number of reads, average read length, average read Phred score, annd proportion of reads with an average Phred score of at least 30. Also includes a blank column in which to calculate estimated coverage based on read number, read length, and expected genome size.

If reads meet QC standards, proceed using the forward read only as input for bacWGS_pipeline.py, which does the following:

1. **Trimming**: By default, reads are trimmed using the ILLUMINACLIP (NexteraPE-PE.fa:2:30:10) and SLIDINGWINDOW (window size 4, Q 20) programs in trimmomatic. This behavior can be disabled using the --no_trim (-t) flag (if, for example, you've already trimmed your input reads using some other method).

2. **Assembly**: Reads are assembled using SPAdes. Once the pipeline completes, the final contig and spades log files will be moved to the directory from which the scrpipt was launched and renamed to include the file stem. Assembly statistics are calculated using Quast and concatenated into a simplified report named assembly_stats.tsv.

3. **Mapping**: By default, reads are mapped back to the assembly in order to calculate coverage depth and insert size using BBMap. These metrics and some assembly statistics from Quast are then merged into a file called AssemblyQC_<datetime>.tsv. This behavior can be disabled using the --no_qc (-q) flag.
 
 **OPTIONAL STEPS (included in bacWGS_pipeline.py)**
 
 4. **Antimicrobial Resistance Gene Screening**: Conduct screening for AMR genes in the NCBI database with ABRicate and produce a nice summary output table that includes the gene name, identity, coverage, and relevant antibiotic class. Enabled by using the -a or --amr flag. The default coverage and identity thresholds are 50% and 85% by default, but these can be adjusted by using the --mincov (-c) and --minid (-i) flags. To use this option, a table linking the accession numbers output by ABRicate to the AMR gene and antibiotic class table must be either located in the default location (/workdir/miniconda3/envs/bacWGS/abricate_ncbi_table_050719.tsv) or provided using the --amr_table (-r) flag.
 
 5. ***Salmonella* Serotype Prediction**: Runs *Salmonella* serotyping using SISTR. Enabled using the --salm_sero (-s) flag.
 
 6. ***E. coli* Serotype Prediction**: Runs *E. coli* serotyping using ectyper. Enabled using the --ecol_sero (-e) flag.
