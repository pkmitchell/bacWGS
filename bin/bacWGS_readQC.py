#!/workdir/miniconda3/envs/bacWGS/bin/python

import sys
import string
import time
from subprocess import call
import argparse
import re

'''
CHANGES
-Changed default resource allocation to be functional for a single isolate
-Added Illuminaclip NexteraPE to default trimming
-Added option for AMR screening
-Added option for Salmonella serotyping
-Added option of E. coli serotyping
-Automated removal of unpaired reads following trimming
-Automated concatenation of assembly statistics into single file
-
-TODO: Echo versions to log
-TODO: Add resource multiplier option? -Rebuild such that resource allocation is per job, then set soft caps (give myself option to override?)
-TODO: More cleanup: trimmed reads, spades corrected reads/other spades output aside from the assembly
-TODO: Replace absolute path program calls with setting environment to include in path in submission script
-TODO: Long-term - use python multithreading rather than perl_fork_univ
-TODO: Add bbmap coverage stats
-TODO: Add fastqc step with QC metrics - Run first/separately, check if okay, then proceed? No, do in biopython, integrated into the job loop
-TODOING: Change thresholds to 85% identity, 50% coverage (set as defaults, allow to be modified)
'''

#Argument Parser
parser = argparse.ArgumentParser(description="Assemble genomes using SPAdes, or generate scripts to do so but do not submit immediately")
parser.add_argument("-l", "--launch", help="Launch job now (Default = Off)", action="store_true")
parser.add_argument("-t", "--no_trim", help="Turn off default read trimming (trimmomatic ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20)", action="store_true")
parser.add_argument("-q", "--no_qc", help="Turn off default mapping/assembly QC output", action="store_true")
parser.add_argument("-a", "--amr", help="Run AMR screeing on assemblies using ABRicate with NCBI database (Default = Off)", action="store_true")
parser.add_argument("-r", "--amr_table", help="Path to AMR lookup table (Default = /workdir/miniconda3/envs/bacWGS/abricate_ncbi_table_050719.tsv)", default="/workdir/miniconda3/envs/bacWGS/abricate_ncbi_table_050719.tsv")
parser.add_argument("-c", "--mincov", help="Minimum percent coverage for AMR gene detection. Must be between 0 and 100 (Default = 50)", type = float, default=50.0)
parser.add_argument("-i", "--minid", help="Minimum percent identity for AMR gene detection. Must be between 0 and 100 (Default = 85)", type = float, default=85.0)
parser.add_argument("-s", "--salm_sero", help="Run Salmonella serotype prediction on assemblies using SISTR (Default = Off)", action="store_true")
parser.add_argument("-e", "--ecol_sero", help="Run E. coli serotype prediction on assemblies using ectyper (Default = Off)", action="store_true")
parser.add_argument("-p", "--threads", help="Maximum number of processors to use (Default = 4)", type = int, default=4)
parser.add_argument("-n", "--jobs", help="Maximum number of jobs to run simultaneously, equivalent to perl_fork_univ.pl argument (Default = 1)", type = int, default=1)
parser.add_argument("-m", "--max_mem", help="Maximum total memory allocation to simultaneous processes in GB (Default = 24)", type = int, default=24)
parser.add_argument("Fastq", help="Input forward fastqs to be assembled", nargs='+')
args = parser.parse_args()


#Determining the number of jobs and threads per job
ncpus = args.threads #The total number of CPUs to be used simultaneously across all jobs
njobs = args.jobs #The number of jobs to be run simultaneously
totjobs = len(args.Fastq) #The total number of input sequences, which is also the total number of jobs to run

#Determining the number of threads per job and number of jobs to run simultaneously using the bounds provided by the user
if njobs > totjobs: 
	njobs = totjobs
	
memper = int(args.max_mem/njobs)

#This is a bit dumb and clunky, might be worth improving at some point to handle different input scenarios
if memper <=10:
	print("Seems like not enough memory for the number of jobs requested.\n")
	print("Reducing number of simultaneous processes")
	memper= 24
	njobs= int(args.max_mem/memper)
elif memper < 24:
	print("Might not be enough memory, trying anyway")

if njobs > ncpus:
	njobs = ncpus
	
threadsper = int(ncpus/njobs)

strT = time.strftime("%Y%m%d_%H%M%S")
#This creates an input file for perl_fork_univ.pl with one line per input sequence
#It is then called by the launch script generated in the next step

fn1="SPAdeSub_" + strT
ns=open(fn1,'w')

prefs = []

for i in range(0, totjobs):
	fs=args.Fastq[i]
	mtch= re.search(r'_R1_001.fastq',fs)
	if mtch:
		stem=re.split(r'_S[0-9][0-9]*_L001_R1_001.fastq', fs)[0]
		rs=fs.replace('_R1_', '_R2_')
	else:
		mtch2 = re.search(r'_1.fastq',fs)
		if mtch2:
			stem=fs.split("_1.fastq")[0]
			rs=fs.replace('_1.fastq', '_2.fastq')
	prefix=stem.split("/")[-1]
	prefs.append(prefix)
	
	if args.no_trim == False:
		ns.write("trimmomatic PE -threads " + str(threadsper) + " -phred33 " + fs + " " + rs + " trim_" + prefix + "_1.fastq.gz unpaired_" + prefix + "_1.fastq.gz trim_" + prefix + "_2.fastq.gz unpaired_" + prefix + "_2.fastq.gz ILLUMINACLIP:/workdir/miniconda3/envs/bacWGS/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 &>trim_" + prefix + ".log &&")
		ns.write(" rm unpaired_" + prefix + "_1.fastq.gz unpaired_" + prefix + "_2.fastq.gz &&")
		ns.write(" spades.py --pe1-1 " + "trim_" + prefix + "_1.fastq.gz --pe1-2 trim_" + prefix + "_2.fastq.gz -o ./" + prefix + "/ --careful -t " + str(threadsper) + " -m " + str(memper) + " &&") 
	else:
		ns.write(" spades.py --pe1-1 " + fs + " --pe1-2 " + rs + " -o ./" + prefix + "/ --careful -t " + str(threadsper) + " -m " + str(memper) + " &&")
	ns.write(" quast.py --threads " + str(threadsper) + " -o " + prefix + " --min-contig 1 -l \"" + prefix + "\" --contig-thresholds 500,1000 " + prefix + "/contigs.fasta && tail -n 1 " + prefix + "/transposed_report.tsv >>assembly_stats_nohead.tsv &&")
	ns.write(" mv " + prefix + "/contigs.fasta ./" + prefix + ".fasta && mv " + prefix + "/spades.log ./" + prefix + "_spades.log")
	if args.no_qc == False:
		ns.write(" && bbmap.sh nodisk ref=" + prefix + ".fasta in1=" + fs + " in2=" + rs + " covstats=" + prefix + "_covstats.txt t=" + str(threadsper) + " &>" + prefix + "_bbmap.log")
	if args.amr:
		ns.write(" && abricate --db ncbi --mincov " + str(args.mincov) + " --minid " + str(args.minid) + " --threads " +str(threadsper) + " " + prefix + ".fasta >" + prefix + "_abricate_ncbi.tsv 2>" + prefix + "_abricate_ncbi.log")
	ns.write(";\n")

ns.close()
 
#This step generates the bash script to set the environment and launch perl_fork_univ.pl, calling the previously generated input file 

fn2="./LaunchSPAdes_" + strT + ".sh"
ss=open(fn2, 'w')

ss.write("#!/bin/bash \n")
ss.write("\n")
ss.write("export PATH=/workdir/miniconda3/bin:$PATH\n\n")
ss.write("source activate bacWGS\n\n")
if args.no_trim == False:
	ss.write("echo \"trimmomatic $(trimmomatic -version)\" >>software_versions_" + strT + ".txt\n")
ss.write("spades.py -v >>software_versions_" + strT + ".txt\n")
ss.write("perl_fork_univ.pl " + fn1 + " " + str(njobs) + " &>" + fn1 + ".log \n")
ss.write("(head -n 1 " + prefs[0] + "/transposed_report.tsv && cat assembly_stats_nohead.tsv) >assembly_stats.tsv && rm assembly_stats_nohead.tsv\n")
if args.no_qc == False:
	ss.write("\nbbmap.sh --version 2>&1 >/dev/null| sed -n -e 's/^.*BBMap version /BBMap /p' >>software_versions_" + strT + ".txt\n")
	ss.write("echo -e \"Isolate\tTotal Length\tContigs\tN50\tMapped Coverage\tInsert Size\" >AssemblyQC_" + strT + ".tsv\n")
	ss.write("for s in " + ' '.join(prefs) + "; do i=$(grep $s assembly_stats.tsv |awk '{print $1,\"\t\",$8,\"\t\",$6,\"\t\",$10}');\n j=$(grep \"Average coverage\" ${s}_bbmap.log |cut -d$'\t' -f 2);\n k=$(grep \"insert size avg\" ${s}_bbmap.log |cut -d$'\t' -f 2);\n echo -e \"$i\t$j\t$k\";done >>AssemblyQC_" + strT + ".tsv\n")
if args.amr:
	ss.write("\nabricate -v >>software_versions_" + strT + ".txt\n")
	abrfiles = [p + "_abricate_ncbi.tsv" for p in prefs]
	ss.write("make_amr_table.py " + args.amr_table + " " + ' '.join(abrfiles) + " >ARG_table_" + strT + "\n")

if args.salm_sero:
	ss.write("\nsistr -V >>software_versions_" + strT + ".txt\n")
	fastas = [p + ".fasta" for p in prefs]
	ss.write("\nsistr --qc -f tab -t " + str(ncpus) + " -o Salm_seros_sistr_" + strT + " " + ' '.join(fastas) + " 2>Salm_seros_sistr_" + strT + ".log\n")
if args.ecol_sero:
	ss.write("\nectyper -V >>software_versions_" + strT + ".txt\n")
	fastas = [p + ".fasta" for p in prefs]
	ss.write("\nectyper -i " + ','.join(fastas) + " --cores " + str(ncpus) + " --verify -o ectyper_out_" + strT + " &>ectyper_" + strT + ".log\n")
ss.write("conda deactivate")
ss.close()
call(["chmod", "u+x", fn2])

#Optionally run the job
if args.launch:
	call([fn2])
else:
	print("To run, call " + fn2)