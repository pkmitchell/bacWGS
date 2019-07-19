#!/workdir/miniconda3/envs/bacWGS/bin/python

import sys
import string
import time
from subprocess import call
import argparse
import re
import multiprocessing as mp

'''
-TODO: Add resource multiplier option? -Rebuild such that resource allocation is per job, then set soft caps (give myself option to override?)
-TODO: Long-term - use python multithreading rather than perl_fork_univ
'''

#Argument Parser
parser = argparse.ArgumentParser(description="Assemble genomes using SKESA (or, optionally, SPAdes), or generate scripts to do so but do not submit immediately")
parser.add_argument("-l", "--launch", help="Launch job now (Default = Off)", action="store_true")
parser.add_argument("-t", "--trim", help="Turn on read trimming (trimmomatic ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20)", action="store_true")
parser.add_argument("-q", "--no_qc", help="Turn off default mapping/assembly QC output", action="store_true")
parser.add_argument("-a", "--amr", help="Run AMR screeing on assemblies using NCBI AMRFinder (Default = Off)", action="store_true")
parser.add_argument("-s", "--salm_sero", help="Run Salmonella serotype prediction on assemblies using SISTR (Default = Off)", action="store_true")
parser.add_argument("-e", "--ecol_sero", help="Run E. coli serotype prediction on assemblies using ectyper (Default = Off)", action="store_true")
parser.add_argument("-p", "--threads", help="Maximum number of processors to use (Default = 4)", type = int, default=4)
parser.add_argument("-n", "--jobs", help="Maximum number of jobs to run simultaneously, equivalent to perl_fork_univ.pl argument (Default = 1)", type = int, default=1)
parser.add_argument("-m", "--max_mem", help="Maximum total memory allocation to simultaneous processes in GB (Default = 24)", type = int, default=24)
parser.add_argument("-d", "--spades", help="Use SPAdes instead of SKESA for assembly (Default = Off)", action="store_true")
parser.add_argument("Fastq", help="Input forward fastqs to be assembled", nargs='+')
args = parser.parse_args()


#Determining the number of jobs and threads per job
ncpus = args.threads #The total number of CPUs to be used simultaneously across all jobs
tot_cpus = mp.cpu_count()
if ncpus > tot_cpus:
	ncpus = tot_cpus

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

fn1="joblist_" + strT
ns=open(fn1,'w')

prefs = []

for i in range(0, totjobs):
	fs=args.Fastq[i]
	rs=None
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
	
	if args.trim == True:
		ns.write("trimmomatic PE -threads " + str(threadsper) + " -phred33 " + fs + " " + rs + " trim_" + prefix + "_1.fastq.gz unpaired_" + prefix + "_1.fastq.gz trim_" + prefix + "_2.fastq.gz unpaired_" + prefix + "_2.fastq.gz ILLUMINACLIP:/workdir/miniconda3/envs/bacWGS/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 &>trim_" + prefix + ".log &&")
		ns.write(" rm unpaired_" + prefix + "_1.fastq.gz unpaired_" + prefix + "_2.fastq.gz &&")
		if args.spades == True:
			ns.write(" spades.py --pe1-1 " + "trim_" + prefix + "_1.fastq.gz --pe1-2 trim_" + prefix + "_2.fastq.gz -o ./" + prefix + "/ --careful -t " + str(threadsper) + " -m " + str(memper) + " &&")
			ns.write(" quast.py --threads " + str(threadsper) + " -o " + prefix + " --min-contig 1 -l \"" + prefix + "\" --contig-thresholds 500,1000 " + prefix + "/contigs.fasta && tail -n 1 " + prefix + "/transposed_report.tsv >>assembly_stats_nohead.tsv &&")
			ns.write(" mv " + prefix + "/contigs.fasta ./" + prefix + ".fasta && mv " + prefix + "/spades.log ./" + prefix + "_spades.log")
		else:
			ns.write(" skesa --cores " + str(threadsper) + " --memory " + str(memper) + " --fastq " + "trim_" + prefix + "_1.fastq.gz trim_" + prefix + "_2.fastq.gz >" + prefix + ".fasta 2>" + prefix + "_skesa.log &&")
			ns.write(" quast.py --threads " + str(threadsper) + " -o " + prefix + " --min-contig 1 -l \"" + prefix + "\" --contig-thresholds 500,1000 " + prefix + ".fasta && tail -n 1 " + prefix + "/transposed_report.tsv >>assembly_stats_nohead.tsv")
	else:
		if args.spades == True:
			ns.write(" spades.py --pe1-1 " + fs + " --pe1-2 " + rs + " -o ./" + prefix + "/ --careful -t " + str(threadsper) + " -m " + str(memper) + " &&")
			ns.write(" quast.py --threads " + str(threadsper) + " -o " + prefix + " --min-contig 1 -l \"" + prefix + "\" --contig-thresholds 500,1000 " + prefix + "/contigs.fasta && tail -n 1 " + prefix + "/transposed_report.tsv >>assembly_stats_nohead.tsv &&")
			ns.write(" mv " + prefix + "/contigs.fasta ./" + prefix + ".fasta && mv " + prefix + "/spades.log ./" + prefix + "_spades.log")
		else:
			ns.write(" skesa --cores " + str(threadsper) + " --memory " + str(memper) + " --fastq " + fs + " " + rs + " >" + prefix + ".fasta 2>" + prefix + "_skesa.log &&")
			ns.write(" quast.py --threads " + str(threadsper) + " -o " + prefix + " --min-contig 1 -l \"" + prefix + "\" --contig-thresholds 500,1000 " + prefix + ".fasta && tail -n 1 " + prefix + "/transposed_report.tsv >>assembly_stats_nohead.tsv")
	if args.no_qc == False:
		ns.write(" && bbmap.sh -Xmx"+ str(memper) + "g nodisk ref=" + prefix + ".fasta in1=" + fs + " in2=" + rs + " covstats=" + prefix + "_covstats.txt t=" + str(threadsper) + " &>" + prefix + "_bbmap.log")
	if args.amr:
		ns.write(" && amrfinder --threads " +str(threadsper) + " -n " + prefix + ".fasta >" + prefix + "_amrfinder.tsv 2>" + prefix + "_amrfinder.log")
	ns.write(";\n")

ns.close()
 
#This step generates the bash script to set the environment and launch perl_fork_univ.pl, calling the previously generated input file 

fn2="./LaunchJobs_" + strT + ".sh"
ss=open(fn2, 'w')

ss.write("#!/bin/bash \n")
ss.write("\n")
ss.write("export PATH=/workdir/miniconda3/bin:$PATH\n\n")
ss.write("source activate bacWGS\n\n")
ss.write("\nperl_fork_univ.pl " + fn1 + " " + str(njobs) + " &>" + fn1 + ".log \n")
ss.write("(head -n 1 " + prefs[0] + "/transposed_report.tsv && cat assembly_stats_nohead.tsv) >assembly_stats.tsv && rm assembly_stats_nohead.tsv\n")
if args.trim == True:
	ss.write("\necho \"trimmomatic $(trimmomatic -version)\" >>software_versions_" + strT + ".txt\n")
	ss.write("mkdir Trimmomatic_output && mv trim* Trimmomatic_output\n")
if args.spades == True:
	ss.write("\nspades.py -v >>software_versions_" + strT + ".txt\n")
	ss.write("mkdir SPAdes_logs && mv *_spades.log SPAdes_logs\n")
else:
	ss.write("\nskesa -v >>software_versions_" + strT + ".txt\n")
	ss.write("mkdir SKESA_logs && mv *_skesa.log SKESA_logs\n")
ss.write("\nquast.py -v >>software_versions_" + strT + ".txt\n")
if args.no_qc == False:
	ss.write("\nbbmap.sh --version 2>&1 >/dev/null| sed -n -e 's/^.*BBMap version /BBMap /p' >>software_versions_" + strT + ".txt\n")
	ss.write("echo -e \"Isolate\tTotal Length\tContigs\tN50\tMapped Coverage\tInsert Size\" >AssemblyQC_" + strT + ".tsv\n")
	ss.write("for s in " + ' '.join(prefs) + "; do i=$(grep $s assembly_stats.tsv |awk '{print $1,\"\t\",$8,\"\t\",$6,\"\t\",$10}');\n j=$(grep \"Average coverage\" ${s}_bbmap.log |cut -d$'\t' -f 2);\n k=$(grep \"insert size avg\" ${s}_bbmap.log |cut -d$'\t' -f 2);\n echo -e \"$i\t$j\t$k\";done >>AssemblyQC_" + strT + ".tsv\n")
	ss.write("mkdir BBMap_output && mv *_bbmap.log BBMap_output && mv *_covstats.txt BBMap_output\n")
if args.amr:
	ss.write("\nconda list|grep \"amrfinder\"|tr -s ' ' |cut -d \" \" -f 1,2 >>software_versions_" + strT + ".txt\n")
	ss.write("echo \"AMRFinder database$(ls -l $(which amrfinder|sed 's/amrfinder/data/')|grep \"latest\"|cut -d \">\" -f 2)\" >>software_versions_" + strT + ".txt\n")
	amrfiles = [p + "_amrfinder.tsv" for p in prefs]
	ss.write("make_amr_table.py " +  ' '.join(amrfiles) + " >ARG_table_" + strT + ".tsv\n")
	ss.write("mkdir AMRFinder_output && mv *amrfinder* AMRFinder_output\n")
if args.salm_sero:
	ss.write("\nsistr -V >>software_versions_" + strT + ".txt\n")
	fastas = [p + ".fasta" for p in prefs]
	ss.write("sistr --qc -f tab -t " + str(ncpus) + " -o Salm_seros_sistr_" + strT + " " + ' '.join(fastas) + " 2>Salm_seros_sistr_" + strT + ".log\n")
	ss.write("cut -d$'\t' -f 7,14,15 Salm_seros_sistr_" + strT + ".tab >SALM_sero_" + strT + ".tsv\n")
	ss.write("mkdir SISTR_output && mv Salm_seros_sistr_* SISTR_output\n")
if args.ecol_sero:
	ss.write("\nectyper -V >>software_versions_" + strT + ".txt\n")
	fastas = [p + ".fasta" for p in prefs]
	ss.write("ectyper -i " + ','.join(fastas) + " --cores " + str(ncpus) + " --verify -o ectyper_out_" + strT + " &>ectyper_" + strT + ".log\n")
	ss.write("cut -d$'\t' -f 1,2,3 ectyper_out_" + strT + "/output.tsv >ECOL_sero_" + strT + ".tsv\n")
	ss.write("mv ectyper_" + strT + ".log ectyper_out_" + strT + "\n")
ss.write("\nconda deactivate")
ss.close()
call(["chmod", "u+x", fn2])

#Optionally run the job
if args.launch:
	call([fn2])
else:
	print("To run, call " + fn2)