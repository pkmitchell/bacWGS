#!/workdir/miniconda3/envs/bacWGS/bin/python

import sys
from subprocess import call
import argparse
import re
import multiprocessing as mp

'''
Run fastqc on a input samples, then parse output to produce table
'''

#Argument Parser
parser = argparse.ArgumentParser(description="Run FastQC and parse output to produce table of read quality metrics")
parser.add_argument("-p", "--threads", help="Maximum number of processors to use (Default = 1)", type = int, default=1)
parser.add_argument("Fastqs", help="Input fastqs to be processed", nargs='+')
args = parser.parse_args()


#Determining the number of jobs and threads per job
ncpus = args.threads #The total number of CPUs to be used simultaneously
totjobs = len(args.Fastqs) #The total number of input sequences, which is also the total number of jobs to run

#Determining the number of threads using the bounds provided by the user
tot_cpus = mp.cpu_count()

if ncpus > tot_cpus:
	ncpus = tot_cpus

if ncpus > totjobs: 
	ncpus = totjobs

#Call FastQC
fqc_prog = open("fastqc_progress.log", 'w')
fqc_cmd = ["fastqc", "--extract", "--threads", str(ncpus)] + args.Fastqs
call(fqc_cmd, stdout=fqc_prog, stderr=fqc_prog)

#Process FastQC output and write QC table
sys.stdout.write("Isolate\tDirection\tReads\tMean Length\tMean Q\tQ30%\t Est. Coverage\n")

fqc_files = [f.split(sep=".")[0]+"_fastqc/fastqc_data.txt" for f in args.Fastqs]

dirdict = {"1":"Forward", "2":"Reverse"} #silly little dictionary to say 1 means forward and 2 means reverse

for file in fqc_files:
	section=""
	qreadsum = 0.0
	qscoresum = 0.0
	q30readsum = 0.0
	lreadsum = 0.0
	lensum = 0.0
	stem=""
	direction=""

	with open(file, 'r') as f:
		for line in f:
			if line[0:2] == ">>":
				section = re.split("	", line[2:])[0]
			elif section not in ["", "Per sequence quality scores", "Sequence Length Distribution", "Basic Statistics"]:
				pass
			elif section == "Per sequence quality scores" and line[0] != "#":
				comps = re.split("	", line)
				qreadsum += float(comps[1])
				qscoresum += float(comps[0]) * float(comps[1])
				if float(comps[0]) >= 30:
					q30readsum += float(comps[1])
			elif section == "Sequence Length Distribution" and line[0] != "#":
				comps = re.split("	", line)
				lreadsum += float(comps[1])
				rng = list(map(float, comps[0].split(sep="-")))
				mp = sum(rng) / len(rng)
				lensum += mp * float(comps[1])
			elif section == "Basic Statistics" and line[0] != "#":
				comps = re.split("	", line)
				if comps[0] == "Filename":
					fname = comps[1].split(sep=".")[0]
					mtch= re.search(r'_R[12]_001$',fname)
					if mtch:
						name_comps=fname.split(sep="_")
						stem=name_comps[0]
						direction=dirdict[name_comps[3][1]]
					else:
						mtch2 = re.search(r'_[12]$',fname)
						if mtch2:
							name_comps=fname.split(sep="_")
							stem=name_comps[0]
							direction=dirdict[name_comps[1][0]]
				elif comps[0] == "Total Sequences":
					totseq = float(comps[1])
	if len(set([qreadsum, lreadsum, totseq])) == 1:
		nreads = totseq
	else:
		nreads == "?"
	alen = lensum/lreadsum
	avgQ = qscoresum/qreadsum
	q30 = q30readsum/qreadsum
	sys.stdout.write(stem + "\t" + direction + "\t" + str(nreads) + "\t" + str(alen) + "\t" + str(avgQ) + "\t" + str(q30) + "\n")
	
sys.stdout.write("\n")
call(["fastqc", "--version"], stdout=sys.stdout)


