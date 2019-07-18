#!/workdir/miniconda3/envs/bacWGS/bin/python

import sys
import csv

'''
amrdict={}
for amr_row in csv.reader(open(sys.argv[1]), csv.excel_tab):
	amrdict[amr_row[0]] = amr_row[4:6]
'''	
amrout = csv.writer(sys.stdout, csv.excel_tab)
amrout.writerow(["Isolate",	"Gene",	"Identity (%)",	"Coverage (%)",	"Class", "Subclass"])

for file in sys.argv[1:len(sys.argv)]:
	linecount=0
	for line in csv.reader(open(file), csv.excel_tab):
		linecount += 1
		outrow = [file.split(sep="_amrfinder")[0]]
		if linecount != 1:
			outrow = outrow + [line[5]] + [line[16]] + [line[15]] + [line[10]] + [line[11]]
			amrout.writerow(outrow)
	if linecount == 1:
		outrow = outrow + ["None", "NA", "NA", "NA", "NA"]
		amrout.writerow(outrow)