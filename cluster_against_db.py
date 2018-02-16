#!/usr/bin/env python

# This script maps reads against a reference database containg allele variants. 
# For more details see script mhc_cluster2.py 

# Meinolf Ottensmann, 2017

import os, argparse, subprocess, inspect, re, multiprocessing, warnings, itertools, math
from natsort import natsorted
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO
    from Bio import SeqIO
    
# Get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='cluster_mhc2.py', usage="%(prog)s [options] -f file.demux.fq\n%(prog)s -h for help menu",
    description='''Clustering mhc sequences into OTUs based on hidden markov model references.''',
    epilog="""Meinolf Ottensmann (2017) https://github.com/mottensmann/mhc_cluster""",
    formatter_class=MyFormatter)

parser.add_argument('-f','--fastq', dest="FASTQ", required=True, help='FASTQ input file')
parser.add_argument('-o','--out', default='out', help='Path and prefix of the output')
parser.add_argument('-ref','--reference', required=True, help='Reference list of otus')
parser.add_argument('-pct','--pct_mapping', default='1.0', help="Identity threshold for mapping reads to OTUs")
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8.exe', help='usearch version to use')
parser.add_argument('-cpus','--cpus', default=4, help='Number of cpus')
args=parser.parse_args()

# make proper output name
args.out = args.out + '_pct_' + args.pct_mapping 

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count
    
def countfastq(input):
    lines = sum(1 for line in open(input))
    count = int(lines) / 4
    return count

# set number of cpus
cpus = args.cpus
cpus = str(cpus)

# Open log file for usearch8 stderr redirect
log_name = args.out + '.log'
if os.path.isfile(log_name):
    os.remove(log_name)
log_file = open(log_name, 'ab')

usearch = args.usearch

try:
    subprocess.call([usearch, '--version'], stdout = log_file, stderr = log_file)
except OSError:
    print "%s not found in your PATH, exiting." % usearch 
    os._exit(1)

# set up log file for prints to the console
console_log = args.out + '.console_log.txt'

file = open(console_log, "w")
file.write('Arguments:\n' + str(args))
file.write('\nLoading records: ' + '{0:,}'.format(countfastq(args.FASTQ)) + ' reads\n') 
file.close() 

# Count input reads
print '\nLoading records: ' + '{0:,}'.format(countfastq(args.FASTQ)) + ' reads'

# Convert to fasta to use HMMER3

raw_reads_fasta = args.out + '.raw_reads.fa'
subprocess.call([usearch, '-fastq_filter', args.FASTQ, '-fastq_qmax', '45', '-fastaout', raw_reads_fasta], stdout = log_file, stderr = log_file)

# to capture some output
FNULL = open(os.devnull, 'w')

# 6.) Map reads back to OTUs
# ###########################

uc_out = args.out + '.mapping.uc'
otu_table = args.out + '.otu_table.txt'

print "CMD: Mapping Reads to OTUs\n%s -usearch_global %s -strand plus -id %s -db %s -uc %s\n" % (usearch, raw_reads_fasta, args.pct_mapping, args.reference, uc_out)
subprocess.call([usearch, '-usearch_global', raw_reads_fasta, '-strand', 'plus', '-id', args.pct_mapping, '-db', args.reference, '-uc', uc_out], stdout = log_file, stderr = log_file)

# 7.) Build OTU table
# ###########################
otu_table = args.out + '.otu_table.txt'
uc2tab = script_path + "/lib/uc2otutab.py"

file = open(console_log, "a")
file.write("\nCMD: Creating OTU Table\npython %s %s > %s" % (uc2tab, uc_out, otu_table)) 
file.close()

print "CMD: Creating OTU Table\npython %s %s > %s" % (uc2tab, uc_out, otu_table)
os.system('%s %s %s %s %s' % ('python', uc2tab, uc_out, '>', otu_table))

# 7.) Count Barcodes
# ###########################

## Fake counts, only to avoid breaking code in downstream analysis

BarcodeCountA = {}
with open(args.FASTQ, 'rU') as input:
    header = itertools.islice(input, 0, None, 4)
    for line in header:
        ID = line.split("=")[-1].split(";")[0]
        if ID not in BarcodeCountA:
            BarcodeCountA[ID] = 1
        else:
            BarcodeCountA[ID] += 1
		
bc_count = args.out + '.barcode.counts.txt'
with open(bc_count, 'w') as output:
    output.write("BarcodeSequence\tRaw_total\tFiltered_total\n")               
    for k,v in natsorted(BarcodeCountA.items(), key=lambda (k,v): v, reverse=True):
        bc_name = str(k)
        allCount = str(BarcodeCountA[k])
        filtCount = str(BarcodeCountA[k])
        output.write("%s\t%s\t%s\n" % (bc_name, allCount, filtCount))

#Print location of files to STDOUT
print "\n------------------------------------------------"
print "OTU Clustering Script has Finished Successfully"
print "------------------------------------------------"
print ("Input fastq format:	%s" % (args.FASTQ))
print ("Input fasta format:	%s" % (raw_reads_fasta))
print ("usearch Mapping file:	%s" % (uc_out))
print ("OTU Table:	%s" % (otu_table))
print ("LogFile:	%s" % (console_log))
print ("Reads per Barcode:	%s" % (bc_count))
print "---------------------------------------------------------------"
