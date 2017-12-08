#!/usr/bin/env python

# Processing Illumina MiSeq reads of MHC II loci
# ==============================================
#
# The main steps of this pipeline are:
# 
# (A) Clustering reads into putatitive alleles:
#   (i) Filtering reads based on expected error < 1 
#   (ii) De-replication of reads, by default discard reads with abundance < 10
#   (iii) Discard reads that are not homologous to the targeted MHC locus (Hidden-Markow-Model)
#   (iV) Cluster remaining reads into Zotus using Unoise3
#   
# (B) Demultiplex seuences using barcodes provided in the header
# 
# (C) Map raw reads to Zotus. By default matches require identical sequences 
# 
# (D) Export results, comprising among other things
#   (i) Zotu table 
#   (ii) list of alleles
#
# Adapted from the script 'mhc-OTU_cluster.py':
# Palmer et al (2016). PeerJ PrePrints 4 (2016): e1662v1
# https://github.com/nextgenusfs/mhc_cluster
# 
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
    description='''Processing Illumina MiSeq reads of MHC II loci.''',
    epilog="""Meinolf Ottensmann (2017) https://github.com/mottensmann/mhc_cluster""",
    formatter_class=MyFormatter)

parser.add_argument('-f','--fastq', dest="FASTQ", required=True, help='FASTQ input file')
parser.add_argument('-o','--out', default='out', help='Path and prefix of the output')
parser.add_argument('-l','--length', default='auto', help='Trim Length')
parser.add_argument('-pct','--pct_mapping', default='1.0', help="Identity threshold for mapping reads to OTUs")
parser.add_argument('-minsize','--minsize', default='10', help='Minimum abundance to keep for clustering')
parser.add_argument('-minfreq','--minfreq', default='False', help='Minimum frequency to keep for clustering')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch10.exe', help='usearch version to use')
parser.add_argument('-alpha', '--alpha', default = '2.0', help='unoise3 alpha parameter')
parser.add_argument('-e', '--error', default = '1.0', help='Expected error used to filter reads for clustering') 
parser.add_argument('-hmm','--hmm', dest="hmm", default='seal_dqb.hmm', help='Hidden markov model reference')
parser.add_argument('-cpus','--cpus', default=1, help='Number of cpus to use')
args=parser.parse_args()

# make proper output name
args.out = args.out + '_pct_' + args.pct_mapping + '_a_' + args.alpha + '_ee_' + args.error

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

# Open log file for usearch stderr redirect
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

if args.minfreq != 'False':
    print "\nWarning:  --minfreq = %s was specified, this will override the --minsize option" % args.minfreq

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


# 1.) Quality filter reads based on expected errors before clustering OTUs
# ########################################################################
    
# Run vsearch fastq filtering step, output to fasta
filter_out = args.out + '.filter.fa'
print "\nCMD: Quality Filtering\n%s -fastq_filter %s -fastq_maxee %s -fastq_qmax 45 -fastaout %s" % (usearch, args.FASTQ, args.error, filter_out)

file = open(console_log, "a")
file.write("\nCMD: Quality Filtering\n%s -fastq_filter %s -fastq_maxee %s -fastq_qmax 45 -fastaout %s" % (usearch, args.FASTQ, args.error, filter_out)) 
file.close() 

subprocess.call([usearch, '-fastq_filter', args.FASTQ, '-fastq_maxee', args.error, '-fastq_qmax', '45', '-fastaout', filter_out], stdout = log_file, stderr = log_file)

print "Output: " + '{0:,}'.format(countfasta(filter_out)) + ' reads\n'

file = open(console_log, "a")
file.write("\nOutput: " + '{0:,}'.format(countfasta(filter_out)) + ' reads\n') 
file.close() 


# 2.) De-replication of reads
# ###########################
derep_out = args.out + '.derep.fa'

file = open(console_log, "a")
file.write("\nCMD: De-replication\n%s --fastx_uniques %s -sizeout -fastaout %s -minuniquesize 1" % (usearch, filter_out, derep_out)) 
file.close()

print "CMD: De-replication\n%s --fastx_uniques %s -sizeout -fastaout %s" % (usearch, filter_out, derep_out)
subprocess.call([usearch, '--fastx_uniques', filter_out, '-sizeout', '-fastaout', derep_out, '-minuniquesize', '1'], stdout = log_file, stderr = log_file)
print "Output: " + '{0:,}'.format(countfasta(derep_out)) + ' reads\n'

# 3.) Filter against hidden markov model
# ######################################

derep_out_hmm = args.out + '.derep.hmm.fa'
hmm = script_path + '/lib/' + args.hmm 

print "CMD: Filter de-replicated reads against HMM \nhmmscan --cpu %s --domtblout %s %s %s\n" % (cpus, derep_out_hmm, hmm, derep_out)

file = open(console_log, "a")
file.write("\nCMD: Filter de-replicated reads against HMM \nhmmscan --cpu %s --domtblout %s %s %s\n" % (cpus, derep_out_hmm, hmm, derep_out)) 
file.close()
 
subprocess.call(['hmmscan', '--cpu', cpus, '--domtblout', derep_out_hmm, hmm, derep_out], stdout = FNULL, stderr = FNULL)

q_list = []
l_list = []
counter = 0
hmmer_results = open(derep_out_hmm, 'r')
for qresult in SearchIO.parse(hmmer_results, "hmmscan3-domtab"):
	q_list.append(qresult.id)
	l_list.append(qresult.seq_len)
	counter += 1
hmmer_results.close()

sum_1 = sum(l_list)
avg_len = sum_1 / counter

if args.length == 'auto':
	trim_len = avg_len
else:
	trim_len = int(args.length)

file = open(console_log, "a")
file.write("%10u sequences passed filter" % counter) 
file.write("%10u bp is average length" % avg_len) 
file.close()

print "%10u sequences passed filter" % counter
print "%10u bp is average length" % avg_len

pass_out = args.out + '.derep.hmm.pass.fa'
pass_handle = open(pass_out, 'wb')
filtered = SeqIO.parse(derep_out, "fasta")
for rec in filtered:
		if rec.id in q_list:
			L = len(rec.seq)
			if L < trim_len:
				Seq = rec.seq + (trim_len - L)*'N'
			else:
				T = trim_len - 1
				Seq = rec.seq[:T]
			pass_handle.write(">%s\n%s\n" % (rec.id, Seq))
pass_handle.close()


# 4.) Cluster into OTUs
# #####################

otu_out = args.out + '.otus.fa'
otu_sizes = args.out + '.otus_size.fa'

file = open(console_log, "a")
file.write("\nCMD: Clustering OTUs\n%s -'-unoise3' %s -sizein -sizeout -otus %s -unoise_alpha %s -minsize %s" % (usearch, pass_out, otu_out, args.alpha, args.minsize)) 
file.close()

if args.minfreq != 'False':
	number_of_reads = countfasta(pass_out)
	cutoff = math.floor(number_of_reads/100*float(args.minfreq))
	args.minsize = str(cutoff)
	print "minfreq %s equls minsize %s" %(args.minfreq, args.minsize)


print "CMD: Clustering OTUs\n%s -'-unoise3' %s -sizein -sizeout -otus %s -unoise_alpha %s -minsize %s" % (usearch, pass_out, otu_out, args.alpha, args.minsize)
subprocess.call([usearch, '-unoise3', pass_out, '-zotus', otu_out, '-unoise_alpha', args.alpha,'-minsize', args.minsize, ], stdout = log_file, stderr = log_file)

# Fix OTUs, remove trailing N's
fix_otus = args.out + '.fixed.otus.fa'
fix_handle = open(fix_otus, 'wb')
fix = open(otu_out, 'rb')
otu_count = 0
for rec in SeqIO.parse(fix, "fasta"):
    otu_count += 1
    Seq = re.sub('[^GATC]', "", str(rec.seq).upper())
    fix_handle.write(">%s\n%s\n" % (rec.id, Seq))
fix_handle.close()
fix.close()

file = open(console_log, "a")
file.write("%10u total OTUs\n" % otu_count) 
file.close()

print "%10u total OTUs\n" % otu_count

# 5.) Add size annotation to OTUs
# ###############################

subprocess.call([usearch, '-otutab', derep_out_hmm, '-db', otu_out, '-dbmatched', otu_sizes,'-sizeout'], stdout = log_file, stderr = log_file)

#-otutabout otutab.txt 


# 6.) Map reads back to OTUs
# ###########################

uc_out = args.out + '.mapping.uc'
otu_table = args.out + '.otu_table.txt'

file = open(console_log, "a")
file.write("\nCMD: Mapping Reads to OTUs\n%s -otutab %s -strand plus -otutabout %s -mapout %s\n" % (usearch, raw_reads_fasta, otu_table, uc_out)) 
file.close()

print "CMD: Mapping Reads to OTUs\n%s -usearch_global %s -strand plus -id %s -db %s -uc %s\n" % (usearch, raw_reads_fasta, args.pct_mapping, fix_otus, uc_out)
subprocess.call([usearch, '-usearch_global', raw_reads_fasta, '-strand', 'plus', '-id', args.pct_mapping, '-db', fix_otus, '-uc', uc_out], stdout = log_file, stderr = log_file)

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

# raw
BarcodeCountA = {}
with open(args.FASTQ, 'rU') as input:
    header = itertools.islice(input, 0, None, 4)
    for line in header:
        ID = line.split("=")[-1].split(";")[0]
        if ID not in BarcodeCountA:
            BarcodeCountA[ID] = 1
        else:
            BarcodeCountA[ID] += 1
# filtered
BarcodeCountB = {}
with open(filter_out, 'rU') as input:
    for line in input:
        if line.startswith('>'):
            ID = line.split("=")[-1].split(";")[0]
            if ID not in BarcodeCountB:
                BarcodeCountB[ID] = 1
            else:
                BarcodeCountB[ID] += 1

bc_count = args.out + '.barcode.counts.txt'
with open(bc_count, 'w') as output:
    output.write("BarcodeSequence\tRaw_total\tFiltered_total\n")               
    for k,v in natsorted(BarcodeCountA.items(), key=lambda (k,v): v, reverse=True):
        bc_name = str(k)
        allCount = str(BarcodeCountA[k])
        filtCount = BarcodeCountB.get(k)
        output.write("%s\t%s\t%s\n" % (bc_name, allCount, filtCount))

#Print location of files to STDOUT
print "\n------------------------------------------------"
print "OTU Clustering Script has Finished Successfully"
print "------------------------------------------------"
print ("Input fastq format:	%s" % (args.FASTQ))
print ("Input fasta format:	%s" % (raw_reads_fasta))
print ("Filtered reads:	%s" % (filter_out))
print ("Dereplicated reads:	%s" % (derep_out))
print ("Dereplicated & hmm filtered:	%s" % (pass_out))
print ("Clustered OTUs:	%s" % (fix_otus))
print ("usearch Mapping file:	%s" % (uc_out))
print ("OTU Table:	%s" % (otu_table))
print ("LogFile:	%s" % (console_log))
print ("Reads per Barcode:	%s" % (bc_count))
print "---------------------------------------------------------------"
