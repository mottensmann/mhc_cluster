#!/usr/bin/env python

# Clustering Illumina MiSeq reads of MHC sequences into 'OTUs' i.e. putative alleles
# ===================================================================================
#
# This scripts does the following main steps:
# 1.) Filter reads against an hidden markov model comprising mhc a multiple sequence alignment
#	of known mhc sequcens from various mammals. (Computationally demanding!)
# 2.) Quality filter reads prior to cluster generation based on expected error rates
# 3.) Find unique sequences and estimate their abundance for the filtered reads
# 4.) Cluster filtered reads into OTUs 
# 5.) Map all reads to the OTUs

# Perequisites are:
# 	A hidden markov model may be constructed based on a multiple sequence alignment and stored in the 
# 	subfolder 'lib'.
#	Example using MUSCLE and HMMER3:
#	Align sequences
# 	muscle -in <fasta file containing sequences to align> -out <output fasta file with ending .afa>
#	Create hidden markov model	
#	hmmbuild <output> <muscle alignment>
#	Create help files
#	hmmpress <hmmbuild output>
#
# Principal outline adapted from the script 'mhc-OTU_cluster.py':
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

parser=argparse.ArgumentParser(prog='cluster_mhc.py', usage="%(prog)s [options] -f file.demux.fq\n%(prog)s -h for help menu",
    description='''Clustering mhc sequences into OTUs based on hidden markov model references.''',
    epilog="""Meinolf Ottensmann (2017) https://github.com/mottensmann/mhc_cluster""",
    formatter_class=MyFormatter)

parser.add_argument('-f','--fastq', dest="FASTQ", required=True, help='FASTQ input file')
parser.add_argument('-o','--out', default='out', help='Path and name of the output')
parser.add_argument('-maxee','--maxee', default='1.0', help='Expected error rate for filtering')
parser.add_argument('-l','--length', default='auto', help='Trim Length')
parser.add_argument('-pct_otu','--pct_otu', default='1', help="Minimum differeces between OTUs")
parser.add_argument('-pct_mapping','--pct_mapping', default='0.97', help="Identity threshold for mapping reads to OTUs")
parser.add_argument('-num_diff','--num_diff', default='False', help="OTU Clustering Number of differences")
parser.add_argument('-minsize','--minsize', default='2', help='Minimum abundance to keep for clustering')
parser.add_argument('-minfreq','--minfreq', default='False', help='Minimum frequency to keep for clustering')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8.exe', help='usearch version to use')
parser.add_argument('-hmm','--hmm', dest="hmm", default='mammalia_dqb.hmm', help='Hidden markov model reference')
parser.add_argument('-hmm_aas','--hmm_aas', dest="hmm_aas", default='mhcII_beta_aas.hmm', help='Hidden markov model reference of amino acid sequences')
parser.add_argument('--run_hmm', action='store_true', help='Filter raw reads against Hidden markov model')
parser.add_argument('-cpus','--cpus', default=4, help='Number of cpus')
parser.add_argument('--translate', action='store_true', help='Translate OTUs')
args=parser.parse_args()

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
log_name = args.out + '.EE' + args.maxee + '.log'
if os.path.isfile(log_name):
    os.remove(log_name)
log_file = open(log_name, 'ab')

usearch = args.usearch

try:
    subprocess.call([usearch, '--version'], stdout = log_file, stderr = log_file)
except OSError:
    print "%s not found in your PATH, exiting." % usearch 
    os._exit(1)

if args.num_diff != 'False':
    print "\nWarning:  --num_diff = %s was specified, this will override the --pct_otu option" % args.num_diff
	
if args.minfreq != 'False':
    print "\nWarning:  --minfreq = %s was specified, this will override the --minsize option" % args.minfreq

	
	
# set up log file for prints to the console
console_log = args.out + '.EE' + args.maxee + '.console_log.txt'

file = open(console_log, "w")
file.write('Arguments:\n' + str(args))
file.write('\nLoading records: ' + '{0:,}'.format(countfastq(args.FASTQ)) + ' reads\n') 
file.close() 

# Count input reads
print '\nLoading records: ' + '{0:,}'.format(countfastq(args.FASTQ)) + ' reads'

# Convert to fasta to use HMMER3

raw_reads_fasta = args.out + '.EE' + args.maxee + '.raw_reads.fa'
subprocess.call([usearch, '-fastq_filter', args.FASTQ, '-fastq_qmax', '45', '-fastaout', raw_reads_fasta], stdout = log_file, stderr = log_file)

# 1. Run HMMER3 to filter contaminant sequences out.
# ##################################################
FNULL = open(os.devnull, 'w')
all_reads_hmm_out = args.out + '.EE' + args.maxee + '.raw_reads.hmm.txt'
if args.run_hmm:
	hmm = script_path + '/lib/' + args.hmm 
	print "CMD: Filter all reads against HMM \nhmmscan --cpu %s --domtblout %s %s %s\n" % (cpus, all_reads_hmm_out, hmm, raw_reads_fasta)
	file = open(console_log, "a")
	file.write("\nCMD: Filter all reads against HMM \nhmmscan --cpu %s --domtblout %s %s %s\n" % (cpus, all_reads_hmm_out, hmm, raw_reads_fasta)) 
	file.close()
 
	subprocess.call(['hmmscan', '--cpu', cpus, '--domtblout', all_reads_hmm_out, hmm, raw_reads_fasta], stdout = FNULL, stderr = FNULL)

# Parse HMMer results
	file = open(console_log, "a")
	file.write("\nCMD: Filtering HMM results\n") 
	file.close()

	print "CMD: Filtering HMM results"
	q_list = []
	l_list = []
	counter = 0
	hmmer_raw_results = open(all_reads_hmm_out, 'r')
	for qresult in SearchIO.parse(hmmer_raw_results, "hmmscan3-domtab"):
		q_list.append(qresult.id)
		l_list.append(qresult.seq_len)
		counter += 1
	hmmer_raw_results.close()

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

	file = open(console_log, "a")
	file.write("\nCMD: Retrieving results and trimming/padding to %i bp for clustering" % trim_len) 
	file.close()
    
	print "\nCMD: Retrieving results and trimming/padding to %i bp for clustering" % trim_len
# Retrieve sequences in the "pass" list by looping through the query list, make same length
	raw_pass_out = args.out + '.EE' + args.maxee + '.raw_reads.hmm.pass.fa'
	pass_handle = open(raw_pass_out, 'wb')
	filtered = SeqIO.parse(raw_reads_fasta, "fasta")
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

	file = open(console_log, "a")
	file.write("\nOutput: " + '{0:,}'.format(countfasta(raw_pass_out)) + ' reads\n') 
	file.close()

	print "Output: " + '{0:,}'.format(countfasta(raw_pass_out)) + ' reads\n'

# 2.) Quality filter reads based on expected errors before clustering OTUs
# ########################################################################
    
# Run vsearch fastq filtering step, output to fasta
filter_out = args.out + '.EE' + args.maxee + '.filter.fa'
print "\nCMD: Quality Filtering\n%s -fastq_filter %s -fastq_maxee %s -fastq_qmax 45 -fastaout %s" % (usearch, args.FASTQ, args.maxee, filter_out)

file = open(console_log, "a")
file.write("\nCMD: Quality Filtering\n%s -fastq_filter %s -fastq_maxee %s -fastq_qmax 45 -fastaout %s" % (usearch, args.FASTQ, args.maxee, filter_out)) 
file.close() 

subprocess.call([usearch, '-fastq_filter', args.FASTQ, '-fastq_maxee', args.maxee, '-fastq_qmax', '45', '-fastaout', filter_out], stdout = log_file, stderr = log_file)

print "Output: " + '{0:,}'.format(countfasta(filter_out)) + ' reads\n'

file = open(console_log, "a")
file.write("\nOutput: " + '{0:,}'.format(countfasta(filter_out)) + ' reads\n') 
file.close() 

# Now, filter against hmm again

hmm_filtered_out = args.out + '.EE' + args.maxee + '.reads.filtered.hmm.txt'
hmm = script_path + '/lib/' + args.hmm 

print "CMD: Filter filtered reads against HMM \nhmmscan --cpu %s --domtblout %s %s %s\n" % (cpus, hmm_filtered_out, hmm, filter_out)

file = open(console_log, "a")
file.write("\nCMD: Filter filtered reads against HMM \nhmmscan --cpu %s --domtblout %s %s %s\n" % (cpus, hmm_filtered_out, hmm, filter_out)) 
file.close()

subprocess.call(['hmmscan', '--cpu', cpus, '--domtblout', hmm_filtered_out, hmm, filter_out], stdout = FNULL, stderr = FNULL)

# Parse HMMer results
file = open(console_log, "a")
file.write("\nCMD: Filtering HMM results\n") 
file.close()

print "CMD: Filtering HMM results"
q_list = []
l_list = []
counter = 0
hmmer_results = open(hmm_filtered_out, 'r')
for qresult in SearchIO.parse(hmmer_results, "hmmscan3-domtab"):
    q_list.append(qresult.id)
    l_list.append(qresult.seq_len)
    counter += 1
hmmer_results.close()
sum_2 = sum(l_list)
avg_len = sum_2 / counter

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

file = open(console_log, "a")
file.write("\nCMD: Retrieving results and trimming/padding to %i bp for clustering" % trim_len) 
file.close()
    
print "\nCMD: Retrieving results and trimming/padding to %i bp for clustering" % trim_len
# Retrieve sequences in the "pass" list by looping through the query list, make same length
filtered_pass_out = args.out + '.EE' + args.maxee + '.hmm.filtered.pass.fa'
pass_handle = open(filtered_pass_out, 'wb')
filtered = SeqIO.parse(filter_out, "fasta")
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

file = open(console_log, "a")
file.write("\nOutput: " + '{0:,}'.format(countfasta(filtered_pass_out)) + ' reads\n') 
file.close()

print "Output: " + '{0:,}'.format(countfasta(filtered_pass_out)) + ' reads\n'

# 3.) De-replication of reads
# ###########################
derep_out = args.out + '.EE' + args.maxee + '.derep.fa'
sort_out = args.out + '.EE' + args.maxee + '.sort.fa'


file = open(console_log, "a")
file.write("\nCMD: De-replication\n%s -derep_fulllength %s -sizeout -fastaout %s" % (usearch, filter_out, derep_out)) 
file.close()

print "CMD: De-replication\n%s -derep_fulllength %s -sizeout -fastaout %s" % (usearch, filter_out, derep_out)
subprocess.call([usearch, '-derep_fulllength', filter_out, '-sizeout', '-fastaout', derep_out], stdout = log_file, stderr = log_file)
print "Output: " + '{0:,}'.format(countfasta(derep_out)) + ' reads\n'

if args.minfreq != 'False':
	number_of_reads = countfasta(derep_out)
	cutoff = math.floor(number_of_reads/100*float(args.minfreq))
	args.minsize = str(cutoff)
	
print "CMD: usearch -sortbysize -minsize %s"  % (args.minsize)

subprocess.call([usearch, '-sortbysize', derep_out, '-minsize', args.minsize, '-fastaout', sort_out], stdout = log_file, stderr = log_file)

print "Output: " + '{0:,}'.format(countfasta(sort_out)) + ' reads\n'

# 4.) Cluster into OTUs
# #####################

otu_out = args.out + '.EE' + args.maxee + '.otus.fa'
file = open(console_log, "a")
file.write("\nCMD: Clustering OTUs\n%s -cluster_otus %s -sizein -sizeout -relabel MHC_ -otus %s" % (usearch, sort_out, otu_out)) 
file.close()

print "CMD: Clustering OTUs\n%s -cluster_otus %s -sizein -sizeout -relabel MHC_ -otu_radius_pct %s -otus %s" % (usearch, sort_out, args.pct_otu, otu_out)
subprocess.call([usearch, '-cluster_otus', sort_out, '-sizein', '-relabel', 'MHC_', '-otu_radius_pct', args.pct_otu, '-otus', otu_out], stdout = log_file, stderr = log_file)

# Fix OTUs, remove trailing N's
fix_otus = args.out + '.EE' + args.maxee + '.fixed.otus.fa'
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

# 5.) Map reads back to OTUs
# ###########################

if args.run_hmm == False:
	raw_pass_out = raw_reads_fasta

uc_out = args.out + '.EE' + args.maxee + '.mapping.uc'
otu_table = args.out + '.EE' + args.maxee + '.otu_table.txt'

file = open(console_log, "a")
file.write("\nCMD: Mapping Reads to OTUs\n%s -otutab %s -strand plus -otutabout %s -mapout %s\n" % (usearch, raw_pass_out, otu_table, uc_out)) 
file.close()

print "CMD: Mapping Reads to OTUs\n%s -usearch_global %s -strand plus -id %s -db %s -uc %s\n" % (usearch, raw_pass_out, args.pct_mapping, fix_otus, uc_out)
subprocess.call([usearch, '-usearch_global', raw_pass_out, '-strand', 'plus', '-id', args.pct_mapping, '-db', fix_otus, '-uc', uc_out], stdout = log_file, stderr = log_file)

# Build OTU table
otu_table = args.out + '.EE' + args.maxee + '.otu_table.txt'
uc2tab = script_path + "/lib/uc2otutab.py"

file = open(console_log, "a")
file.write("\nCMD: Creating OTU Table\npython %s %s > %s" % (uc2tab, uc_out, otu_table)) 
file.close()

print "CMD: Creating OTU Table\npython %s %s > %s" % (uc2tab, uc_out, otu_table)
os.system('%s %s %s %s %s' % ('python', uc2tab, uc_out, '>', otu_table))

#translate to protein space (optional)
if args.translate:
    print "\nCMD: Translating to Amino Acid Sequence\n"
    trans_out = args.out + '.EE' + args.maxee + '.proteins.fa'
    trans_file = open(trans_out, "w")
    trans_in = open(fix_otus, 'r')
    for rec in SeqIO.parse(trans_in, "fasta"):
        trans_file.write(">%s;frame_1\n" % (rec.id))
        trans_file.write("%s\n" % (rec.seq.translate(to_stop=True)))
        trans_file.write(">%s;frame_2\n" % (rec.id))
        trans_file.write("%s\n" % (rec.seq[1:].translate(to_stop=True)))
        trans_file.write(">%s;frame_3\n" % (rec.id))
        trans_file.write("%s\n" % (rec.seq[2:].translate(to_stop=True)))
    trans_in.close()
    trans_file.close()
 
    #HMM against translated amino acids
    trans_hmm = args.out + '.EE' + args.maxee + '.proteins.hmm.txt'
    hmm_prot = script_path + '/lib/' + args.hmm_aas 
	
    print "\nCMD: Running HMMER3 MHC HMM model (using %s cpus)\nhmmscan --cpu %s --tblout %s %s %s" % (cpus, cpus, trans_hmm, hmm_prot, trans_out)
    subprocess.call(['hmmscan', '--cpu', cpus, '--domtblout', trans_hmm, hmm_prot, trans_out], stdout = FNULL, stderr = FNULL)
    
    #now filter results for best hit, and get alignment coordinates
    hmmer_prots = open(trans_hmm, 'r')
    hit_list = {}
    for qresult in SearchIO.parse(hmmer_prots, "hmmscan3-domtab"):
        hits = qresult.hits
        if len(hits) > 0:
            hit_list['%s' % hits[0].query_id] = []
            hit_list['%s' % hits[0].query_id].append(hits[0].hsps[0].query_start)
            hit_list['%s' % hits[0].query_id].append(hits[0].hsps[0].query_end)
    hmmer_prots.close()

    #now filter the fasta file to get only hits that have MHC domain
    pass_prot = args.out + '.EE' + args.maxee + '.proteins.pass.fa'
    pass_file = open(pass_prot, 'wb')
    SeqRecords = SeqIO.to_dict(SeqIO.parse(trans_out, "fasta"))
    for key in sorted(hit_list.keys(), key=natural_sort_key):
        start = hit_list[key][0]
        end = hit_list[key][1]
        subseq = SeqRecords[key][start:end].seq
        pass_file.write(">%s;%s-%s\n%s\n" % (key, start, end, subseq))
    pass_file.close()
    
    #finally concatenate duplicated sequences
    concat_out = args.out + '.EE' + args.maxee + '.proteins.unique.fa'
    keep = {}
    total_count = 0
    concat_file = open(concat_out, 'wb')
    Seq = SeqIO.parse(pass_prot, "fasta")
    for rec in Seq:
        total_count += 1
        sequence=str(rec.seq).upper()
        if sequence not in keep:
            keep[sequence]=rec.id
        else:
            keep[sequence]+="|"+rec.id
    flipKeep = {y:x for x,y in keep.iteritems()}
    unique_count = 0
    for key in sorted(flipKeep, key=natural_sort_key):
        unique_count += 1
        concat_file.write(">"+key+"\n"+flipKeep[key]+"\n")
    concat_file.close()
    print "%10u total proteins" % total_count
    print "%10u unique proteins" % unique_count

    
#output reads per Barcode for original, filtered, and hmm passed files
#now loop through data and find barcoded samples, counting each.....
BarcodeCountA = {}
with open(args.FASTQ, 'rU') as input:
    header = itertools.islice(input, 0, None, 4)
    for line in header:
        ID = line.split("=")[-1].split(";")[0]
        if ID not in BarcodeCountA:
            BarcodeCountA[ID] = 1
        else:
            BarcodeCountA[ID] += 1
BarcodeCountB = {}
with open(filter_out, 'rU') as input:
    for line in input:
        if line.startswith('>'):
            ID = line.split("=")[-1].split(";")[0]
            if ID not in BarcodeCountB:
                BarcodeCountB[ID] = 1
            else:
                BarcodeCountB[ID] += 1
BarcodeCountC = {}
with open(raw_pass_out, 'rU') as input:
    for line in input:
        if line.startswith('>'):
            ID = line.split("=")[-1].split(";")[0]
            if ID not in BarcodeCountC:
                BarcodeCountC[ID] = 1
            else:
                BarcodeCountC[ID] += 1
bc_count = args.out + '.EE' + args.maxee + '.barcode.counts.txt'
with open(bc_count, 'w') as output:
    output.write("Barcode\tDemux_total\tFilter_total\tHMM_total\n")               
    for k,v in natsorted(BarcodeCountA.items(), key=lambda (k,v): v, reverse=True):
        bc_name = str(k)
        allCount = str(BarcodeCountA[k])
        filtCount = BarcodeCountB.get(k)
        hmmCount = BarcodeCountC.get(k)
        output.write("%s\t%s\t%s\t%s\n" % (bc_name, allCount, filtCount, hmmCount))

#Print location of files to STDOUT
print "\n------------------------------------------------"
print "OTU Clustering Script has Finished Successfully"
print "------------------------------------------------"
print ("Input fastq format:	%s" % (args.FASTQ))
print ("Input fasta format:	%s" % (raw_reads_fasta))
print ("HMM passed raw reads:	%s" % (raw_pass_out))
print ("Filtered reads:	%s" % (filter_out))
print ("HMM passed filtered reads:	%s" % (filtered_pass_out))
print ("Dereplicated reads:	%s" % (derep_out))
print ("Reads sorted by frequency:	%s" % (sort_out))
print ("Clustered OTUs:	%s" % (fix_otus))
if args.translate:
    print ("Translated OTUs:	%s" % (pass_prot))
    print ("Translated Unique:	%s" % (concat_out))
print ("usearch Mapping file:	%s" % (uc_out))
print ("OTU Table:	%s" % (otu_table))
print ("LogFile:	%s" % (console_log))
print ("Reads per Barcode:	%s" % (bc_count))
print "---------------------------------------------------------------"

 
