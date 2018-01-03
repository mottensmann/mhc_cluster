`Clustering MHC class II sequences`

_This repository was forked from https://github.com/nextgenusfs/mhc_cluster.git` and adapted for processing Fur seal MHC II sequences obtained from Illumina MiSeq reads_

### Using cluster_mhc2.py

The script **cluster_mhc2.py** allows to cluster high-quality Illumina reads into putative alleles using the [Unoise3](http://drive5.com/usearch/manual/cmd_unoise3.html) approach [(Edgar 2016)](https://www.biorxiv.org/content/early/2016/10/15/081257).
The expected input of this pipeline is a single _fastq_ file that may contain barcodes for individual samples within the header as shown below. Barcodes are required for demultiplexing reads prior to mapping reads to the generated list of [`Zotus'](https://drive5.com/usearch/manual/pipe_otus.html), i.e. alleles.

```
@MISEQ:279:000000000-AVVMJ:1:1101:14590:18831:N:0:barcodelabel=CAGAGAGGAAGGAGTA
```
For further information and a description of parameters invoke the helpfile:

```bash
cluster_mhc2.py -h
```

### Prerequisites

The subfolder **lib** contains required functions and may contain hidden Markov models created from multiple sequence alignment of previously characterised MHC genes using [HMMER3](http://hmmer.org/).   
Example using [MUSCLE](http://www.drive5.com/muscle/manual/) and [HMMER3](hmmer.org) to create these files:

```bash
##	Align sequences
muscle -in input_sequences.fasta -out aligned_sequences.afa
#	Create hidden markov model	
hmmbuild hmm aligned_sequences.afa
#	Create auxiliary files
hmmpress hmm
```

### Dependencies

* [USEARCH10](http://www.drive5.com/usearch)
* [VSEARCH](https://github.com/torognes/vsearch)

```bash
wget https://github.com/torognes/vsearch/releases/download/v2.4.4/vsearch-2.4.4-linux-x86_64.tar.gz
tar xzf vsearch-2.4.4-linux-x86_64.tar.gz
```

* [HMMER3 v3.1b2](hmmer.org)

```bash
tar zxf hmmer-3.1b2.tar.gz
cd hmmer-3.1b2
./configure
make
make check
```

* [Biophyton](http://biopython.org/wiki/Download)

```bash
pip install biopython
```
*Installing Biophyton can easily cause some problems. See troubleshooting options here: https://askubuntu.com/questions/677566/biopython-installation* 
 
* [natsort](https://pypi.python.org/pypi/natsort)

```bash
pip install natsort
```

### References

Edgar, R.C., 2016. UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing. bioRxiv, p.081257.
