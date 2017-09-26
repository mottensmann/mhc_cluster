**Clustering mhc class II sequences**

*This repository was forked from https://github.com/nextgenusfs/mhc_cluster.git` and subsequently modified. See details there.*

The script `cluster_mhc.py` allows to cluster mhc class II sequences into `'OTUs'` that represent putatitive alleles. Input sequences are `FASTQ` files generated on a Illumina MiSeq run. The following steps may be executed outside the clustering:

* Extracting barcodes from sequences (`QIIME extract_barcodes.py`)
* Merging paired-end reads (e.g. `vsearch ---fastq_mergepairs`)
* Stripping primer sequences and extracting full length exon sequences (e.g. `fastx_trimmer`)
* Adding barcodes to sequence headers (e.g. `custom R scripts`)

### Using `cluster_mhc.py`

```
## helpfile
cluster_mhc.py -h
## usage
cluster_mhc.py -f reads.fastq -o output_name
```

The subfolder `lib` contains required functions and may contain hidden marko models created from multiple sequence alignment of published mhc alleles using `HMMER3` for both nucleotide and amino acid sequences.  

### Dependent software

The installation of some of these may cause some trouble and requires different attempts depending on system and the rights of the user (*root* etc...)

* USEARCH8 (http://www.drive5.com/usearch)
* VSEARCH (https://github.com/torognes/vsearch)
```
wget https://github.com/torognes/vsearch/releases/download/v2.4.4/vsearch-2.4.4-linux-x86_64.tar.gz
tar xzf vsearch-2.4.4-linux-x86_64.tar.gz
```

* HMMER3 v3.1b2 (hmmer.org)
Compiling from source
```
tar zxf hmmer-3.1b2.tar.gz
cd hmmer-3.1b2
./configure
make
make check
```

* Biophyton (http://biopython.org/wiki/Download)
```
pip install biopython
```
 *Potential troubleshooting option: https://askubuntu.com/questions/677566/biopython-installation* 
 
* natsort
```
pip install natsort
```
