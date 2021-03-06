
`Clustering MHC class II alleles`
=================================

Suite of functions for the clustering of MHC II sequences into putative allele sequences. So far the pipeline has been applied to Antarctic fur seal (*Arctocephalus gazella*) MHC II DQB & DRB sequences as well as MHC II sequences of the Zebra finch (*Taeniopygia guttata*)

------------------------------------------------------------------------

This repository was forked from `https://github.com/nextgenusfs/mhc_cluster.git` [1]

<img src="figures/Fig_6_clustering.png" width="317" style="display: block; margin: auto;" />

The script `cluster_mhc2.py` allows to cluster high-quality Illumina reads into putative alleles using the [Unoise3](http://drive5.com/usearch/manual/cmd_unoise3.html) approach [2](<https://www.biorxiv.org/content/early/2016/10/15/081257>). Input sequences are `FASTQ` files generated on a Illumina MiSeq run. The following steps must be executed outside the clustering:

### Using cluster\_mhc2.py

Input files are expected to be in `fastq` format containing the barcodes for individual samples within the header as shown below. Barcodes are required for demultiplexing reads prior to mapping reads to the generated list of alleles.

    @MISEQ:279:000000000-AVVMJ:1:1101:14590:18831:N:0:barcodelabel=CAGAGAGGAAGGAGTA

For further information and a description of parameters see the helpfile:

``` bash
## first version 
cluster_mhc2.py -h
## latest version introducing vsearch 
cluster_mhc3.py -h 
```

The subfolder `lib` contains required functions and may contain hidden Markov models created from multiple sequence alignment of previously characterised MHC genes using [HMMER3](http://hmmer.org/).
Example using [MUSCLE](http://www.drive5.com/muscle/manual/) and [HMMER3](hmmer.org):

``` bash
##  Align sequences
muscle -in input_sequences.fasta -out aligned_sequences.afa
#   Create hidden markov model  
hmmbuild hmm aligned_sequences.afa
#   Create auxiliary files
hmmpress hmm
```

### Dependencies

-   [USEARCH10](http://www.drive5.com/usearch)
-   [VSEARCH](https://github.com/torognes/vsearch)

``` bash
wget https://github.com/torognes/vsearch/archive/v2.9.1.tar.gz
tar xzf v2.9.1.tar.gz
cd vsearch-2.9.1
./autogen.sh
./configure
make
make install  # as root or sudo make install
```

-   [HMMER3 v3.1b2](hmmer.org)

``` bash
tar zxf hmmer-3.1b2.tar.gz
cd hmmer-3.1b2
./configure
make
make check
```

-   [Biophyton](http://biopython.org/wiki/Download)

``` bash
pip install biopython
```

\*Installing Biophyton can easily cause some problems. See troubleshooting options here: <https://askubuntu.com/questions/677566/biopython-installation*>

-   [natsort](https://pypi.python.org/pypi/natsort)

``` bash
pip install natsort
```

### References

[1] Palmer JM, Berkman LK, Marquardt PE, Donner DM, Jusino MA, Lindner DL. Preliminary characterization of little brown bats (Myotis lucifugus) immune MHC II DRB alleles using next-generation sequencing. PeerJ PrePrints. 2016 Jan 21;4:e1662v1.

[2] Edgar, R.C., 2016. UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing. bioRxiv, p.081257.
