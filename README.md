**Clustering Illumina MiSeq paired-end reads of MHC sequences**

Clustering of paired-end reads is conducted using the script `mhc_clustering.py`.

```
mhc_clustering.py -h
```

The subfolder `lib` contains required functions and may contain hidden marko models created from multiple sequence alignment of published mhc alleles using **HMMER3**.  

*This repository was forked from https://github.com/nextgenusfs/mhc_cluster.git` and subsequently modified. See details there.*

**Dependent software and installation**
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
