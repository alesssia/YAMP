# Yet Another Metagenomic Pipeline (YAMP)

Thanks to the increased cost-effectiveness of high-throughput technologies, the number of studies focusing on microorganisms (bacteria, archaea, microbial eukaryotes, fungi, and viruses) and their connections with human health and diseases has surged, and, consequently, a plethora of approaches and software has been made available for their study, making it difficult to select the best methods and tools. 

Here we present Yet Another Metagenomic Pipeline (YAMP) that, starting from the raw sequencing data and having a strong focus on quality control, allows, within hours, the data processing up to the functional annotation (please refer to the YAMP [wiki](https://github.com/alesssia/YAMP/wiki) for more information).

YAMP is constructed on [Nextflow](https://www.nextflow.io/docs/latest/index.html), a framework based on the dataflow programming model, which allows writing workflows that are highly parallel, easily portable (including on distributed systems), and very flexible and customisable, characteristics which have been inherited by YAMP. New modules can be added easily and the existing ones can be customised -- even though we have already provided default parameters deriving from our own experience.

YAMP is accompanied by a set of (customisable) containers, that saves the users from the hassle of installing the required software, increasing, at the same time, the reproducibility of the YAMP results (see [Using Docker or Singularity](#using-docker-or-singularity)). 


## Table of contents

- [Citation](#Citation)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Other requirements](#other-requirements)
- [Usage](#usage)
- [Using Docker or Singularity](#using-docker-or-singularity)
- [Troubleshooting](#Troubleshooting)
- [Acknowledgements](#acknowledgements)
- [License](#license)
- [Changelog](#changelog)


## Citation

Please cite YAMP as:

> Visconti A., Martin T.C., and Falchi M., *"YAMP: a containerised workflow enabling reproducibility in metagenomics research"*, GigaScience (2018), [https://doi.org/10.1093/gigascience/giy072](https://doi.org/10.1093/gigascience/giy072)


## Dependencies

To run YAMP you will need to install Nextflow  (version 20.10 or higher), as explained [here](https://www.nextflow.io/docs/latest/getstarted.html). Please note that Nextflow requires BASH and [Java 7+](http://www.oracle.com/technetwork/java/javase/downloads/index.html). Both should be already available in most of the POSIX compatible systems (Linux, Solaris, OS X, etc). 

If you are using the containerised version of YAMP (as we strongly suggest), you should also install [Docker](https://www.docker.com) or [Singularity](http://singularity.lbl.gov/), as explained [here](https://docs.docker.com/engine/installation/) and [here](http://singularity.lbl.gov/docs-installation), respectively.

Once you have either Docker or Singularity up and running, you will not need to install anything additional tools. All the pieces of software are already specified in the YAMP pipeline and will be downloaded during the first run. Please refer to [Using Docker or Singularity](#using-docker-or-singularity) for more details. 

**For expert users only.** 
If you do not want to use the containerised version of YAMP, you must install the following pieces of software:

- fastQC v0.11.2+ ([http://www.bioinformatics.babraham.ac.uk/projects/fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc))
- BBmap v38.87+ ([https://sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap))
- Biobakery v3+ ([https://huttenhower.sph.harvard.edu/humann/](https://huttenhower.sph.harvard.edu/humann/))
- QIIME v2.0+ ([http://qiime.org](http://qiime.org))
- MultiQC v1.9-1+ ([https://multiqc.info/](https://multiqc.info/))

All of them should be in the system path with execute and read permission.

Following the links, you will find detailed instructions on how to install them, as explained by their developers. 
Notably, many of these tools are also available in [bioconda](https://anaconda.org/bioconda/).


## Installation

Clone the YAMP repository in a directory of your choice:

```
git clone https://github.com/alesssia/YAMP.git
```

A detailed description of what is included in the repository is available at [this wiki page](https://github.com/alesssia/YAMP/wiki/Folders-layout), while a description of the pre-set `config` files (used in our tutorials) is available [this wiki page](https://github.com/alesssia/YAMP/wiki/How-to-use-Nextflow-profiles),


## External databases

YAMP requires a set of databases that are queried during its execution. Some of them are already available with YAMP, others should be automatically downloaded either the first time you use the tool (MetaPhlAn), or using specialised scripts (HUMAnN), or should be created by the user. Specifically, you will need:

- A FASTA file listing the adapter sequences to remove in the trimming step. A basic version is provided in this repository (`./assests/data/adapters.fa`), but please note that this file may need to be customised.
- Two FASTA file describing synthetic contaminants. Basic versions are provided in this repository (`./assests/data/sequencing_artifacts.fa.gz` and `./assests/data/phix174_ill.ref.fa.gz`), but please note that both may need to be customised.
- A FASTA file describing the contaminant (pan)genome. This file should be created by the users according to the contaminants present in their dataset. When analysing human metagenomes, we recommend always including the human genome. 
- the BowTie2 database files for MetaPhlAn. These files should be downloaded the first time you run MetaPhlAn. 
- the ChocoPhlAn and UniRef databases for HUMAnN. Both can be downloaded directly by HUMAnN. Please refer to their [webpage](https://huttenhower.sph.harvard.edu/humann) for details. 

More details on these files and how to use and get them are available on our [wiki](https://github.com/alesssia/YAMP/wiki/Getting-started). Please read it carefully before using YAMP.


## Usage

The simplest way to use YAMP (after having satisfied all the dependencies and requirements) is by using the following command:

```
nextflow run YAMP.nf --reads1 myfile_R1.fq.gz --reads2 myfile_R2.fq.gz --prefix my_sample 
   --outdir output_folder --mode complete -profile base,docker
```

or you can run a test with the following command:

```
nextflow run YAMP.nf -profile test,docker
```

More information on the YAMP parameters and running mode are available in the [YAMP wiki](https://github.com/alesssia/YAMP/wiki), where there are also several tutorials.


## Using Docker or Singularity

YAMP takes advantage of a multi-image scenario. This means that each process will specify which container should be used, along with its version (as explained [here](https://github.com/alesssia/YAMP/wiki/Multi-image-scenario)).
YAMP also provides a `docker` and a `singularity` profile that can be used to tell Nextflow to enable the use of Docker/Singularity (as explained [here](https://github.com/alesssia/YAMP/wiki/How-to-use-Nextflow-profiles)), for instance using the following commands:

```
nextflow run YAMP.nf --reads1 myfile_R1.fq.gz --reads2 myfile_R2.fq.gz --prefix my_sample 
   --outdir output_folder --mode complete -profile base,docker
```


```
nextflow run YAMP.nf --reads1 myfile_R1.fq.gz --reads2 myfile_R2.fq.gz --prefix my_sample 
   --outdir output_folder --mode complete -profile base,singularity
```

Please note that Nextflow is not included in the Docker container and should be installed as explained [here](https://www.nextflow.io/docs/latest/getstarted.html).


## Troubleshooting

We have listed all known issues and solutions on this [wiki page](https://github.com/alesssia/YAMP/wiki/Troubleshooting). Please report any issue using the [GitHub platform](https://github.com/alesssia/YAMP/issues).


## Acknowledgements

Alessia would like to thank:
- Paolo Di Tommaso, for helping her in using Nextflow properly and his infinite patience;
- Brian Bushnell, for his helpful suggestions about how to successfully use the BBmap suite in a metagenomics context and for providing several useful resources;
- The [nf-core](https://nf-co.re/) team (especially [Phil Ewels](https://github.com/ewels) and [Harshil Patel](https://github.com/drpatelh)) for the resources they have provided (some of the code from version 0.9.5 is yours!) and for the nice discussions;
- All the users for their valuable feedbacks (especially Richard Davies [@richardjdavies](https://github.com/richardjdavies) and [Flavia Flaviani](https://github.com/flacchy)).


## License

YAMP is licensed under GNU GPL v3.


## Changelog

### 0.9.5.2 / 2021-04-16

Fixes:
* Modified default settings to avoid an error with Docker

### 0.9.5.1 / 2021-04-13

Fixes:
* Fixed a bug in `nextflow.config` 
* Fixed a typo in this README  (_Usage_ and _Using Docker or Singularity_ sections)

### 0.9.5 / 2021-03-23

Enhancements:
- Code refactoring
- Multi-image scenario is now used
- New process added to remove synthetic contaminants (separated from trimming)
- New process added to index contaminant (pan)genome
- New process added to deal with paired-end reads provided as input in `characterisation` mode
- Temporary files are now discarded 
- Profile files are now provided
- Testing and demo data are now provided

### 0.9.4.3 / 2019-12-06

Fixes:
* Solved problem with compressed files in `characterisation` mode
* Fixed warnings (Nextflow v19.10.0)

### 0.9.4.2 / 2018-09-14

Fixes:
* Solved problem in loading data in `complete` mode

Notes:
* YAMP now requires Nextflow version 0.29.x or higher

### 0.9.4.1 / 2018-04-24

Enhancements:
* QC'd files are now compressed (fq.gz) before being saved when `keepQCtmpfile` is true

Fixes:
* Solved problem in loading data when using single library layout
* Solved problem in loading data in 'characterisation` mode

### 0.9.4 / 2017-12-07

Enhancements:
* Improved logs
* Version and help message printed upon request

### 0.9.3.1 / 2017-10-04

Enhancements:
* Users no longer need to specify the number of threads and the maximum amount of memory -- both values are now read from the `nextflow.config` file

### 0.9.3 / 2017-08-30
 
 Enhancements:
 * YAMP can now handle both paired-end and single-end reads
 * The de-duplication step is now optional and can be skipped (default: true)
 
### 0.9.2 / 2017-07-10 

Enhancements:
* YAMP can now be run in three *"modes"* : < QC, characterisation, complete >.





