# System Requirements

PlasmidExtractor has been developed and tested on Linux-based systems (more specifically Ubuntu and Mint). Any Linux/Unix-based system (including MacOS) should be able to run the pipeline, while users of Windows-based systems will have to use Docker in order to have the pipeline run.

In terms of system specs the more CPUs your system has, the better, as almost every step of the pipeline is multi-threaded. RAM requirements can be somewhat heavy, with some parts of the pipeline using up to 20 GB due to the size of the database that PlasmidExtractor uses (a low-memory version of the database is on the to-do list). 32 GB or RAM is recommended, though machines with 24 GB may work.

### Installation Using Docker

As PlasmidExtractor has a fair number of dependencies, the easiest way to get it installed and working on your machine is using docker. Instructions on docker installation can be found [here](https://docs.docker.com/engine/installation/).

Once you have docker installed, you can load the PlasmidExtractor image by booting a terminal and typing:

`docker pull olcbioinformatics/cowpig`.

The pipeline can then be run with: 

`docker run -i -v /path/to/your/sequences/:/sequences olcbioinformatics/cowpig python3 /home/Extractor.py -i /sequences -o /sequences/output -sdb /home/new_database.fasta`

Where `/path/to/your/sequences/` is the folder with the FASTQ files you want analyzed. A folder called output will be created in your input folder that will contain your results.

Currently, the docker image is unable to create the output visualizations - this will be addressed in the future.

### Installation From Source

To install from source, you will first need to clone the GitHub repository. Open a terminal, navigate to where you want to download PlasmidExtractor, and type:

`git clone https://github.com/lowandrew/Plasmid_Assembler.git`

With that done, you will need to make sure that you have all of the depedencies for PlasmidExtractor installed and present on your $PATH. The dependencies for PlasmidExtractor are:

- [samtools >= 1.6](http://www.htslib.org/download/)
- [bcftools >= 1.6](http://www.htslib.org/download/)
- [ncbi-blast >= 2.2.31](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [bedtools >= 2.25.0](http://bedtools.readthedocs.io/en/latest/content/installation.html)
- [BBtools >= 37.23](https://jgi.doe.gov/data-and-tools/bbtools/)
- [mash >= 2.0](https://github.com/marbl/Mash/releases)
- [kmc >= 3.0](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download)
- [python >= 3.5](https://www.python.org/downloads/)

For instructions on how to add programs to your $PATH, click [here](https://askubuntu.com/questions/60218/how-to-add-a-directory-to-the-path).

You will also need to install the python packages that PlasmidExtractor needs to run. The packages need can be found in `requirements.txt`. To install all of them in one go, type:

`pip3 install -r requirements.txt`

If all of your dependencies are properly installed, you should now be able to run PlasmidExtractor. If installation of a dependency has not worked, you should get a _ModuleNotFoundError_ for the depedency that has not been able to run properly.

For visualization of sourmash results, you may need to install the python3-tk library. On a debian-based system, you can do this with: `apt-get install python3-tk`
