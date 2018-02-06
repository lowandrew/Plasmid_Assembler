# System Requirements

PlasmidExtractor has been developed and tested on Linux-based systems (more specifically Ubuntu and Mint). Any Linux/Unix-based system (including MacOS) should be able to run the pipeline, while users of Windows-based systems will have to use Docker in order to have the pipeline run.

In terms of system specs the more CPUs your system has, the better, as almost every step of the pipeline is multi-threaded. RAM requirements can be somewhat heavy, with some parts of the pipeline more than 20 GB due to the size of the database that PlasmidExtractor uses. 32 GB or RAM is recommended. If your machine does not meet these requirements, PlasmidExtractor can be run in a low-memory mode, where memory usage should peak at ~7GB, enabling the pipeline to run on a much wider range of machines.


## Databases

The databases necessary for PlasmidExtractor are hosted on FigShare, and need to be downloaded no matter what installation method you choose.

To download them, use the following command in the folder you want to download the databases to:

`wget https://ndownloader.figshare.com/files/9827323  && tar xf 9827323 && rm 9827323`

This will create a folder in you current working directory called databases which has all that you need to run PlasmidExtractor.

## Installation Using Conda

Installing using conda is recommended, as it takes care of installing external program depedencies as well as the PlasmidExtractor script itself.

If you don't already have conda installed, you can find instructions on the [conda ReadTheDocs page.](https://conda.io/docs/user-guide/install/index.html)

With conda installed, all you'll need to do is download the PlasmidExtractor environment, which has been uploaded to the Anaconda cloud. To do this, use this command:

`conda env create lowandrew/plasmidextractor`

You should then be able to enter the conda environment:

`source activate plasmidextractor`

Within the environment, typing `PlasmidExtractor.py` will run the PlasmidExtractor pipeline - see the [Usage](usage.md) section for instructions on how to run PlasmidExtractor.

## Installing Using Pip

#### Executable

PlasmidExtractor can also be installed using pip. Use of a virtual environment for PlasmidExtractor is highly recommended. To create a virtualenv:

- Create an empty directory (i.e. `mkdir ~/Virtual_Environments/PlasmidExtractor`)
- Virtualenv that directory (`virtualenv -p /usr/bin/python3 ~/Virtual_Environments/PlasmidExtractor`)
- Activate the virtualenv (`source ~/Virtual_Environments/PlasmidExtractor/bin/activate`)
- Install PlasmidExtractor - this should also install any python packages necessary for PlasmidExtractor to run. (`pip install plasmidextractor`)

With this done, you'll need to make sure that any necessary dependencies are installed.

#### Mandatory Dependencies 
With pip installation done, you will need to make sure that you have all of the depedencies for PlasmidExtractor installed and present on your $PATH. The dependencies for PlasmidExtractor are:

- [samtools >= 1.6](http://www.htslib.org/download/)
- [bcftools >= 1.6](http://www.htslib.org/download/)
- [ncbi-blast >= 2.2.31](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [bedtools >= 2.25.0](http://bedtools.readthedocs.io/en/latest/content/installation.html)
- [BBtools >= 37.23](https://jgi.doe.gov/data-and-tools/bbtools/)
- [mash >= 2.0](https://github.com/marbl/Mash/releases)
- [kmc >= 3.0](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download)
- [python >= 3.5](https://www.python.org/downloads/)
- [tabix >= 1.2.1](https://github.com/samtools/htslib)

For instructions on how to add programs to your $PATH, click [here](https://askubuntu.com/questions/60218/how-to-add-a-directory-to-the-path).

If all of your dependencies are properly installed, you should now be able to run PlasmidExtractor. If installation of a dependency has not worked, you will get an error message printed to screen telling you which depedencies are not installed. Note that the check only looks for an executable and does _not_ check the version, so you may need to check that you have the correct versions installed manually. 

#### Optional Dependencies

- [pilon >= 1.22](https://github.com/broadinstitute/pilon/releases)

If you have Pilon installed, PlasmidExtractor will try to use it in order to polish the 
consensus sequence created. In order to have this work, you need to have an executable called `pilon`
accessible somewhere on your $PATH. 

To do this, download pilon from the above link, use chmod to make the jar file executable (`chmod 755 pilon-1.22.jar`),
and then rename the executable to `pilon` and add it to your $PATH.

