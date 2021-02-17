===================
PointFinder
===================

This project documents PointFinder service


Documentation
=============

## What is it?

The PointFinder service contains one python script *PointFinder.py* which is the script of the lates
version of the PointFinder service. The method detects chromosomal mutations predictive of drug resistance based on WGS data.

## Content of the repository
1. PointFinder.py     - the program
2. test.fsa     - test fasta file


## Installation

Setting up PointFinder
```bash
# Go to wanted location for resfinder
cd /path/to/some/dir
# Clone and enter the pointfinder directory
git clone https://bitbucket.org/genomicepidemiology/pointfinder.git
cd pointfinder
```

Installing up the PointFinder database
```bash
git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git
```

Installing dependencies:

Biopython: http://biopython.org/DIST/docs/install/Installation.html
Blastn: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
or 
KMA
```bash
git clone https://bitbucket.org/genomicepidemiology/kma.git 
```

If using the KMA mapping, the databases .fsa should be index using kma index, please read the KMA readme for more information


## Usage 

The program can be invoked with the -h option to get help and more information of the service.

```bash
Usage: python3 PointFinder.py [options]

Options:

    -h HELP
                    Prints a message with options and information to the screen
    -p DATABASE
                    The path to where you have located the database folder
    -m 
                    The path to the location of blast-2.2.26 if it is not added
                    to the users path (see the install guide in 'README.md')
    -i INFILE
                    Your input file which needs to be preassembled partial
                    or complete genomes in fasta format
    -o OUTFOLDER
                    The folder you want to have your output files places.
                    If not specified the program will create a folder named
                    'Output' in which the result files will be stored.
    -s SPECIES
                    The PointFinder scheme you want to use. The options can be found in
                    the 'config' file in the database folder
    -u UNKNOWN
                    Incude if you want to have a report of all differences to the
                    reference gene including all that has not at this point been
                    associated with resistance (use '-u 1' to include this feature). 
```


## Web-server

A webserver implementing the methods is available as a feature of the ResFinder service at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/ResFinder-3.0/


## The Latest Version

The latest version can be found at
https://bitbucket.org/genomicepidemiology/pointfinder/overview

## Documentation

The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/pointfinder/overview.


Citation
=======

When using the method please cite:

Will be added soon

License
=======

Copyright (c) 2014, Ole Lund, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
