# pipelord

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14048334.svg)](https://doi.org/10.5281/zenodo.14048334)

## About

This repo contains a genomic analysis pipeline used for analysing mostly *E. coli* genomes, though it also sees use in analysing *Salmonella* genomes and some other species as well, however note that using this tool for analysing other species may cause certain tools to either break or will otherwise just limit the usefullness of this tool in the first place.

## Disclaimer
While I have made every effort to ensure these scripts produce reliable and accurate output, it is important that you spot check your outputs to ensure they are correct :)

The documentation is still a work in progress but if you need any help don't hestitate to post an issue :)

## Modules/Workflows

This pipeline has three main modules:

1. Quality control
2. Genotyping
3. Phylogenetics

All individual rules/analyses can be turned on or off in the configuration file.

### 1. Quality control

Reads are first filtered with fastp before being assembled with shovill. Assembly stats are analysed using assembly-stats. Species ID is determined using Kraken and Bracken and a contamination and completeness screen is carried out using CheckM. GUNC is also run, though usually I disable this as the run-time can be a bit long and I haven't found the output especially helpful yet for my purposes, though this might change.

This information is then consolidated into a tab delimited QC report.

### 2. Genotyping

The word genotyping here is used a bit loosely but this step performs a variety of analyses using:

* Abricate to screen for virulence, AMR, IS elements and other custom gene databases. New databases can be added and added to the configuration file to be included for use.
* Antibiotic resistance gene detection is performed using abritamr
* Point mutations associated with AMR are detected using CGE's pointfinder (there is some redundancy here with abritamr which now has this functionality)
* MLST is performed using mlst
* Plasmid MLST is performed using pMLST
* Fimbrial adhesins are typed using fimtyper
* Phylogroup detection is performed with Clermontyping and EZClermont
* Serotyping is performed with ECtyper
* Salmonella pathogenicity islands are detected using SPIfinder
* Genome annotation is performed with Prokka
* cgMLST is performed with cgMLSTfinder - this will be changed to chewBBACA
* Plasmid screening is performed using a slightly older version of abricate (screening of large reference genes seemed to break after 0.9.8).

Most of these rules output tab delimited summary files relevant for each rule/analysis, some utilise some downstream scripts in R but not many of these.

### 3. Phylogenetics

This workflow runs a pangenome derived phylogenetic analysis utilising:

* Prokka
* Panaroo
* SNP-sites
* IQtree

Pairwise SNP distances are also determined using SNP-dists.

## Installation and setup

Because there are many dependencies to install across all of the different analytical tools the first time you run the tool it can take some time to install everything via conda. Because of this I have added a script to install environments for the three modules `setup_environments.sh`

```
# Clone the git repository and enter it
git clone https://github.com/maxlcummins/pipelord.git
cd pipelord

# Setup and installation - this may take some time
sh setup_environments.sh
```

## Configuration

Configuration of the snakemake workflow can be done via the configuration file in `config`. It is a good idea to just make a new one though by copying the template and editing it accordingly, then specifying it in the snakemake command (e.g. `--configfile config/new_template.yaml`).

## Usage

Note: Change your config file accordingly if you created a new one.

### QC Workflow

```
#Run job in the background
nohup snakemake --resources mem_mb={max_memory_per_rule} -j 20 -p --use-conda --configfile config/config_template.yaml --use-conda -s workflows/QC_workflow.smk 1> QC.log 2> QC.log &

#Save the job ID so you can cancel it later if need be
PROCESS_ID=$!
echo "JOB_ID =" "$PROCESS_ID" > QC_jobID.txt
```

### Genotyping

```
#Run job in the background
nohup snakemake --resources mem_mb={max_memory_per_rule} -j 20 -p --use-conda --configfile config/config_template.yaml --use-conda -s workflows/Snakefile 1> genotype.log 2> genotype.log &

#Save the job ID so you can cancel it later if need be
PROCESS_ID=$!
echo "JOB_ID =" "$PROCESS_ID" > genotype_jobID.txt
```
### Pangenomic Phylogenetics

```
#Run job in the background
nohup snakemake --resources mem_mb={max_memory_per_rule} -j 20 -p --use-conda --configfile config/config_template.yaml --use-conda -s workflows/Treebuild_panaroo.smk 1> treebuild.log 2> treebuild.log &

#Save the job ID so you can cancel it later if need be
PROCESS_ID=$!
echo "JOB_ID =" "$PROCESS_ID" > treebuild_jobID.txt
```

## Database update

To update your databases you can use a command like the following - make sure the db dir path is correct and change the database appropriately

```
abricate-get_db --dbdir resources/dbs/abricate --db vfdb --force
```
