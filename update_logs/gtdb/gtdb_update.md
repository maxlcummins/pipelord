# Logs of commands used to add the new kraken DB

From the directory `/home/malcummi/Data/pipelord2_new/pipelord2_0`

First I made a new config_template file called config/config_template_gtdb.yaml where I added the path to the new database:

```/projects/AusGEM/databases/gtdb_r202```

Then I added a flag to the species_id.smk snakefile that requires users to change a field in the config to confirm they have enough RAM

```
if config['krakendb'] == "/projects/AusGEM/databases/gtdb_r202":
    if config['confirm_memory_available'] != True:
        print("Please confirm you have enough memory to run gtdb. You will need over 250 GB (!)")
        print("Otherwise pick a different kraken database (e.g. /projects/AusGEM/databases/kraken/bacteria)")
    else sys.exit()
```

## Run the original pipeline to generate a bracken summary

```
nohup snakemake --use-conda --conda-frontend mamba -T 10 --rerun-incomplete -s workflow/Snakefile --configfile config/config_template.yaml -j 50 -k -p results/test/QC_workflow/summaries/bracken_report.txt > update_logs/gtdb/kraken_normal_bacteria_db.out.log 2> update_logs/gtdb/kraken_normal_bacteria_db.err.log &
# Job ID = 57872
# Node = Jupiter14
```

Note in the config I accidentally set the prefix to test_gtdb when it should have just been test. Before running the gtdb job I need to change the name of the output folder for the above job to test.

## Run the original pipeline to generate a bracken summary

After renaming a folder as per above I ran the following

```
nohup snakemake --use-conda --conda-frontend mamba -T 10 --rerun-incomplete -s workflow/Snakefile --configfile config/config_template_gtdb.yaml -j 50 -k -p results/test_gtdb/QC_workflow/summaries/bracken_report.txt > update_logs/gtdb/kraken_gtdb_bacteria_db.out.log 2> update_logs/gtdb/kraken_gtdb_bacteria_db.err.log &
# Job ID = 37610
# Node = Jupiter14
```

As anticipated the gtdb was not processed for use with bracken. This will probably take some time...

I added a rule to ensure the database is processed for bracken before continuing.

## Run Bracken test with new db

```
nohup snakemake --use-conda --conda-frontend mamba -T 10 --rerun-incomplete -s workflow/Snakefile --configfile config/config_template_gtdb.yaml -j 24 -k -p results/test_gtdb/QC_workflow/summaries/bracken_report.txt > update_logs/gtdb/kraken_gtdb_bacteria_db.out2.log 2> update_logs/gtdb/kraken_gtdb_bacteria_db.err2.log &
# Job ID = 49151
# Node = Jupiter11
```

The above job kept failing throwing this error message:

`ERROR: Database library /projects/AusGEM/databases/gtdb_r202/library does not exist`

Advice from Torsten was to try to use kraken2-inspect on the database and then download the sequences using [ncbi-acc-download](https://github.com/kblin/ncbi-acc-download).

So I ran the following:

```
nohup kraken2-inspect --db /projects/AusGEM/databases/gtdb_r202 --threads 30 > update_logs/gtdb/kraken2_inspect.out.log 2> update_logs/gtdb/kraken2_inspect.err.log &
# Job ID = 23573
# Node = Jupiter11
```

This gave me a list of accession numbers from the database but the library folder contains a lot more than just the sequences of the genomes in the database. I think using bracken with GTDB is going to be very complicated so I have decided to put things on hold for the time being and instead try and set up the tool GTDBTK for species ID. This will do the job of kraken and bracken assuming that the checkM contamination screen is sufficient. Not sure if thats the case though and I should check with Torsten.
