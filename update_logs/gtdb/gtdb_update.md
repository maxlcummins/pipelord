# Logs of commands used to add the new kraken DB

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

nohup snakemake --use-conda --conda-frontend mamba -T 10 --rerun-incomplete -s workflow/Snakefile --configfile config/config_template.yaml -j 50 -k -p results/test/QC_workflow/summaries/bracken_report.txt > kraken_normal_bacteria_db.out.log 2> kraken_normal_bacteria_db.err.log &
[1] 57872

Note in the config I accidentally set the prefix to test_gtdb when it should have just been test. Before running the gtdb job I need to change the name of the output folder for the above job to test.