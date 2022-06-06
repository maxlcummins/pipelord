# Log of commands and thought process for abritamr integration

Abritamr offers some advantages to our current AMR workflow so the plan is to substitute our AMR workflow for Abritamr.

## Update process

1. Install Abritamr and dependencies

As per the (github)[https://github.com/MDU-PHL/abritamr], with some small changes, install abritamr.

Instead of the recommended `conda create -n abritamr -c bioconda ncbi-amrfinder`, install using `conda create -n abritamr -c bioconda ncbi-amrfinderplus`.
This is necessary as the old conda link is dead.

Next I activated the environment and installed the databases for amrfinder, though I will need to write this step into the snakefile.


```
#Activate abritamr
conda activate abritamr
#Install amfinder databases
amrfinder -U
#Install abritamr
pip3 install abritamr
```

I actually had a bit of difficulty installing abritamr, I repeatedly got an error indicating that the software didnt even exist. I tried various approaches including using pip rather than pip3 but if i remember correctly what ended up working was using `conda install -c bioconda abritamr` and then updating from version 0.0.27 (?) to 1.0

2. Update the config file

Added the following lines to the config_template.yaml

```
genotype_modules:
#...
#An option to enable or disable abritamr
    run_abritamr: Yes
```

and

```
### abritamr config
#Options:
#Acinetobacter_baumannii,Campylobacter,Enterococcus_faecalis,Enterococcus_faecium
#Escherichia,Klebsiella,Salmonella,Staphylococcus_aureus,Staphylococcus_pseudintermedius
#Streptococcus_agalactiae,Streptococcus_pneumoniae,Streptococcus_pyogenes,Vibrio_cholerae
abritamr_species: "Escherichia"
```

3. Run Abritamr Snakemake workflow
```
nohup snakemake -p -j16 --use-conda --conda-frontend mamba results/test/abritamr/summary_matches.txt > update_logs/abritamr/abritamr_test_run.out.log 2> update_logs/abritamr/abritamr_test_run.err.log & 
# Job ID = 214959
# Mon May 30 16:18:44 AEST 2022
# Node = jupiter4
```

Now Abritamr will run but it will just vomit everything into the current WD. This is a problem I dont exactly know how to fix. The prefix flag only works for individual directories it seems.

Made some changes so that the summary files go where they are supposed to. Now I think all of the files for the individual outputs still get vomited into CWD. Lets try again and see:

```
nohup snakemake -p -j16 --use-conda --conda-frontend mamba results/test/abritamr/summary_matches.txt > update_logs/abritamr/abritamr_test_run.out2.log 2> update_logs/abritamr/abritamr_test_run.err2.log & 

# Job ID = 13165
# Mon Jun  6 12:09:52 AEST 2022
# Node = jupiter10
```

Still wasnt working for the suspected reason. I modified the snakefile to try and send the outputs to the results file by adding a parameter to the fofn to prefix the sample names. This time I just tried to generate the fofn though.

```
nohup snakemake -p -j16 --use-conda --conda-frontend mamba results/test/abritamr/fofns/full_fofn.txt > update_logs/abritamr/abritamr_test_run.out3.log 2> update_logs/abritamr/abritamr_test_run.err3.log & 

# Job ID = 23690
# Mon Jun  6 12:09:52 AEST 2022
# Node = jupiter10
```

With the fofn generated as intended, now I try and run the full abritamr workflow:

```
nohup snakemake -p -j16 --use-conda --conda-frontend mamba results/test/abritamr/summary_matches.txt > update_logs/abritamr/abritamr_test_run.out4.log 2> update_logs/abritamr/abritamr_test_run.err4.log &
# Job ID = 26141
# Mon Jun  6 12:09:52 AEST 2022
# Node = jupiter10
```