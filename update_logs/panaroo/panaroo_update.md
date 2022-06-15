# Panaroo module development

## 1. Create a paneroo environment

```
#Create our environment
conda create -n panaroo -c conda-forge -c bioconda -c defaults panaroo

#Activate our environment
source activate panaroo

#Export our environment
conda env export > workflow/envs/panaroo.yaml
```

## 2. Create our snakefile

For this I modified the Roary workflow at workflows/Treebuild.smk. I replaced all instances of roary with paneroo, added a rule for a QC step from paneroo as well as a rule to download a database which it requires to run. I also modified the rule all to specify what files to make and changed the ruleorder to have it run the QC step first.

This new workflow is at workflows/Treebuild_paneroo.smk. I might rename it... We will see.


## 3. Config modification

Next I modified `config/config_template.yaml` to include newlines that specify paneroo parameters.

## 4. Test runs



Run our first test
```
nohup snakemake --use-conda --conda-frontend mamba -T 10 --rerun-incomplete -s workflow/Treebuild_panaroo.smk --configfile config/config_template.yaml -j 8 -p > test_panaroo.out.log 2> test_panaroo.err.log & 
# JobID = 447
# Jupiter11
```

Run our second test
```
nohup snakemake --use-conda --conda-frontend mamba --rerun-incomplete -s workflow/Treebuild_panaroo.smk --configfile config/config_template.yaml -j 8 -p > test_panaroo.out2.log 2> test_panaroo.err2.log & 
# JobID = 32089
# Jupiter11
```

I realised it kept breaking because one of the AVC genomes is intentionally crap quality. Running it again this time with the QC module in the rule_all

```
nohup snakemake --use-conda --conda-frontend mamba --rerun-incomplete -s workflow/Treebuild_panaroo.smk --configfile config/config_template.yaml -j 8 -p > test_panaroo.out3.log 2> test_panaroo.err3.log & 
# JobID = 63171
# Jupiter11
```

Now seems to be working but not sure rule order is correctly being followed. Rerunning to enable the full alignment process.

```
nohup snakemake --use-conda --conda-frontend mamba --rerun-incomplete -s workflow/Treebuild_panaroo.smk --configfile config/config_template.yaml -j 8 -p > test_panaroo.out4.log 2> test_panaroo.err4.log &
# JobID = 13958
# Jupiter11
```