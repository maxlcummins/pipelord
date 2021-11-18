# pipelord2.0

## Disclaimer
While I have made every effort to ensure these scripts produce reliable and accurate output, it is important that you spot check your outputs to ensure they are correct.

## Installation and setup

```
# Clone the git repository and enter it
git clone https://github.com/maxlcummins/pipelord2_0.git
cd pipelord2_0

# Setup and installation - this may take some time
snakemake -j --use-conda
```

## Configuration

Configuration of the snakemake workflow can be done via the configuration file in `misc`

## Usage

```
nohup snakemake --resources mem_mb=450000 -j 20 -p --use-conda --configfile config/config_template.yaml --use-conda 1> snakemake_out.log 2> snakemake_err.log &

PROCESS_ID=$!
echo "JOB_ID =" "$PROCESS_ID" > mostrecentjobid.txt
```
