# pipelord2.0

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
nohup snakemake -p -j --use-conda 1> snakemake_out.log 2> snakemake_err.log &
```
