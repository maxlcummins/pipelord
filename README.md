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
<<<<<<< HEAD
nohup snakemake -p -j --use-conda 1> snakemake_out.log 2> snakemake_err.lo
=======
nohup snakemake -p -j --use-conda 1> snakemake_out.log 2> snakemake_err.log &
```
>>>>>>> 13f967b2655d51f55d21da940ba14699d4e69c0f
