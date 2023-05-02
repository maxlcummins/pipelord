#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas as pd
import os.path
import re

# Directing python to the input from snakemake

mobsuite = snakemake.input["mobsuite"]

# Create an empty list

dataframes = []

# Initialise a for loop for processing our pMLST outputs

for file in mobsuite:

    # Read it in as df

    df = pd.read_csv(file, sep="\t", header=0, on_bad_lines="skip")

    # Append our dataframes
    dataframes.append(df)

# Concatenate our dfs into a single one
combined = pd.concat(dataframes)

# Write our output file to text

combined.to_csv(snakemake.output[0], sep="\t", index=False)
