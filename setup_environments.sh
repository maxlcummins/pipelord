#!/usr/bin/env bash

if test config/config_template.yaml; then
    echo "Detected pipelord config file at \`config/config_template.yaml\`"
    echo "Continuing with set-up"
elif test pipelord/config/config_template.yaml && mv pipelord; then
    echo "Detected pipelord config file at \`pipelord/config/config_template.yaml\`"
    echo "Moved into the pipelord directory with \`cd pipelord\`"
    echo "Continuing with set-up"
else
    SCRIPT=$(realpath "$0")
    SCRIPTPATH=$(dirname "$SCRIPT")
    echo "Cannot locate pipelord directory. Please change directory to $SCRIPTPATH"
fi

if which conda; then
    echo "Located conda"
    echo "Attepting to load environment \`snakemake\`..."
    conda activate snakemake
elif test ~/bin/my_python_path; then
    echo "Attempting to load wrapper script \`~/bin/my_python_path..."
    my_conda_path
    conda activate snakemake && echo "Successfully loaded snakemake environment"
else
    echo "Unable to load wrapper script or load snakemake environment."
    echo "Check the pipelord installation instructions at https://github.com/maxlcummins/pipelord"
    exit
fi 

if which mamba; then

    echo "Setting up environments for QC Module... Hang tight!"
    snakemake --use-conda --conda-frontend mamba -s workflow/QC_workflow.smk

    echo "Setting up environments for Genotyping Module... Hang tight!"
    snakemake --use-conda --conda-frontend mamba -s workflow/Snakefile

    echo "Setting up environments for Pangenomic Module... Hang tight!"
    snakemake --use-conda --conda-frontend mamba -s workflow/Treebuild_panaroo.smk

    echo "Successfully installed dependencies for all three modules!"
    echo "You can now perform a test run with the command:"
    echo "snakemake -j12 --use-conda"
elif which conda; then
    echo "Mamba not detected - will install environments using conda instead"
    echo "Setting up environments for QC Module... Hang tight!"
    snakemake --use-conda -s workflow/QC_workflow.smk

    echo "Setting up environments for Genotyping Module... Hang tight!"
    snakemake --use-conda -s workflow/Snakefile

    echo "Setting up environments for Pangenomic Module... Hang tight!"
    snakemake --use-conda -s workflow/Treebuild_panaroo.smk
    
    echo "Successfully installed dependencies for all three modules!"
    echo "You can now perform a test run with the command:"
    echo "snakemake -j12 --use-conda"
fi