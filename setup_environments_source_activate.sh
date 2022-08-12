#!/usr/bin/env bash

if test config/config_template.yaml; then
    printf "\n"
    printf "Detected pipelord config file at \`config/config_template.yaml\`"
    printf "\n"
    printf "Continuing with set-up"
elif test pipelord/config/config_template.yaml && mv pipelord; then
    printf "\n"
    printf "Detected pipelord config file at \`pipelord/config/config_template.yaml\`"
    printf "\n"
    printf "Moved into the pipelord directory with \`cd pipelord\`"
    printf "\n"
    printf "Continuing with set-up"
else
    SCRIPT=$(realpath "$0")
    SCRIPTPATH=$(dirname "$SCRIPT")
    printf "Cannot locate pipelord directory. Please change directory to $SCRIPTPATH"
fi

if which conda > /dev/null; then
    printf "\n"
    printf "Located conda"
    printf "\n"
    printf "Attepting to load environment \`snakemake\`..."
    source activate snakemake
elif test ~/bin/my_python_path; then
    printf "Attempting to load wrapper script \`~/bin/my_python_path..."
    my_conda_path
    source activate snakemake && printf "Successfully loaded snakemake environment"
else
    printf "\n"
    printf "Unable to load wrapper script or load snakemake environment."
    printf "\n"
    printf "Check the pipelord installation instructions at https://github.com/maxlcummins/pipelord"
    exit
fi 

if which mamba; then

    printf "Setting up environments for QC Module... Hang tight!"
    printf "\n"
    snakemake --use-conda --conda-frontend mamba --conda-create-envs-only -j6 -s workflow/QC_workflow.smk
    
    printf "\n"
    printf "Setting up environments for Genotyping Module... Hang tight!"
    snakemake --use-conda --conda-frontend mamba --conda-create-envs-only -j6 -s workflow/Snakefile
    
    
    printf "\n"
    printf "Setting up environments for Pangenomic Module... Hang tight!"
    printf "\n"
    snakemake --use-conda --conda-frontend mamba --conda-create-envs-only -j6 -s workflow/Treebuild_panaroo.smk
    printf "\n"

    printf "\n"
    printf "Successfully installed dependencies for all three modules!"
    printf "\n"
    printf "You can now perform a test run with the command:"
    printf "\n"
    printf "snakemake -j12 --use-conda"
    printf "\n"

    printf "Successfully installed dependencies for all three modules using mamba!"
    printf "\n"
    printf "You can now perform a test run with the command:"
    printf "\n"
    printf "snakemake -j12 --use-conda"
    printf "\n"
    printf "Make sure to activate your snakemake environment first though..."
elif which conda; then
    printf "\n"
    printf "Mamba not detected - will install environments using conda instead"
    printf "\n"
    printf "Setting up environments for QC Module... Hang tight!"
    snakemake --use-conda --conda-create-envs-only -j6 -s workflow/QC_workflow.smk

    printf "\n"
    printf "Setting up environments for Genotyping Module... Hang tight!"
    printf "\n"
    snakemake --use-conda --conda-create-envs-only -j6 -s workflow/Snakefile

    printf "\n"
    printf "Setting up environments for Pangenomic Module... Hang tight!"
    printf "\n"
    snakemake --use-conda --conda-create-envs-only -j6 -s workflow/Treebuild_panaroo.smk
    printf "\n"

    printf "Successfully installed dependencies for all three modules!"
    printf "\n"
    printf "You can now perform a test run with the command:"
    printf "\n"
    printf "snakemake -j12 --use-conda"
    printf "\n"
    printf "Make sure to activate your snakemake environment first though..."
fi