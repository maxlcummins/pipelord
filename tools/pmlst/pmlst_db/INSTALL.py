#!/usr/bin/python3
import shutil, os, sys

# This scripts installs the pMLST database for using KMA
# KMA should be installed before running this script
# The scripts assumes that it is placed together with the pMLST scheme fasta files
# First clone the repository: git clone https://bitbucket.org/genomicepidemiology/pmlst_db.git

# Check if executable kma_index program is installed, if not promt the user for path

interactive = True
if len(sys.argv) >= 2:
    kma_index = sys.argv[1]
    if "non_interactive" in sys.argv:
        interactive = False
else:
    kma_index = "kma_index"

while shutil.which(kma_index) is None:
    if not interactive:
        sys.exit("KMA index program, {}, does not exist or is not executable".format(kma_index))
    ans = input("Please input path to executable kma_index program or enter 'q'/'quit' to exit:")
    if ans == "q" or ans == "quit":
        print("Exiting!\n\n \
               Please install executable KMA programs in order to install this database.\n\n \
               KMA can be obtained from bitbucked:\n\n\
               git clone https://bitbucket.org/genomicepidemiology/kma.git\n\n\
               KMA programs must afterwards be compiled:\n\n\
               gcc -O3 -o kma KMA.c -lm -lpthread\n\
               gcc -O3 -o kma_index KMA_index.c -lm")
        sys.exit()

    kma_index = ans

    if shutil.which(kma_index) is None:
        print("Path, {}, is not an executable path. Please provide absolute path\n".format(ans))


# Index databases


# Use config_file to go through database files
dirname = os.path.dirname(sys.argv[0])

config_file = open(os.path.join(dirname, "config"), "r")
for line in config_file:
    if line.startswith("#"):
        continue
    else:
        line = line.rstrip().split("\t")
        scheme = line[0] 
        f = os.path.join(dirname, scheme) + ".fsa"
        if not os.path.isfile(f):
            sys.exit("Database file '{}' does not exists".format(f))
        # for each dir index the fasta files
        os.system("{} -i {} -o {}".format(kma_index, f, os.path.join(dirname, scheme)))

config_file.close() 

print("Done")

