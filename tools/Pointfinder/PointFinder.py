#!/usr/bin/env python3
import sys, os, re, math, pprint
import argparse
from cgecore.blaster import Blaster

##########################################################################
# FUNCTIONS
##########################################################################

def get_gene_seqs(database_path, gene):
    """
    This function takes the database path and a gene name as inputs and 
    returns the gene sequence contained in the file given by the gene name
    """
    gene_path = database_path + "/" + gene + ".fsa"
    gene_seq = ""
    # Open fasta file
    with open(gene_path) as gene_file:
        header = gene_file.readline()
        for line in gene_file:
            seq = line.strip()
            gene_seq += seq
    return gene_seq

def get_db_mutations(mut_db_path, gene_list, res_stop_codons):
    """
    This function opens the file resistenss-overview.txt, and reads the
    content into a dict of dicts. The dict will contain information about
    all known mutations given in the database. This dict is returned.
    """

    # Open resistens-overview.txt
    try:
        drugfile = open(mut_db_path, "r")
    except:
        sys.exit("Wrong path: %s"%(mut_db_path))
    
    # Initiate variables
    known_mutations = dict()
    drug_genes = dict()
    known_stop_codon = dict()
    indelflag = False
    stopcodonflag = False
   
    # Go throug mutation file line by line
    for line in drugfile:
        # Ignore headers and check where the indel section starts
        if line.startswith("#"):
            if "indel" in line.lower():
                indelflag = True
            elif "stop codon" in line.lower():
                stopcodonflag = True
            else:
                stopcodonflag = False
            continue
        # Ignore empty lines 
        if line.strip() == "":
            continue

        # Strip data entries
        mutation = [data.strip() for data in line.strip().split("\t")]

        # Extract all info on the line (even though it is not all used)
        gene_ID = mutation[0]

        # Only consider mutations in genes found in the gene list
        if gene_ID in gene_list:
            gene_name = mutation[1]
            #no_of_mut = int(mutation[2])
            mut_pos = int(mutation[2])
            ref_codon = mutation[3]
            ref_aa = mutation[4]
            alt_aa = mutation[5].split(",")
            res_drug = mutation[6].replace("\t", " ")
            pmid = mutation[7].split(",")

            # Check if resistance is known to be caused by a stop codon in the gene
            if ("*" in alt_aa and res_stop_codons != 'specified') or (res_stop_codons == 'specified' and stopcodonflag == True):
                if gene_ID not in known_stop_codon:
                    known_stop_codon[gene_ID] = {"pos": [], "drug": res_drug}
                known_stop_codon[gene_ID]["pos"].append(mut_pos)

            # Add genes associated with drug resistance to drug_genes dict
            drug_lst = res_drug.split(",")
            for drug in drug_lst:
                drug = drug.upper()
                if drug not in drug_genes:
                    drug_genes[drug] = []
                if gene_ID not in drug_genes[drug]:
                    drug_genes[drug].append(gene_ID)

            # Initiate empty dict to store relevant mutation information
            mut_info = dict()
            
            # Save need mutation info with pmid cooresponding to the amino acid change
            for i in range(len(alt_aa)):
                try:
                    mut_info[alt_aa[i]] = {"gene_name": gene_name, "drug": res_drug, "pmid": pmid[i]}
                except IndexError:
                    mut_info[alt_aa[i]] = {"gene_name": gene_name, "drug": res_drug, "pmid": "-"}

            # Add all possible types of mutations to the dict
            if gene_ID not in known_mutations:
                known_mutations[gene_ID] = {"sub" : dict(), "ins" : dict(), "del" : dict()}

            # Check for the type of mutation
            if indelflag == False:
                mutation_type = "sub"
            else:
                mutation_type = ref_aa

    	    # Save mutations positions with required information given in mut_info
            if mut_pos not in known_mutations[gene_ID][mutation_type]:
                known_mutations[gene_ID][mutation_type][mut_pos] = dict() 
            for aa in alt_aa:
                known_mutations[gene_ID][mutation_type][mut_pos][aa] = mut_info[aa]

    drugfile.close()

    # Check that all genes in the gene list has known mutations
    for gene in gene_list:
        if gene not in known_mutations:
            known_mutations[gene] = {"sub" : dict(), "ins" : dict(), "del" : dict()}
    return known_mutations, drug_genes, known_stop_codon

def KMA(inputfile_1, gene_list, kma_db, out_path, sample_name, min_cov, mapping_path):
    """
    This function is called when KMA is the method of choice. The 
    function calls kma externally and waits for it to finish. 
    The kma output files with the prefixes .res and .aln are parsed 
    throught to obtain the required alignment informations. The subject 
    and query sequences as well as the start and stop position, 
    coverage, and subject length are stored in a results directory
    which is returned in the end.
    """
    
    # Get full path to input of output files
    inputfile_1 = os.path.abspath(inputfile_1)

    kma_outfile = os.path.abspath(out_path + "/kma_out_" + sample_name)

    kma_cmd = "%s -i %s -t_db %s -o %s -1t1 -gapopen -5 -gapextend -2 -penalty -3 -reward 1"%(mapping_path, inputfile_1, kma_db, kma_outfile) # -ID 90
    
    # Call KMA
    os.system(kma_cmd)
    if os.path.isfile(kma_outfile + ".aln") == False:
        os.system(kma_cmd)

    # Fetch kma output files
    align_filename = kma_outfile + ".aln"
    res_filename = kma_outfile + ".res"
    
    results = dict()

    # Open KMA result file
    with open(res_filename, "r") as res_file:
        header = res_file.readline()

        # Parse through each line
        for line in res_file:
            data = [data.strip() for data in line.split("\t")]
            gene = data[0]

            # Check if gene one of the user specified genes  
            if gene not in gene_list:
                continue

            # Store subject length and coverage
            sbjct_len = int(data[3])
            identity = float(data[6])
            coverage = float(data[7])

            # Result dictionary assumes that more hits can occur
            if gene not in results:
                hit = '1'
                results[gene] = dict()

            # Gene will only be there once with KMA
            else:
                hit = str(len(results[gene])) +1

            results[gene][hit] = dict()
            results[gene][hit]['sbjct_length'] = sbjct_len
            results[gene][hit]['coverage'] = coverage / 100
            results[gene][hit]["sbjct_string"] = []
            results[gene][hit]["query_string"] = []
            results[gene][hit]["homology"] = []
            results[gene][hit]['identity'] = identity

    # Open KMA alignment file   
    with open(align_filename, "r") as align_file:
        hit_no = dict()
        gene = ""

        # Parse through alignments
        for line in align_file:

            # Check when a new gene alignment start
            if line.startswith("#"):
                gene = line[1:].strip()
                if gene not in hit_no:
                    hit_no[gene] = str(1)
                else:
                    hit_no[gene] += str(int(hit_no[gene]) + 1)

            else:
                # Check if gene is one of the user specified genes             
                if gene in results:
                    if hit_no[gene] not in results[gene]:
                        sys.exit("Unexpected database redundency")
                    line_data = line.split("\t")[-1].strip()
                    if line.startswith("template"):
                        results[gene][hit_no[gene]]["sbjct_string"] += [line_data]
                    elif line.startswith("query"):
                        results[gene][hit_no[gene]]["query_string"] += [line_data]
                    else:
                        results[gene][hit_no[gene]]["homology"] += [line_data]
    
    # Concatinate all sequences lists and find subject start and subject end
    seq_start_search_str = re.compile("^-*(\w+)")
    seq_end_search_str = re.compile("\w+(-*)$")
    for gene in gene_list:
        if gene in results:
            for hit in results[gene]:
                results[gene][hit]['sbjct_string'] = "".join(results[gene][hit]['sbjct_string'])
                results[gene][hit]['query_string'] = "".join(results[gene][hit]['query_string'])
                results[gene][hit]['homology'] = "".join(results[gene][hit]['homology'])
 
   	        
                seq_start_object = seq_start_search_str.search(results[gene][hit]['query_string'])
                sbjct_start = seq_start_object.start(1) + 1

                seq_end_object = seq_end_search_str.search(results[gene][hit]['query_string'])
                sbjct_end = seq_end_object.start(1) + 1

                results[gene][hit]['query_string'] = results[gene][hit]['query_string'][sbjct_start-1:sbjct_end-1]
                results[gene][hit]['sbjct_string'] = results[gene][hit]['sbjct_string'][sbjct_start-1:sbjct_end-1]


                #if sbjct_start:
                results[gene][hit]["sbjct_start"] = sbjct_start
                results[gene][hit]["sbjct_end"] = sbjct_end
        else:
           results[gene] = ""

    return results

def find_best_sequence(hits_found, specie_path, gene, silent_N_flag):
    """
    This function takes the list hits_found as argument. This contains all 
    hits found for the blast search of one gene. A hit includes the subjct 
    sequence, the query, and the start and stop position of the allignment 
    corresponding to the subject sequence. This function finds the best 
    hit by concatinating sequences of found hits. If different overlap 
    sequences occurr these are saved in the list alternative_overlaps. The 
    subject and query sequence of the concatinated sequence to gether with 
    alternative overlaps and the corresponding start stop
    positions are returned.
    """

    # Get information from the fisrt hit found	
    all_start = hits_found[0][0]
    current_end = hits_found[0][1]
    final_sbjct = hits_found[0][2]
    final_qry = hits_found[0][3] 
    sbjct_len = hits_found[0][4]  

    alternative_overlaps = []
    
    # Check if more then one hit was found within the same gene
    for i in range(len(hits_found)-1):

        # Save information from previous hit
        pre_block_start = hits_found[i][0]
        pre_block_end = hits_found[i][1]
        pre_sbjct = hits_found[i][2]
        pre_qry = hits_found[i][3]

	# Save information from next hit
        next_block_start = hits_found[i+1][0]
        next_block_end = hits_found[i+1][1]
        next_sbjct = hits_found[i+1][2]
        next_qry = hits_found[i+1][3]

        # Check for overlapping sequences, collaps them and save alternative overlaps if any
        if next_block_start <= current_end:

            # Find overlap start and take gaps into account     
            pos_count = 0
            overlap_pos = pre_block_start
            for i in range(len(pre_sbjct)):

                # Stop loop if overlap_start position is reached
                if overlap_pos == next_block_start:
                    overlap_start = pos_count
                    break
                if pre_sbjct[i] != "-":
                    overlap_pos += 1
                pos_count += 1
            
            # Find overlap length and add next sequence to final sequence 
            if len(pre_sbjct[overlap_start:]) > len(next_sbjct):
                #  <--------->
                #     <--->
                overlap_len = len(next_sbjct)
                overlap_end_pos = next_block_end
            else:
                #  <--------->
                #        <--------->
                overlap_len = len(pre_sbjct[overlap_start:])
                overlap_end_pos = pre_block_end

                # Update current end
                current_end = next_block_end

                # Use the entire pre sequence and add the last part of the next sequence
                final_sbjct += next_sbjct[overlap_len:]
                final_qry += next_qry[overlap_len:]
                
            # Find query overlap sequences
            pre_qry_overlap = pre_qry[overlap_start : (overlap_start + overlap_len)] # can work for both types of overlap
            next_qry_overlap = next_qry[:overlap_len]
            sbjct_overlap = next_sbjct[:overlap_len]

            # If alternative query overlap excist save it
            if pre_qry_overlap != next_qry_overlap:
                print("OVERLAP WARNING:")
                print(pre_qry_overlap, "\n", next_qry_overlap)

                # Save alternative overlaps
                alternative_overlaps += [(next_block_start, overlap_end_pos, sbjct_overlap, next_qry_overlap)]
        
        elif next_block_start > current_end:
            #  <------->
            #              <-------> 
            gap_size = next_block_start - current_end - 1
            final_qry += "N"*gap_size
            if silent_N_flag:
                final_sbjct += "N"*gap_size
            else:
                ref_seq = get_gene_seqs(specie_path, gene)
                final_sbjct += ref_seq[pre_block_end:pre_block_end+gap_size]

            current_end = next_block_end
            final_sbjct += next_sbjct
            final_qry += next_qry
    
    # Calculate coverage
    no_call = final_qry.upper().count("N")
    coverage = (current_end - all_start +1 - no_call) / float(sbjct_len)
    
    # Calculate identity
    equal = 0
    not_equal = 0
    for i in range(len(final_qry)):
        if final_qry[i].upper() != "N":
            if final_qry[i].upper() == final_sbjct[i].upper():
                equal += 1
            else:
                not_equal += 1
    identity = equal/float(equal + not_equal)

    return final_sbjct, final_qry, all_start, current_end, alternative_overlaps, coverage, identity 


def find_mismatches(gene, sbjct_start, sbjct_seq, qry_seq, alternative_overlaps = []):
    """
    This function finds mis matches between two sequeces. Depending on the
    the sequence type either the function find_codon_mismatches or 
    find_nucleotid_mismatches are called, if the sequences contains both 
    a promoter and a coding region both functions are called. The function 
    can also call it self if alternative overlaps is give. All found mis 
    matches are returned
    """

    # Initiate the mis_matches list that will store all found mis matcehs
    mis_matches = []

    # Find mis matches in RNA genes
    if gene in RNA_gene_list:
        mis_matches += find_nucleotid_mismatches(sbjct_start, sbjct_seq, qry_seq)
    else:
        # Check if the gene sequence is with a promoter
        regex = r"promoter_size_(\d+)(?:bp)"
        promtr_gene_objt = re.search(regex, gene)

        # Check for promoter sequences
        if promtr_gene_objt:

            # Get promoter length
            promtr_len = int(promtr_gene_objt.group(1))

            # Extract promoter sequence, while considering gaps	
            # --------agt-->----
            #    ---->?
            if sbjct_start <= promtr_len:

                #Find position in sbjct sequence where promoter ends
                promtr_end = 0
                nuc_count = sbjct_start - 1
                for i in range(len(sbjct_seq)): 
                    promtr_end += 1
                    if sbjct_seq[i] != "-":
                        nuc_count += 1
                    if nuc_count == promtr_len:
                        break    

                # Check if only a part of the promoter is found
                #--------agt-->----
                # ----
                promtr_sbjct_start = -1
                if nuc_count < promtr_len:
                    promtr_sbjct_start = nuc_count - promtr_len

                # Get promoter part of subject and query 
                sbjct_promtr_seq = sbjct_seq[:promtr_end]
                qry_promtr_seq = qry_seq[:promtr_end]

                
                # For promoter part find nucleotide mis matches
                mis_matches += find_nucleotid_mismatches(promtr_sbjct_start, sbjct_promtr_seq, qry_promtr_seq, promoter = True)
                
                # Check if gene is also found
                #--------agt-->----
                #     -----------           
                if (sbjct_start + len(sbjct_seq.replace("-", ""))) > promtr_len:
                    sbjct_gene_seq = sbjct_seq[promtr_end:]
                    qry_gene_seq = qry_seq[promtr_end:]
                    sbjct_gene_start = 1

                    # Find mismatches in gene part
                    mis_matches += find_codon_mismatches(sbjct_gene_start, sbjct_gene_seq, qry_gene_seq)
            
            # No promoter, only gene is found
            #--------agt-->----
            #            ----- 
            else:
                sbjct_gene_start = sbjct_start - promtr_len
            
                # Find mismatches in gene part
                mis_matches += find_codon_mismatches(sbjct_gene_start, sbjct_seq, qry_seq)
            
        else:
            # Find mismatches in gene
            mis_matches += find_codon_mismatches(sbjct_start, sbjct_seq, qry_seq)

    # Find mismatches in alternative overlaps if any
    for overlap in alternative_overlaps:
        mis_matches += find_mismatches(gene, overlap[0], overlap[2], overlap[3])

    return mis_matches


def find_nucleotid_mismatches(sbjct_start, sbjct_seq, qry_seq, promoter = False):
    """ 
    This function takes two alligned sequence (subject and query), and the 
    position on the subject where the alignment starts. The sequences are 
    compared one nucleotide at a time. If mis matches are found they are 
    saved. If a gap is found the function find_nuc_indel is called to find
    the entire indel and it is also saved into the list mis_matches. If 
    promoter sequences are given as arguments, these are reversed the and 
    the absolut value of the sequence position  used, but when mutations
    are saved the negative value and det reverse sequences are saved in 
    mis_mathces.
    """

    # Initiate the mis_matches list that will store all found mis matcehs
    mis_matches = []
    
    sbjct_start = abs(sbjct_start)
    seq_pos = sbjct_start

    # Set variables depending on promoter status
    factor = 1
    mut_prefix = "r."
    if promoter == True:
        factor = (-1)
        mut_prefix = "n."
        # Reverse promoter sequences
        sbjct_seq = sbjct_seq[::-1]
        qry_seq = qry_seq[::-1]    
    
    # Go through sequences one nucleotide at a time
    shift = 0
    for index in range(sbjct_start - 1, len(sbjct_seq)):
        mut_name = mut_prefix
        mut = ""

        # Shift index according to gaps
        i = index + shift
        
        # If the end of the sequence is reached, stop
        if i == len(sbjct_seq):
            break
        
        sbjct_nuc = sbjct_seq[i]
        qry_nuc = qry_seq[i]
        
        # Check for mis matches
        if sbjct_nuc.upper() != qry_nuc.upper():
            
            # check for insertions and deletions
            if sbjct_nuc == "-" or qry_nuc == "-":
                if sbjct_nuc == "-":
                    mut = "ins"
                    indel_start_pos = (seq_pos -1) *factor
                    indel_end_pos = seq_pos * factor
                    indel = find_nuc_indel(sbjct_seq[i:], qry_seq[i:])
                else:
                    mut = "del"
                    indel_start_pos = seq_pos * factor
                    indel = find_nuc_indel(qry_seq[i:], sbjct_seq[i:]) 
                    indel_end_pos = (seq_pos + len(indel) - 1) * factor  
                    seq_pos += len(indel) - 1
                
                # Shift the index to the end of the indel
                shift += len(indel) - 1
                                     
                # Write mutation name, depending on sequnce
                if len(indel) == 1 and mut == "del":
                    mut_name += str(indel_start_pos) + mut + indel
                else:
                    if promoter == True:

                        # Reverse the sequence and the start and end positions
                        indel = indel[::-1]
                        temp = indel_start_pos
                        indel_start_pos = indel_end_pos
                        indel_end_pos = temp
    
                    mut_name += str(indel_start_pos) + "_" +str(indel_end_pos) + mut + indel  
                
                mis_matches += [[mut, seq_pos * factor, seq_pos * factor, indel, mut_name, mut, indel]]
            
            # Check for substitutions mutations
            else:
                mut = "sub"
                mut_name += str(seq_pos * factor) + sbjct_nuc + ">" + qry_nuc
                mis_matches += [[mut, seq_pos * factor, seq_pos * factor, qry_nuc, mut_name, sbjct_nuc, qry_nuc]]

        # Increment sequence position
        if mut != "ins":
            seq_pos += 1

    return mis_matches


def find_nuc_indel(gapped_seq, indel_seq):
    """
    This function finds the entire indel missing in from a gapped sequence 
    compared to the indel_seqeunce. It is assumes that the sequences start
    with the first position of the gap.
    """
    ref_indel = indel_seq[0]
    for j in range(1,len(gapped_seq)):
        if gapped_seq[j] == "-":
            ref_indel += indel_seq[j]
        else:
            break
    return ref_indel


def aa(codon):
    """
    This function converts a codon to an amino acid. If the codon is not 
    valid an error message is given, or else, the amino acid is returned.
    """
    codon = codon.upper()
    aa = {"ATT": "I", "ATC": "I", "ATA": "I",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "TTA": "L", "TTG": "L",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "TTT": "F", "TTC": "F",
            "ATG": "M",
            "TGT": "C", "TGC": "C",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
            "TAT": "Y", "TAC": "Y",
            "TGG": "W",
            "CAA": "Q", "CAG": "Q",
            "AAT": "N", "AAC": "N",
            "CAT": "H", "CAC": "H",
            "GAA": "E", "GAG": "E",
            "GAT": "D", "GAC": "D",
            "AAA": "K", "AAG": "K",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
            "TAA": "*", "TAG": "*", "TGA": "*"} 

    # Translate valid codon 
    try:
        amino_a = aa[codon]
    except KeyError:
        amino_a = "?"
    return amino_a


def get_codon(seq, codon_no, start_offset):
    """
    This function takes a sequece and a codon number and returns the codon
    found in the sequence at that position 
    """
    seq = seq.replace("-","")
    codon_start_pos = int(codon_no - 1)*3 - start_offset
    codon = seq[codon_start_pos:codon_start_pos + 3]
    return codon


def name_insertion(sbjct_seq, codon_no, sbjct_nucs, aa_alt, start_offset):
    """
    This function is used to name a insertion mutation based on the HGVS 
    recommendation. 
    """
    start_codon_no = codon_no - 1
    if len(sbjct_nucs) == 3:
        start_codon_no = codon_no
    start_codon = get_codon(sbjct_seq, start_codon_no, start_offset)
    end_codon = get_codon(sbjct_seq, codon_no, start_offset)
    pos_name = "p.%s%d_%s%dins%s"%(aa(start_codon), start_codon_no, aa(end_codon), codon_no, aa_alt)
    return pos_name


def name_deletion(sbjct_seq, sbjct_rf_indel, sbjct_nucs, codon_no, aa_alt, start_offset, mutation = "del"):
    """
    """
    del_codon = get_codon(sbjct_seq, codon_no, start_offset)
    pos_name = "p.%s%d"%(aa(del_codon), codon_no)
    if len(sbjct_rf_indel) == 3:
        return pos_name + mutation
    end_codon_no = codon_no + math.ceil(len(sbjct_nucs)/3) - 1
    end_codon = get_codon(sbjct_seq, end_codon_no, start_offset)
    pos_name += "_%s%d%s"%(aa(end_codon), end_codon_no, mutation)
    if mutation == "delins":
        pos_name += aa_alt
    return pos_name


def name_indel_mutation(sbjct_seq, indel, sbjct_rf_indel, qry_rf_indel, codon_no, mut, start_offset):
    """
    This function serves to name the individual mutations dependently on 
    the type of the mutation.
    """
    # Get the subject and query sequences without gaps
    sbjct_nucs = sbjct_rf_indel.replace("-", "")
    qry_nucs = qry_rf_indel.replace("-", "")

    # Translate nucleotides to amino acids
    aa_ref = ""
    aa_alt = ""
    for i in range(0, len(sbjct_nucs), 3):
        aa_ref += aa(sbjct_nucs[i:i+3])
    for i in range(0, len(qry_nucs), 3):
        aa_alt += aa(qry_nucs[i:i+3])

    # Identify the gapped sequence 
    if mut == "ins":
        gapped_seq = sbjct_rf_indel
    else:
        gapped_seq = qry_rf_indel
    gap_size = gapped_seq.count("-")

    # Write mutation names
    if gap_size < 3 and len(sbjct_nucs) ==3 and len(qry_nucs) == 3:

        # Write mutation name for substitution mutation
        mut_name = "p.%s%d%s"%(aa(sbjct_nucs), codon_no, aa(qry_nucs))

    elif len(gapped_seq) == gap_size:
        if mut == "ins":

            # Write mutation name for insertion mutation
            mut_name = name_insertion(sbjct_seq, codon_no, sbjct_nucs, aa_alt, start_offset)
            aa_ref = mut
        else:

            # Write mutation name for deletion mutation
            mut_name = name_deletion(sbjct_seq, sbjct_rf_indel, sbjct_nucs, codon_no, aa_alt, start_offset, mutation = "del")
            aa_alt = mut

    # Check for delins - mix of insertion and deletion
    else:

        # Write mutation name for a mixed insertion and deletion mutation
        mut_name = name_deletion(sbjct_seq, sbjct_rf_indel, sbjct_nucs, codon_no, aa_alt, start_offset, mutation = "delins")

    # Check for frameshift
    if gapped_seq.count("-")%3 != 0:
        # Add the frameshift tag to mutation name
        mut_name += " - Frameshift"

    return mut_name, aa_ref, aa_alt


def get_inframe_gap(seq, nucs_needed = 3):
    """
    This funtion takes a sequnece starting with a gap or the complementary 
    seqeuence to the gap, and the number of nucleotides that the seqeunce
    should contain in order to maintain the correct reading frame. The 
    sequence is gone through and the number of non-gap characters are 
    counted. When the number has reach the number of needed nucleotides 
    the indel is returned. If the indel is a 'clean' insert or deletion 
    that starts in the start of a codon and can be divided by 3, then only 
    the gap is returned.
    """
    nuc_count = 0
    gap_indel  = ""
    nucs = ""
    for i in range(len(seq)):

        # Check if the character is not a gap
        if seq[i] != "-":

            # Check if the indel is a 'clean' 
            # i.e. if the insert or deletion starts at the first nucleotide in the codon and can be divided by 3
            if gap_indel.count("-") == len(gap_indel) and gap_indel.count("-") >= 3 and len(gap_indel) != 0:
                return gap_indel
            nuc_count += 1
        gap_indel += seq[i]

        # If the number of nucleotides in the indel equals the amount needed for the indel, the indel is returned.
        if nuc_count == nucs_needed:
            return gap_indel

    # This will only happen if the gap is in the very end of a sequence
    return gap_indel


def get_indels(sbjct_seq, qry_seq, start_pos):
    """
    This function uses regex to find inserts and deletions in sequences
    given as arguments. A list of these indels are returned. The list 
    includes, type of mutations(ins/del), subject codon no of found 
    mutation, subject sequence position, insert/deletions nucleotide 
    sequence, and the affected qry codon no.
    """

    seqs = [sbjct_seq, qry_seq]
    indels = []
    gap_obj = re.compile(r"-+")
    for i in range(len(seqs)):
        for match in gap_obj.finditer(seqs[i]):
            pos = int(match.start())
            gap = match.group()

            # Find position of the mutation corresponding to the subject sequence
            sbj_pos = len(sbjct_seq[:pos].replace("-","")) + start_pos
    
            # Get indel sequence and the affected sequences in sbjct and qry in the reading frame
            indel = seqs[abs(i-1)][pos:pos+len(gap)]                   

            # Find codon number for mutation
            codon_no = int(math.ceil((sbj_pos)/3))
            qry_pos = len(qry_seq[:pos].replace("-","")) + start_pos
            qry_codon = int(math.ceil((qry_pos)/3))
            if i == 0:
                mut = "ins"
            else:
                mut = "del"
            
            indels.append( [mut, codon_no, sbj_pos, indel, qry_codon])

    # Sort indels based on codon position and sequence position
    indels = sorted(indels, key = lambda x:(x[1],x[2]))
    
    return indels


def find_codon_mismatches(sbjct_start, sbjct_seq, qry_seq):
    """
    This function takes two alligned sequence (subject and query), and
    the position on the subject where the alignment starts. The sequences 
    are compared codon by codon. If a mis matches is found it is saved in 
    'mis_matches'. If a gap is found the function get_inframe_gap is used 
    to find the indel sequence and keep the sequence in the correct 
    reading frame. The function translate_indel is used to name indel 
    mutations and translate the indels to amino acids
    The function returns a list of tuples containing all needed informations
    about the mutation in order to look it up in the database dict known 
    mutation and the with the output files the the user.  
    """
    mis_matches = []
    
    # Find start pos of first codon in frame, i_start
    codon_offset = (sbjct_start-1) % 3
    i_start = 0
    if codon_offset != 0:
        i_start = 3 - codon_offset
    sbjct_start = sbjct_start + i_start
    
    # Set sequences in frame
    sbjct_seq = sbjct_seq[i_start:]
    qry_seq = qry_seq[i_start:]
    
    # Find codon number of the first codon in the sequence, start at 0
    codon_no = int((sbjct_start-1) / 3) # 1,2,3 start on 0
  
    # s_shift and q_shift are used when gaps appears
    q_shift = 0
    s_shift = 0
    mut_no = 0
    
    # Find inserts and deletions in sequence
    indel_no = 0
    indels = get_indels(sbjct_seq, qry_seq, sbjct_start)

    # Go through sequence and save mutations when found
    for index in range(0, len(sbjct_seq), 3):
        # Count codon number
        codon_no += 1
        
        # Shift index according to gaps
        s_i = index + s_shift
        q_i = index + q_shift

        # Get codons
        sbjct_codon = sbjct_seq[s_i:s_i+3]
        qry_codon =  qry_seq[q_i:q_i+3]
        
        if len(sbjct_seq[s_i:].replace("-","")) + len(qry_codon[q_i:].replace("-","")) < 6:
            break

        # Check for mutations
        if sbjct_codon.upper() != qry_codon.upper():

            # Check for codon insertions and deletions and frameshift mutations
            if "-" in sbjct_codon or "-" in qry_codon:

                # Get indel info
                try:
                    indel_data = indels[indel_no]
                except IndexError:
                    print(sbjct_codon, qry_codon)
                    print(indels)
                    print(gene, indel_data, indel_no)
                mut = indel_data[0]
                codon_no_indel = indel_data[1]                
                seq_pos = indel_data[2] + sbjct_start - 1
                indel = indel_data[3]
                indel_no +=1
                
                # Get the affected sequence in frame for both for sbjct and qry 
                if mut == "ins":
                    sbjct_rf_indel = get_inframe_gap(sbjct_seq[s_i:], 3)
                    qry_rf_indel = get_inframe_gap(qry_seq[q_i:], int(math.floor(len(sbjct_rf_indel)/3) *3))                    
                else:
                    qry_rf_indel = get_inframe_gap(qry_seq[q_i:], 3)
                    sbjct_rf_indel = get_inframe_gap(sbjct_seq[s_i:], int(math.floor(len(qry_rf_indel)/3) *3))
                                        
                mut_name, aa_ref, aa_alt = name_indel_mutation(sbjct_seq, indel, sbjct_rf_indel, qry_rf_indel, codon_no, mut, sbjct_start - 1)

                # Set index to the correct reading frame after the indel gap 
                shift_diff_before = abs(s_shift - q_shift)
                s_shift += len(sbjct_rf_indel) - 3
                q_shift += len(qry_rf_indel) - 3
                shift_diff = abs(s_shift - q_shift)

                if shift_diff_before != 0 and shift_diff %3 == 0:

                    if s_shift > q_shift:
                        nucs_needed = int((len(sbjct_rf_indel)/3) *3) + shift_diff
                        pre_qry_indel = qry_rf_indel
                        qry_rf_indel = get_inframe_gap(qry_seq[q_i:], nucs_needed)
                        q_shift += len(qry_rf_indel) - len(pre_qry_indel)
                    elif q_shift > s_shift:
                        nucs_needed = int((len(qry_rf_indel)/3)*3) + shift_diff
                        pre_sbjct_indel = sbjct_rf_indel
                        sbjct_rf_indel = get_inframe_gap(sbjct_seq[s_i:], nucs_needed)
                        s_shift += len(sbjct_rf_indel) - len(pre_sbjct_indel)

                    
                    mut_name, aa_ref, aa_alt = name_indel_mutation(sbjct_seq, indel, sbjct_rf_indel, qry_rf_indel, codon_no, mut, sbjct_start - 1) 

                    if "Frameshift" in mut_name:
                        mut_name = mut_name.split("-")[0] + "- Frame restored"

                mis_matches += [[mut, codon_no_indel, seq_pos, indel, mut_name, sbjct_rf_indel, qry_rf_indel, aa_ref, aa_alt]]

                # Check if the next mutation in the indels list is in the current codon
                # Find the number of individul gaps in the evaluated sequence
                no_of_indels = len(re.findall("\-\w", sbjct_rf_indel)) + len(re.findall("\-\w", qry_rf_indel))
                if no_of_indels > 1:

                    for j in range(indel_no, indel_no + no_of_indels - 1):
                        try:
                            indel_data = indels[j]
                        except IndexError:
                            sys.exit("indel_data list is out of range, bug!")
                        mut = indel_data[0]
                        codon_no_indel = indel_data[1]                
                        seq_pos = indel_data[2] + sbjct_start - 1
                        indel = indel_data[3]
                        indel_no +=1
                        mis_matches += [[mut, codon_no_indel, seq_pos, indel, mut_name, sbjct_rf_indel, qry_rf_indel, aa_ref, aa_alt]]

                # Set codon number, and save nucleotides from out of frame mutations                
                if mut == "del":
                    codon_no += int((len(sbjct_rf_indel) - 3)/3)
                # If evaluated insert is only gaps codon_no should not increment
                elif sbjct_rf_indel.count("-") == len(sbjct_rf_indel):
                    codon_no -= 1
            
            # Check of point mutations
            else:
                mut = "sub"
                aa_ref = aa(sbjct_codon)
                aa_alt = aa(qry_codon)
                
                if aa_ref != aa_alt:
                    # End search for mutation if a premature stop codon is found
                    mut_name = "p." + aa_ref + str(codon_no) + aa_alt

                    mis_matches += [[mut, codon_no, codon_no, aa_alt, mut_name, sbjct_codon, qry_codon, aa_ref, aa_alt]]

            # If a Premature stop codon occur report it an stop the loop
            try:
                if mis_matches[-1][-1] == "*":
                    mut_name += " - Premature stop codon"
                    mis_matches[-1][4] = mis_matches[-1][4].split("-")[0] + " - Premature stop codon"
                    break
            except IndexError:
                pass

    # Sort mutations on position
    mis_matches = sorted(mis_matches, key = lambda x:x[1])
    
    return mis_matches

def write_output(gene, gene_name, mis_matches, known_mutations, known_stop_codon, unknown_flag, GENES):
    """
    This function takes a gene name a list of mis matches found betreewn subject and query of
    this gene, the dictionary of known mutation in the point finder database, and the flag telling 
    weather the user wants unknown mutations to be reported.
    All mis matches are looked up in the known mutation dict to se if the mutation is known, 
    and in this case what drug resistence it causes.
    The funtions returns a 3 strings that are used as output to the users.
    One string is only tab seperated and contains the mutations listed line by line. 
    If the unknown flag is set to true it will contain both known and unknown mutations. 
    The next string contains only known mutation and are given in in a format that is easy to
    convert to HTML. The last string is the HTML tab sting from the unknown mutations.
    """
    RNA = False
    known_header = "Mutation\tNucleotide change\tAmino acid change\tResistance\tPMID\n"
    unknown_header = "Mutation\tNucleotide change\tAmino acid change\n"
    if gene in RNA_gene_list:
        RNA = True
        known_header = "Mutation\tNucleotide change\tResistance\tPMID\n"
        unknown_header = "Mutation\tNucleotide change\n"

    known_lst = []
    unknown_lst = []
    all_results_lst = [] 
    output_mut = []
    stop_codons = []

    # Go through each mutation    
    for i in range(len(mis_matches)):
        m_type = mis_matches[i][0]
        pos = mis_matches[i][1] # sort on pos?
        look_up_pos = mis_matches[i][2]
        look_up_mut = mis_matches[i][3]
        mut_name = mis_matches[i][4]
        nuc_ref = mis_matches[i][5]
        nuc_alt = mis_matches[i][6]
        ref =  mis_matches[i][-2]
        alt = mis_matches[i][-1]

        # First index in list indicates if mutation is known
        output_mut += [[]]
        #output_mut[i] = [0]

    	# Define output vaiables
        codon_change = nuc_ref + " -> " + nuc_alt
        aa_change = ref + " -> " + alt
        if RNA == True:
            aa_change = "RNA mutations"
        elif pos < 0:
            aa_change = "Promoter mutations"
        
        # Check if mutation is known
        gene_mut_name, resistence, pmid = look_up_known_muts(known_mutations, known_stop_codon, gene, look_up_pos, look_up_mut, m_type, gene_name, mut_name)
        gene_mut_name = gene_mut_name + " " + mut_name

        output_mut[i] = [gene_mut_name, codon_change, aa_change, resistence, pmid]
        
        # Add mutation to output strings for known mutations 
        if resistence != "Unknown":
            if RNA == True:
                # don't include the amino acid change field for RNA mutations
                known_lst += ["\t".join(output_mut[i][:2]) + "\t" + "\t".join(output_mut[i][3:])]
            else:
                known_lst += ["\t".join(output_mut[i])]
            all_results_lst += ["\t".join(output_mut[i])]

        # Add mutation to output strings for unknown mutations 
        else:
            if RNA == True:
                unknown_lst += ["\t".join(output_mut[i][:2])]
            else:
                unknown_lst += ["\t".join(output_mut[i][:3])]
            if unknown_flag == True:
                all_results_lst += ["\t".join(output_mut[i])]

        # Check that you do not print two equal lines (can happen it two indels occure in the same codon)
        if len(output_mut) > 1:
            if output_mut[i] == output_mut[i-1]:
                if resistence != "Unknown":
                    known_lst = known_lst[:-1]
                    all_results_lst = all_results_lst[:-1]
                else:
                    unknown_lst = unknown_lst[:-1]
                    if unknown_flag == True:
                        all_results_lst = all_results_lst[:-1]
        if "Premature stop codon" in mut_name:
            sbjct_len = GENES[gene]['sbjct_len']
            qry_len = pos * 3 
            prec_truckat = round(((float(sbjct_len) - qry_len )/ float(sbjct_len)) * 100, 2) 
            perc = "%"
            stop_codons.append("Premature stop codon in %s, %.2f%s lost"%(gene, prec_truckat, perc))

    # Creat final strings
    all_results = "\n".join(all_results_lst)
    total_known_str = "" 
    total_unknown_str = ""

    # Check if there are only unknown mutations
    resistence_lst = [res for mut in output_mut for res in mut[3].split(",")]

    # Save known mutations
    unknown_no = resistence_lst.count("Unknown")
    if unknown_no < len(resistence_lst):
        total_known_str = known_header + "\n".join(known_lst)
    else:
        total_known_str = "No known mutations found in %s"%gene_name

    # Save unknown mutations
    if unknown_no > 0:
        total_unknown_str = unknown_header + "\n".join(unknown_lst)
    else:
        total_unknown_str = "No unknown mutations found in %s"%gene_name

    return all_results, total_known_str, total_unknown_str, resistence_lst + stop_codons

def check_path(args_input, input_type, error_type = "Input Error", no_return = False):
    """
    
    """
    args_lst = []
    if type(args_input) == str:
        args_lst.append(args_input)
    elif type(args_input) == list:
        args_lst = args_input
    for infile in args_lst:
        if not os.path.exists(infile):
            sys.exit("%s: %s, %s does not exist"%(error_type, input_type, infile))
        elif os.stat(infile).st_size == 0:
            print("Warning: %s, %s is empty"%(input_type, infile))
    if no_return: 
        return None
    return args_input

def get_file_content(file_path, fst_char_only = False):
    """
    
    """
    try:
        with open(file_path, "r") as infile:
            line_lst = []
            for line in infile:
                line = line.strip()
                if line != "":
                    line_lst.append(line)
                if fst_char_only:
                    return line_lst[0][0]
    except IOError as error:
        sys.exit(error)
    return line_lst

def look_up_known_muts(known_mutations, known_stop_codon, gene, pos, found_mut, mut, gene_name, mut_name):
    """
    
    """
    resistence = "Unknown"
    pmid = "-"
    found_mut = found_mut.upper()

    if mut == "del":
        for i, i_pos in enumerate(range(pos, pos+len(found_mut))):
            if i_pos in known_mutations[gene]["del"]:
                for known_indel in known_mutations[gene]["del"][i_pos]:
                    if found_mut[i:len(known_indel)+i] == known_indel and len(found_mut)%3 == len(known_indel)%3:
                        resistence = known_mutations[gene]["del"][i_pos][known_indel]['drug']
                        pmid = known_mutations[gene]["del"][i_pos][known_indel]['pmid']
                        gene_name =  known_mutations[gene]["del"][i_pos][known_indel]['gene_name']
                        break 
    else:
        if pos in known_mutations[gene][mut]:
            if found_mut in known_mutations[gene][mut][pos]:
                resistence = known_mutations[gene][mut][pos][found_mut]['drug']
                pmid = known_mutations[gene][mut][pos][found_mut]['pmid']
                gene_name =  known_mutations[gene][mut][pos][found_mut]['gene_name']

    if resistence == "Unknown" and ("*" in found_mut or "*" in mut_name) and gene in known_stop_codon:
        if res_stop_codons == 'early':
            if max(known_stop_codon[gene]["pos"]) > pos:
                resistence = known_stop_codon[gene]["drug"]
                print("Resistance added:", resistence)
        elif res_stop_codons == 'all' or res_stop_codons == 'specified':
            resistence = known_stop_codon[gene]["drug"]
            print("Resistance added:", resistence)

    return gene_name, resistence, pmid
            
def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = argparse.ArgumentParser(description="This program predicting resistance associated with chromosomal mutations based on WGS data", 
                                 prog="PointFinder.py")

# positional arguments
parser.add_argument("-i", "--inputfiles", nargs ='+', help="Input file, for blast 1 fasta file, for KMA fastq file(s)", required=True) 
parser.add_argument("-o", "--out_path", help="Path to existing output directory", required=True) 
parser.add_argument("-s", "--species", help="Species for point mutation detetion", required=True) 
parser.add_argument("-p", "--databasePath", dest="db_path",help="Path to the databases", default='/databases', required=True)
parser.add_argument("-m", "--method", choices=["kma", "blastn"], required=True)
parser.add_argument("-m_p", "--method_path", help="Path to kma or blastn program depending on the chosen methode", required=True) 

# optional arguments
parser.add_argument("-n", "--no_Ns", dest="no_N_output",help="Silence the output where Ns are found", action='store_false')
parser.add_argument("-t", "--threshold", dest="threshold",help="Blast threshold for identity",type=restricted_float, default=0.9)
parser.add_argument("-l", "--min_cov", dest="min_cov",help="Minimum coverage", type=restricted_float, default=0.6)
parser.add_argument("-u", "--unknown_mut", dest="unknown_mutations",help="Show all mutations found even if it's unknown to the resistance database", action='store_true')
parser.add_argument("-g", "--specific_genes", nargs ='+', dest="specific_genes",help="Specify genes existing in the database to search for - if none is specified all genes are used")
parser.add_argument("-r", "--stop_codons", dest="res_stop_codons",help="predict res on premature stop codons", choices=['early','all','specified'])

args = parser.parse_args()

# If no arguments are given print usage message and exit
if len(sys.argv) == 1:
    sys.exit("Usage: " + parser.usage)

# Check if valid user inputs is provided
inputfiles = check_path(args.inputfiles, "Input file")
out_path = check_path(args.out_path, "Output directory")
db_path = check_path(args.db_path, "Database directory")

# Check database installation checks
specie = args.species
specie_path = db_path  + "/" + specie
specie_path = check_path(db_path  + "/" + specie, "Incorrectly installed database", error_type = "Database Error")
db_gene_lst_path = check_path(specie_path + "/genes.txt", "Gene list file", error_type = "Database Error")
db_RNAgene_lst_path = check_path(specie_path + "/RNA_genes.txt", "RNA gene list file", error_type = "Database Error")

gene_list = get_file_content(db_gene_lst_path)
RNA_gene_list = get_file_content(db_RNAgene_lst_path)

# Creat user defined gene_list if applied  
if args.specific_genes:
    genes_specified = []
    for gene in args.specific_genes:

        # Check that the genes are valid
        if gene not in gene_list:
            sys.exit("Input Error: Specified gene not recognised (%s)\nChoose one or more of the following genes:\n%s"%(gene, "\n".join(gene_list)))
        genes_specified.append(gene)

    # Change the gene_list to the user defined gene_list
    gene_list = genes_specified

# Check path for mapping program, KMA or BLASTN
method = args.method
if method == "blastn":
    if len(inputfiles) != 1:
        sys.exit("Input Error: Blast was chosen as mapping method only 1 input file requied, not %s"%(len(inputfiles)))

    # Check if input file is in fasta format
    infile_fst_char = get_file_content(inputfiles[0], fst_char_only = True)
    if infile_fst_char != ">":
        sys.exit("Input Error: Input file is not in fasta format, first character found in file: %s" %(infile_fst_char))

    # Check that all gene sequences exist in the database
    ref_gene_paths = [specie_path + "/" + gene + ".fsa" for gene in gene_list]
    check_path(ref_gene_paths, "Invalid reference database for Blast", error_type = "Database Error", no_return = True)
elif method == "kma":

    # Check that kma indexed database exist
    kma_db = db_path + "/" + specie + "/" + specie
    kma_db_files = [kma_db + ".length.b"]
    #kma_db_files = [kma_db + ".b", kma_db + ".length.b", kma_db + ".name", kma_db + ".seq.b", kma_db + ".index.b", kma_db + ".comp.b"]
    check_path(kma_db_files, "Invalid indexed KMA database", error_type = "Database Error", no_return = True)
else:
    sys.exit("Input Error: No valid mapping method chosen, choose between -k /path/to/kma or -b /path/to/blastn")

mapping_path = check_path(args.method_path, "Invalid %s path"%method)

# Set coverage and identity thresholds 	
min_cov =  args.min_cov 
threshold = float(args.threshold)

# Set a flag for specifing if unknown mutations are reported
unknown_flag = args.unknown_mutations
silent_N_flag = args.no_N_output

if args.res_stop_codons:
    res_stop_codons = args.res_stop_codons
else:
    res_stop_codons = 'off'

###############################################################################
### MAIN
###############################################################################

# Open resistens-overview file and extract mutation information 
known_mutations, drug_genes, known_stop_codon = get_db_mutations(specie_path + "/resistens-overview.txt", gene_list, res_stop_codons)

# Get sample name
filename = inputfiles[0].split("/")[-1]
sample_name = filename.split(".")[0]#.split("_")[0]
if sample_name == "":
    sample_name = filename

# Call BLAST or KMA
if method == "blastn":

    # Call blast and parse output
    min_cov_blast = 0.01

    res = Blaster(inputfiles[0], gene_list, specie_path, out_path, min_cov_blast, threshold, mapping_path, cut_off=False, max_target_seqs=100)

    results = res.results
    qry_align = res.gene_align_query
    homo_align = res.gene_align_homo
    sbjct_align = res.gene_align_sbjct

else:
    # run KMA
    inputfiles = ' '.join(inputfiles)
    results = KMA(inputfiles, gene_list, kma_db, out_path, sample_name, min_cov, mapping_path)

GENES = dict()

# Find gene hit with the largest coverage
for gene, hits in results.items():
    if gene == "excluded":
        continue

    # Save all hits in the list 'hits_found'
    hits_found = []
    GENES[gene] = dict()

    # Check for hits in blast results
    if type(hits) is dict:
        GENES[gene]['found'] = 'partially'

        # Check coverage for each hit, if coverage is 100% save gene directly else save hit info. if only one hit, save gene directly
        for hit in hits.keys(): 
            hit_coverage = results[gene][hit]['coverage']

            # Append tuble with subject start and end positions to the list 'hits_found'
            hits_found += [(results[gene][hit]['sbjct_start'],results[gene][hit]['sbjct_end'], results[gene][hit]['sbjct_string'], results[gene][hit]['query_string'], results[gene][hit]['sbjct_length'])] 

            # If coverage is 100% change found to yes               
            if hit_coverage == 1.0:
                GENES[gene]['found'] = 'yes'               

        # Sort positions found
        hits_found = sorted(hits_found, key = lambda x:x[0])

        # Find best hit by concatenating sequences if more hits exist 
        final_sbjct, final_qry, all_start, all_end, alternative_overlaps, total_coverage, total_identity  = find_best_sequence(hits_found, specie_path, gene, silent_N_flag)
        
        # Save blast output of gene in GENES
        if total_coverage >= min_cov and total_identity >= threshold:
            GENES[gene]['coverage'] = total_coverage
            GENES[gene]['identity'] = total_identity
            GENES[gene]['sbjct_string'] = final_sbjct
            GENES[gene]['query_string'] = final_qry
            GENES[gene]['sbjct_start'] =  all_start
            GENES[gene]['sbjct_end'] = all_end
            GENES[gene]['sbjct_len'] = results[gene][hit]['sbjct_length']
            GENES[gene]['alternative_overlaps'] = alternative_overlaps
            GENES[gene]['mis_matches'] = []
        else:
            # Gene not found above given coverage
            GENES[gene]['coverage'] = total_coverage
            GENES[gene]['identity'] = total_identity
            
            if total_coverage < min_cov:
                GENES[gene]['found'] = 'Gene found with coverage (%f) below  minimum coverage threshold: %s' %(total_coverage, min_cov)
            else:
                GENES[gene]['found'] = 'Gene found with identity (%f) below  minimum identity threshold: %s' %(total_identity, threshold)
    else:
        # Gene not found!
        GENES[gene]['found'] = 'Gene not found'
        GENES[gene]['coverage'] = 0


# Find known mutations and write output files
# Output filenames
output_files = ["results.tsv", "HTMLtable.txt", "prediction.txt"]
   
# Initiate output stings with header
output_strings = ["Mutation\tNucleotide change\tAmino acid change\tResistance\tPMID", 
                  "Chromosomal point mutations - Results\nSpecies: %s\nMapping methode: %s\n\n\nKnown Mutations\n"%(specie, method), ""]
total_unknown_str = ""
unique_drug_list = []

# Find mutation in gene if gene is found
for gene in GENES:
    # Start writing output string (to HTML tab file)
    gene_name = gene
    regex = r"promoter_size_(\d+)(?:bp)"
    promtr_gene_objt = re.search(regex, gene)
    if promtr_gene_objt:
        gene_name = gene.split("_")[0]
    output_strings[1] += "\n%s\n"%(gene_name)        

    # Check if gene is found
    if GENES[gene]['found'] == 'yes' or GENES[gene]['found'] == 'partially':
        sbjct_start = GENES[gene]['sbjct_start']
        sbjct_seq = GENES[gene]['sbjct_string'] 
        qry_seq = GENES[gene]['query_string']
        alternative_overlaps = GENES[gene]['alternative_overlaps']

        # Find and save mis_matches in gene
        GENES[gene]['mis_matches'] = find_mismatches(gene, sbjct_start, sbjct_seq, qry_seq, alternative_overlaps)

	# If gene isn't found write the reason saved in GENES[gene]['found']
    if GENES[gene]['found'] != 'yes' and GENES[gene]['found'] != 'partially':
        output_strings[1] += GENES[gene]['found'] + "\n"
    else:
        # Check if any mutations was found           
        if len(GENES[gene]['mis_matches']) < 1:
            output_strings[1]  += "No mutations found in %s"%(gene_name)
            if GENES[gene]['coverage'] != 1:
                output_strings[1] += " (coverage: %.2f)"%(GENES[gene]['coverage'] * 100)
            output_strings[1]  += "\n" 
        else:
            # Write mutations found to output file             
            total_unknown_str += "\n%s\n"%(gene_name)
            all_results, total_known, total_unknown, drug_list = write_output(gene, gene_name, GENES[gene]['mis_matches'], known_mutations, known_stop_codon, unknown_flag, GENES)

            # Add results to output strings
            if all_results != "":
                output_strings[0] += "\n" + all_results
            output_strings[1] += total_known + "\n"
            
            # Add unknown mutations the total results of unknown mutations
            total_unknown_str += total_unknown + "\n"

            # Add drugs to druglist
            for drug in drug_list:
                unique_drug_list.append(drug.upper())

# Add unknown results to all results
if unknown_flag == True:
    output_strings[1] += "\n\nUnknown Mutations \n" + total_unknown_str 

# Make Resistance Prediction File

# Add header with all drug names
drug_lst = [drug for drug in drug_genes.keys()]

output_strings[2] = "Sample ID\t"
output_strings[2] += "\t".join(drug_lst) + "\n"

# Go throug all drugs in the database and see if prediction can be called.
pred_output = [sample_name]
for drug in drug_lst:

    # Check if resistance to drug was found
    if drug.upper() in unique_drug_list:
        pred_output.append("1")
    else:
        # Check at all genes associated with the drug resistance where found
        all_genes_found = True
        try:
            for gene in drug_genes[drug]:
                if GENES[gene]['found'] != 'yes' and GENES[gene]['found'] != 'partially':
                   all_genes_found = False
            if all_genes_found == False:
                pred_output.append("?")
            else:
                pred_output.append("0")
        except KeyError:
            pass

output_strings[2] += "\t".join(pred_output) + "\n"

# Write output files    
for i in range(len(output_files)):
    file_fullpath = out_path + "/" + sample_name + "_" + method + "_" + output_files[i]
    with open(file_fullpath, "w") as file_:
        file_.write(output_strings[i])

print(output_strings[1])
