import csv
import time
import re
from Bio import SeqIO
import numpy as np

#---------------------
#   Constants
#---------------------
BASE_MAP = {'A':np.int64(1), 'T':np.int64(2), 'G':np.int64(4), 'C':np.int64(8)}
UNKNOWN_CONDITION_BARCODE_CSV_FILE = "unknown_condition_barcode.csv"


POPCOUNT_TABLE16 = np.array([0] * 2**16)
for index in xrange(len(POPCOUNT_TABLE16)):
    POPCOUNT_TABLE16[index] = (index & 1) + POPCOUNT_TABLE16[index >> 1]

#---------------------
#   Functions
#---------------------

# This encodes a nucleotide sequence (16 bases or less) as a 64 bit integer.
# Each nucleotide is encoded as a 4-bit value (see BASE_MAP).  The order
# of encoded values in the 64 bit integer is reversed from the order of
# the nucleotides in the sequence; e.g., the first nucelotide in the sequence
# corresponds to the least signficant 4 bits.
def encode(seq):
    seq_len = len(seq)
    return sum([np.left_shift(BASE_MAP[x], 4*i) for i,x in enumerate(seq)])

# This function employs a vectorized popcount operation to determine the gene barcode
# that has the minimum number of mismatches with the sgRNA sequence.  It returns
# the minimum mismatch count and the index of the matching sequence.  Note that the index is that
# of the first of potentially multiple gene barcodes that have the minimum number of mismatches. 
def find_min_mismatches(sgRNA_sequence, gene_barcode_match_len, gene_barcode_array):
    # This can be further optimized if necessary
    sgRNA_search_seq_A = encode(sgRNA_sequence[0:min(gene_barcode_match_len, 16)])
    x = np.bitwise_xor(sgRNA_search_seq_A, gene_barcode_array[0])

    if gene_barcode_match_len <= 4:
        y4 = POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]
        return np.min(y4)/2, np.argmin(y4)
    elif gene_barcode_match_len <= 8:
        y4 = POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]        
        y8 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 16), 0xffff)]
        y = y4 + y8
        return np.min(y)/2, np.argmin(y)
    elif gene_barcode_match_len <= 12:
        y4 = POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]        
        y8 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 16), 0xffff)]
        y12 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 32), 0xffff)]
        y = y4 + y8 + y12
        return np.min(y)/2, np.argmin(y)
    elif gene_barcode_match_len <= 16:
        y4 = POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]        
        y8 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 16), 0xffff)]
        y12 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 32), 0xffff)]
        y16 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 48), 0xffff)]
        y = y4 + y8 + y12 + y16
        return np.min(y)/2, np.argmin(y)
    else:
        y4 = POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]        
        y8 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 16), 0xffff)]
        y12 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 32), 0xffff)]
        y16 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 48), 0xffff)]
        sgRNA_search_seq_B = encode(sgRNA_sequence[16:gene_barcode_match_len])
        x = np.bitwise_xor(sgRNA_search_seq_B, gene_barcode_array[1])
        if gene_barcode_match_len <= 20:
            y20 = POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]
            y = y4 + y8 + y12 + y16 + y20
            return np.min(y)/2, np.argmin(y)
        elif gene_barcode_match_len <= 24:
            y20 = POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]            
            y24 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 16), 0xffff)]
            y = y4 + y8 + y12 + y16 + y20 + y24
            return np.min(y)/2, np.argmin(y)
        elif gene_barcode_match_len <= 28:
            y20 = POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]            
            y24 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 16), 0xffff)]
            y28 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 32), 0xffff)]
            y = y4 + y8 + y12 + y16 + y20 + y24 + y28
            return np.min(y)/2, np.argmin(y)
        else:
            y20 = POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]            
            y24 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 16), 0xffff)]
            y28 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 32), 0xffff)]
            y32 = POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 48), 0xffff)]
            y = y4 + y8 + y12 + y16 + y20 + y24 + y28 + y32
            return np.min(y)/2, np.argmin(y)

    assert(False)

def fastqToGct(gene_barcodes_csv_file, condition_barcodes_csv_file, fastq_file, gct_file, vector_marker_sequence, condition_barcode_len, gene_barcode_match_len, max_mismatches, max_num_reads_to_process):
    # first create a dictionary template whose keys are the 
    # gene barcodes
    
    gene_barcode_dictionary_template = {}
    gene_barcode_list = []
    short_gene_barcode_count = 0
    with open(gene_barcodes_csv_file, "rU") as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            if len(row[0]) < gene_barcode_match_len:
                short_gene_barcode_count += 1
                continue
            gene_barcode_dictionary_template[row[0][0:gene_barcode_match_len]] = 0
            gene_barcode_list.append((row[0], row[1]))

    gene_barcode_array = np.zeros((2,len(gene_barcode_list)), dtype=np.int64)
    for i, gene_barcode_entry in enumerate(gene_barcode_list):
        # note that gene barcodes are assumed to be no longer than 32 nucleotides
        gene_barcode_array[0,i] = encode(gene_barcode_entry[0][0:min(gene_barcode_match_len,16)])
        if gene_barcode_match_len > 16:
            gene_barcode_array[1,i] = encode(gene_barcode_entry[0][16:gene_barcode_match_len])


    # initialize a read counts dictionary whose keys are the 
    # condition barcodes and values are 
    # clones of gene_barcode_dictionary_template
    read_counts_dictionary = {}
    condition_barcode_list = []
    with open(condition_barcodes_csv_file, "rU") as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            read_counts_dictionary[row[0]] = gene_barcode_dictionary_template.copy()
            condition_barcode_list.append((row[0], row[1]))

    fastq_file_handle = open(fastq_file, "rU")

    read_count = 0
    unknown_condition_barcode_count = 0
    unknown_condition_barcode_dictionary = {}
    known_condition_barcode_count = 0
    vector_marker_not_present_count = 0
    vector_marker_invalid_position_count = 0
    exact_gene_barcode_match_count = 0
    inexact_gene_barcode_match_count = 0
    unknown_gene_barcode_count = 0
    
    start_time = time.time()

    for read_record in SeqIO.parse(fastq_file_handle, "fastq"):
        if max_num_reads_to_process > 0 and read_count == max_num_reads_to_process:
            break

        read_count +=1

        # get read's condition barcode.  Barcode assumed to be at the start of the read and
        # 6 nucleotides long.
        # NOTE: we may need to relax assumptions regarding the length and position of condition
        # barcode
        condition_barcode = str(read_record.seq[0:condition_barcode_len])
        if condition_barcode not in read_counts_dictionary:
            unknown_condition_barcode_count += 1
            if condition_barcode not in unknown_condition_barcode_dictionary:
                unknown_condition_barcode_dictionary[condition_barcode] = 1
            else:
                unknown_condition_barcode_dictionary[condition_barcode] += 1
        else:
            known_condition_barcode_count += 1

            # get sgRNA/shRNA sequence from read
            read_sequence = str(read_record.seq[condition_barcode_len:])
            # look for vector marker sequence
            ends = [m.end() for m in re.finditer(vector_marker_sequence, read_sequence)]
            if len(ends) == 0:
                vector_marker_not_present_count += 1
            elif ends[0] + gene_barcode_match_len > len(read_sequence):
                vector_marker_invalid_position_count +=1
            else:
                # first look for exact match (across first gene_barcode_match_len bases)
                sgRNA_sequence = read_sequence[ends[0]:ends[0] + gene_barcode_match_len]
                if len(sgRNA_sequence) < gene_barcode_match_len:
                    # drop truncated reads that do not contain sufficent bases to
                    # conduct gene barcode search 
                    unknown_gene_barcode_count += 1
                elif sgRNA_sequence in read_counts_dictionary[condition_barcode]:
                    # found an exact match
                    read_counts_dictionary[condition_barcode][sgRNA_sequence] += 1
                    exact_gene_barcode_match_count += 1
                elif max_mismatches == 0:
                    # drop read if not doing mismatch-tolerant search
                    unknown_gene_barcode_count += 1
                else:
                    # here we search for a gene barcode that has up to a small number (max_mismatches) of
                    # base mismatches with sgRNA_sequence; we DO NOT allow indels.
                    # IMPORTANT: this current implementation can hopefully be made faster if we
                    # deterimine accounting for reads with "fuzzy" matching impacts the results
                    min_mismatches, argmin_mismatches = find_min_mismatches(sgRNA_sequence, gene_barcode_match_len, gene_barcode_array)

                    if min_mismatches <= max_mismatches:
                        read_counts_dictionary[condition_barcode][gene_barcode_list[argmin_mismatches][0][0:gene_barcode_match_len]] +=1
                        inexact_gene_barcode_match_count += 1
                    else:
                        unknown_gene_barcode_count += 1
                    

    fastq_file_handle.close()
    finish_time = time.time()
    elapsed_time_seconds = finish_time - start_time
    elapsed_time_minutes = elapsed_time_seconds/float(60)
    print("elapsed time (sec) = %i" % elapsed_time_seconds)
    print("elapsed time (min) = %i\n" % elapsed_time_minutes)

    print("short_gene_barcode_count: %i\n" % (short_gene_barcode_count))

    print("total number of reads: %i" % (read_count))
    print("number of reads with known barcode: %i" % (known_condition_barcode_count))
    print("number of reads with unknown barcode: %i" % (unknown_condition_barcode_count))
    percent_reads_with_unknown_barcode = 100 * unknown_condition_barcode_count/float(read_count)
    print("percentage of reads dropped due to unrecognized barcode: %4.2f%%" % (percent_reads_with_unknown_barcode))
    print("number of reads with no vector marker: %i" % (vector_marker_not_present_count))
    print("number of reads with misplaced vector marker: %i" % (vector_marker_invalid_position_count))
    percent_reads_missing_valid_vector_marker = 100 * (vector_marker_not_present_count + vector_marker_invalid_position_count)/float(read_count)
    print("percentage of reads dropped due to missing or misplaced vector marker: %4.2f%%" % (percent_reads_missing_valid_vector_marker))
    print("number of reads with exact gene barcode match: %i" % (exact_gene_barcode_match_count))
    if max_mismatches != 0:
        print("number of reads with inexact (<= %i mismatches) gene barcode match: %i" % (max_mismatches,inexact_gene_barcode_match_count))
    print("number of reads with unknown gene barcode: %i" % (unknown_gene_barcode_count))
    percent_reads_with_unknown_gene_barcode = 100 * unknown_gene_barcode_count/float(read_count)
    print("percentage of reads dropped due to unknown gene barcode: %4.2f%%" % (percent_reads_with_unknown_gene_barcode))
    
    # create gct file from read_counts_dictionary

    handle = open(gct_file, "w")

    # first two rows of gct file
    handle.write("#1.2\n")
    handle.write("%i\t%i\n" % (len(gene_barcode_list), len(condition_barcode_list)))

    # header row of gct file
    header_row ="Name\tDescription"
    for condition_barcode in iter(condition_barcode_list):
        header_row += "\t%s" % condition_barcode[1]
    handle.write("%s\n" % header_row)

    for gene_barcode in iter(gene_barcode_list):
        results_row = "%s\t%s" % (gene_barcode[0], gene_barcode[1])
        for condition_barcode in iter(condition_barcode_list):
            results_row += "\t%i" % read_counts_dictionary[condition_barcode[0]][gene_barcode[0][0:gene_barcode_match_len]]
        handle.write("%s\n" % results_row)

    handle.close()

    print("GCT file written to %s" % (gct_file))

    if unknown_condition_barcode_count > 0:
        with open(UNKNOWN_CONDITION_BARCODE_CSV_FILE, 'wb') as csvfile:
            file_writer = csv.writer(csvfile)
            for barcode in sorted(unknown_condition_barcode_dictionary, key=unknown_condition_barcode_dictionary.get, reverse=True):
                file_writer.writerow([barcode, unknown_condition_barcode_dictionary[barcode]])
