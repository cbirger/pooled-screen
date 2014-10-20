import csv
import numpy as np
import time
from Bio import SeqIO
import re
from operator import itemgetter

#---------------------
#   Classes
#---------------------
class PooledScreenError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

class PooledScreen:
    """Base class for analysis of pooled screen data."""
    
    def __init__(self, conditions_file, reference_file, reads_file, construct_barcode_prefix, prefix_search_start_index, prefix_search_end_index,
                 construct_barcode_match_length, read_mismatches, max_num_reads_to_process, output_file_basename):

        self.conditions_file = conditions_file
        self.reference_file = reference_file
        self.reads_file = reads_file
        self.construct_barcode_prefix = construct_barcode_prefix
        self.prefix_search_start_index = prefix_search_start_index
        self.prefix_search_end_index = prefix_search_end_index
        self.construct_barcode_match_length = construct_barcode_match_length
        self.read_mismatches = read_mismatches
        self.max_num_reads_to_process = max_num_reads_to_process
        self.output_file_basename = output_file_basename

        self.__read_conditions_file()
        self.__read_reference_file()
        if self.construct_barcode_match_length > self.construct_barcode_length:
            raise PooledScreenError("condition barcode match length (%i) is greater than construct barcode length (%i)" % (self.construct_barcode_match_length, self.construct_barcode_length))
        read_counts_template = self.__create_read_counts_template()
        
        # initialize data structures for tallying exact and fuzzy read counts
        self.exact_match_read_counts_dict = {}
        self.fuzzy_match_read_counts_dict = {}
        for sample_barcode in self.sample_barcodes_list:
            self.exact_match_read_counts_dict[sample_barcode] = read_counts_template.copy()
        if self.read_mismatches > 0:
            for sample_barcode in self.sample_barcodes_list:
                self.fuzzy_match_read_counts_dict[sample_barcode] = read_counts_template.copy()

        # used for encoding sequence of 16 nucleotide bases as 64-bit integer 
        self.BASE_MAP = {'A':np.int64(1), 'T':np.int64(2), 'G':np.int64(4), 'C':np.int64(8)}

        # POPCOUNT_TABLE16 maps a 16-bit integer value (table index) to the
        # number of "on" bits in its binary representation
        self.POPCOUNT_TABLE16 = np.array([0] * 2**16)
        for index in xrange(len(self.POPCOUNT_TABLE16)):
            self.POPCOUNT_TABLE16[index] = (index & 1) + self.POPCOUNT_TABLE16[index >> 1]

        if self.read_mismatches > 0:
            self.__create_fuzzy_matching_array()

    def __read_conditions_file(self):

        self.conditions_dict = {}
        self.sample_barcodes_list = []
        self.sample_barcode_length = None
        with open(self.conditions_file, "rU") as csv_file:
            reader = csv.reader(csv_file)

            for row in reader:
                if row[0] in self.conditions_dict:
                    raise PooledScreenError("Duplicate sample barcode: (%s, %s)" % (row[0], row[1]))
        
                if self.sample_barcode_length is None:
                    self.sample_barcode_length = len(row[0])
                else:
                    if len(row[0]) != self.sample_barcode_length:
                        raise PooledScreenError("Sample barcode lengths differ: (%s, %s)" % (row[0], row[1]))

                self.conditions_dict[row[0]] = row[1]
                self.sample_barcodes_list.append(row[0])

    def __read_reference_file(self):

        self.constructs_dict = {}
        self.construct_barcodes_list = []
        self.construct_barcode_length = None

        with open(self.reference_file, "rU") as csv_file:
            reader = csv.reader(csv_file)

            for row in reader:
                if row[0] in self.constructs_dict:
                    raise PooledScreenError("Duplicate construct barcode: (%s, $s)" % (row[0], row[1]))
                
                if self.construct_barcode_length is None:
                    self.construct_barcode_length = len(row[0])
                else:
                    if len(row[0]) != self.construct_barcode_length:
                        raise PooledScreenError("Construct barcode lengths differ: (%s, %s)" % (row[0], row[1]))

                self.constructs_dict[row[0]] = row[1]
                self.construct_barcodes_list.append(row[0])

    def __create_read_counts_template(self):
        read_counts_template = {}
        for construct_barcode in self.construct_barcodes_list:
            construct_barcode_match_seq = construct_barcode[:self.construct_barcode_match_length]
            if construct_barcode_match_seq in read_counts_template:
                raise PooledScreenError("Shortened match sequence creating ambiguous mappings: %s" % construct_barcode_match_seq)
            else:
                read_counts_template[construct_barcode_match_seq] = 0

        return read_counts_template

    # This encodes a nucleotide sequence (16 bases or less) as a 64 bit integer.
    # Each nucleotide is encoded as a 4-bit value (see BASE_MAP).  The order
    # of encoded values in the 64 bit integer is reversed from the order of
    # the nucleotides in the sequence; e.g., the first nucelotide in the sequence
    # corresponds to the least signficant 4 bits.
    def __encode(self, seq):
        return sum([np.left_shift(self.BASE_MAP[x], 4*i) for i,x in enumerate(seq)])

    def __create_fuzzy_matching_array(self):
        self.fuzzy_matching_array = np.zeros((2, len(self.construct_barcodes_list)), dtype=np.int64)
        for i, construct_barcode_entry in enumerate(self.construct_barcodes_list):
            # note that construct_barcode_match_length is assumed to be no greater than 32 bases
            self.fuzzy_matching_array[0,i] = self.__encode(construct_barcode_entry[0:min(self.construct_barcode_match_length,16)])
            if self.construct_barcode_match_length > 16:
                self.fuzzy_matching_array[1,i] = self.__encode(construct_barcode_entry[16:self.construct_barcode_match_length])

    def __find_closest_construct_barcode(self, barcode_read_seq):
    
        barcode_read_seq_encoded_A = self.__encode(barcode_read_seq[0:16])
        x = np.bitwise_xor(barcode_read_seq_encoded_A, self.fuzzy_matching_array[0])
    
        y4 = self.POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]        
        y8 = self.POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 16), 0xffff)]
        y12 = self.POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 32), 0xffff)]
        y16 = self.POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 48), 0xffff)]
    
        distances = (y4 + y8 + y12 + y16) / 2

        if self.construct_barcode_match_length > 16:
            barcode_read_encoded_B = self.__encode(barcode_read_seq[16:self.construct_barcode_match_length])
            x = np.bitwise_xor(barcode_read_encoded_B, self.fuzzy_matching_array[1])
            y20 = self.POPCOUNT_TABLE16[np.bitwise_and(x, 0xffff)]            
            y24 = self.POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 16), 0xffff)]
            y28 = self.POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 32), 0xffff)]
            y32 = self.POPCOUNT_TABLE16[np.bitwise_and(np.right_shift(x, 48), 0xffff)]
            distances = distances + ((y20 + y24 + y28 + y32) / 2)
        
        min_distance = np.min(distances)
        return min_distance, np.sum(distances == min_distance), np.argmin(distances) 

    def tally_results(self):

        # Now ready to start parsing the fastq input file and count pairings of
        # condition and construct barcodes!

        start_time = time.time()

        self.read_count = 0
        self.short_read_count = 0

        self.known_sample_barcode_count = 0
        self.unknown_sample_barcode_count = 0
        self.unknown_sample_barcode_dict = {}

        self.construct_barcode_prefix_not_found_count = 0
        self.exact_matched_construct_barcode_count = 0
        self.fuzzy_matched_construct_barcode_count = 0
        self.unmatched_construct_barcode_count = 0

        min_read_length = self.prefix_search_end_index + len(self.construct_barcode_prefix) + self.construct_barcode_match_length

        reads_file_handle = open(self.reads_file, "rU")
        for read_record in SeqIO.parse(reads_file_handle, "fastq"):
            
            if self.max_num_reads_to_process > 0 and self.read_count == self.max_num_reads_to_process:
                break
    
            self.read_count +=1
            
            read_sequence = str(read_record.seq)
            if len(read_sequence) < min_read_length:
                self.short_read_count += 1
                continue

            # get read's sample barcode. 
            # NOTE:for now we assume is positioned at start of the read; we may want to relax
            #      this assumption in the future
            sample_barcode = read_sequence[0:self.sample_barcode_length]

            if sample_barcode in self.conditions_dict:
                self.known_sample_barcode_count += 1

                # look for construct barcode prefix
                matchObj = re.search(self.construct_barcode_prefix, read_sequence[self.prefix_search_start_index:self.prefix_search_end_index + len(self.construct_barcode_prefix)])
                if matchObj is not None:
                    # get construct barcode
                    construct_barcode_read = (read_sequence[self.prefix_search_start_index + matchObj.end():])[:self.construct_barcode_match_length]
                    if construct_barcode_read in self.exact_match_read_counts_dict[sample_barcode]:
                        self.exact_match_read_counts_dict[sample_barcode][construct_barcode_read] += 1
                        self.exact_matched_construct_barcode_count += 1
                    elif self.read_mismatches == 0:
                        self.unmatched_construct_barcode_count += 1
                    else:
                        # do fuzzy matching
                        closest_distance, num_closest, argmin_closest = self.__find_closest_construct_barcode(construct_barcode_read)
                        if closest_distance <= self.read_mismatches and num_closest == 1:
                            self.fuzzy_match_read_counts_dict[sample_barcode][self.construct_barcodes_list[argmin_closest][0:self.construct_barcode_match_length]] += 1
                            self.fuzzy_matched_construct_barcode_count += 1
                        else:
                            self.unmatched_construct_barcode_count += 1
                else:
                    # did not find construct_barcode_prefix
                    self.construct_barcode_prefix_not_found_count += 1
            else:
                # unrecognized sample barcode
                self.unknown_sample_barcode_count += 1
                if sample_barcode not in self.unknown_sample_barcode_dict:
                    self.unknown_sample_barcode_dict[sample_barcode] = 1
                else:
                    self.unknown_sample_barcode_dict[sample_barcode] += 1
        
        reads_file_handle.close()
        
        finish_time = time.time()
        elapsed_time_seconds = finish_time - start_time
        elapsed_time_minutes = elapsed_time_seconds/float(60)

        print("elapsed time (sec) = %i" % elapsed_time_seconds)
        print("elapsed time (min) = %i\n" % elapsed_time_minutes)

    def __write_gct_file(self, gct_filename, read_counts_dict_list):
    
        # create gct file from one or more read_counts_dictionaries
        handle = open(gct_filename, "w")

        # first two rows of gct file
        handle.write("#1.2\n")
        handle.write("%i\t%i\n" % (len(self.construct_barcodes_list), len(self.sample_barcodes_list)))

        # header row of gct file
        header_row ="Name\tDescription"
        for sample_barcode in self.sample_barcodes_list:
            header_row += "\t%s" % self.conditions_dict[sample_barcode]
        handle.write("%s\n" % header_row)

        for construct_barcode in self.construct_barcodes_list:
            results_row = "%s\t%s" % (construct_barcode, self.constructs_dict[construct_barcode])
            for sample_barcode in self.sample_barcodes_list:
                read_count = sum([x[sample_barcode][construct_barcode[0:self.construct_barcode_match_length]] for x in read_counts_dict_list])
                results_row += "\t%i" % read_count
            handle.write("%s\n" % results_row)

        handle.close()

    def write_gct_files(self):

        if self.read_mismatches == 0:
            self.__write_gct_file(self.output_file_basename + ".gct", [self.exact_match_read_counts_dict])
        else:
            self.__write_gct_file(self.output_file_basename + "_EXACT.gct", [self.exact_match_read_counts_dict])
            self.__write_gct_file(self.output_file_basename + "_FUZZY.gct", [self.fuzzy_match_read_counts_dict])
            self.__write_gct_file(self.output_file_basename + ".gct", [self.exact_match_read_counts_dict, self.fuzzy_match_read_counts_dict])
        
    def write_unknown_sample_barcodes_file(self):
        if self.unknown_sample_barcode_count > 0:
            unknown_sample_barcode_list = self.unknown_sample_barcode_dict.items()
            unknown_sample_barcode_list_sorted = sorted(unknown_sample_barcode_list, key=itemgetter(1), reverse=True)
            handle = open(self.output_file_basename + "_UNKNOWN_SAMPLE_BARCODES.txt", 'w')
            writer = csv.writer(handle, delimiter = '\t')
            for x in unknown_sample_barcode_list_sorted:
                writer.writerow(x)
            handle.close()

    def write_summary_file(self):
    
        handle = open(self.output_file_basename + "_SUMMARY.txt", 'w')

        #input parameters
        handle.write("INPUT PARAMETERS\n\n")
        handle.write("Conditions file: %s\n" % self.conditions_file)
        handle.write("Reference file: %s\n" % self.reference_file)
        handle.write("Reads file: %s\n" % self.reads_file)
        handle.write("Construct barcode prefix: %s\n" % self.construct_barcode_prefix)
        handle.write("Prefix search start index: %i\n" % self.prefix_search_start_index)
        handle.write("Prefix search end index: %i\n" % self.prefix_search_end_index)
        handle.write("Construct barcode match length: %i\n" % self.construct_barcode_match_length)
        handle.write("Read mismatches: %i\n" % self.read_mismatches)
        handle.write("Max num reads to process: %i\n" % self.max_num_reads_to_process)
        handle.write("Output file basename: %s\n\n" % self.output_file_basename)

        #results
        handle.write("SUMMARY STATISTICS\n\n")
        handle.write("Total number of reads: %i\n" % self.read_count)
        successfully_matched_read_count = self.exact_matched_construct_barcode_count + self.fuzzy_matched_construct_barcode_count
        handle.write("Total number of successfully counted reads: %i\n" % successfully_matched_read_count)
        successfully_matched_percentage = 100 * successfully_matched_read_count/float(self.read_count)
        handle.write("Percentage of reads successfully matched: %4.2f%%\n\n" % successfully_matched_percentage)

        # this block shouldn't be run read_mismatches == 0
        if self.read_mismatches > 0:
            handle.write("Total number of successfully counted reads exactly matching construct barcode: %i\n" % self.exact_matched_construct_barcode_count)
            handle.write("Total number of successfully counted reads fuzzy matching construct barcode: %i\n" % self.fuzzy_matched_construct_barcode_count)
            exact_matched_percentage = 100 * self.exact_matched_construct_barcode_count/float(successfully_matched_read_count)
            handle.write("Percentage of successfully counted reads exactly matching construct barcode: %4.2f%%\n\n" % exact_matched_percentage)

        handle.write("Total number of reads with unknown sample barcodes: %i\n" % self.unknown_sample_barcode_count)
        handle.write("Total number of reads with known sample barcodes and no construct barcode prefix: %i\n" % self.construct_barcode_prefix_not_found_count)
        handle.write("Total number of reads with known sample barcodes and no matching construct barcode: %i\n" % self.unmatched_construct_barcode_count)

        handle.close()


