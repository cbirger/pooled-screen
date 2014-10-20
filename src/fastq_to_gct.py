import argparse
import pooled_screen as ps
import sys

#---------------------
#   Functions
#---------------------

def fastqToGct(conditions_file, reference_file, reads_file,
               construct_barcode_prefix, prefix_search_start_index, prefix_search_end_index,
               construct_barcode_match_length, read_mismatches,
               max_num_reads_to_process, output_file_basename):

    ps_instance = ps.PooledScreen(conditions_file, reference_file, reads_file,
                                  construct_barcode_prefix, prefix_search_start_index, prefix_search_end_index,
                                  construct_barcode_match_length, read_mismatches,
                                  max_num_reads_to_process, output_file_basename)

    ps_instance.tally_results()
    ps_instance.write_gct_files()
    ps_instance.write_unknown_sample_barcodes_file()
    ps_instance.write_summary_file()
    
#---------------------
#   Script
#---------------------

# to do - error handling

parser = argparse.ArgumentParser("PooledScreen.fastqToGct parser.")

parser.add_argument('cond')
parser.add_argument('ref')
parser.add_argument('reads')
parser.add_argument('prefix')
parser.add_argument('start_index', type=int)
parser.add_argument('end_index', type=int)
parser.add_argument('--match_length', dest='match_length', type=int)
parser.add_argument('read_mismatches', type=int)
parser.add_argument('max_num_reads', type=int)
parser.add_argument('basename')

args = parser.parse_args()

conditions_file = args.cond
reference_file = args.ref
reads_file = args.reads
construct_barcode_prefix = args.prefix
prefix_search_start_index = args.start_index
prefix_search_end_index = args.end_index
if args.match_length is not None:
    construct_barcode_match_length = args.match_length
else:
    construct_barcode_match_length = 21
read_mismatches = args.read_mismatches
max_num_reads_to_process = args.max_num_reads
output_file_basename = args.basename


try:
    fastqToGct(conditions_file, reference_file, reads_file,
               construct_barcode_prefix, prefix_search_start_index, prefix_search_end_index,
               construct_barcode_match_length, read_mismatches,
               max_num_reads_to_process, output_file_basename)
except ps.PooledScreenError as e:
    print "PooledScreenError: %s" % e
    sys.exit(1)
except:
    raise
