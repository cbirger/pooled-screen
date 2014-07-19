import argparse

from ps_fastqToGct import fastqToGct

parser = argparse.ArgumentParser("PooledScreen.fastqToGct parser.")

parser.add_argument('ref', help='csv file containing gene barcodes and names')
parser.add_argument('cond', help='csv file containing condition barcodes and names')
parser.add_argument('reads', help='fastq file containing ngs reads')
parser.add_argument('gct', help='gct filename')
parser.add_argument('vector_marker', help='vector marker sequence that always precedes gene barcodes')
parser.add_argument('condition_barcode_len', help="length of condition barcode", type=int)
parser.add_argument('gene_barcode_match_len', help="length of gene barcode for matching", type=int)
parser.add_argument('max_mismatches', help='max number of mismatches allowed when matching sgRNA or shRNA sequence to gene barcode', type=int)
parser.add_argument('max_num_reads_to_process', help='max number of reads in fastq file to process', type=int)

args = parser.parse_args()

gene_barcodes_csv_file = args.ref
condition_barcodes_csv_file = args.cond
fastq_file = args.reads
gct_file = args.gct
vector_marker_sequence = args.vector_marker
condition_barcode_len = args.condition_barcode_len
gene_barcode_match_len = args.gene_barcode_match_len
max_mismatches = args.max_mismatches
max_num_reads_to_process = args.max_num_reads_to_process

print("gene_barcodes_csv_file: %s" % (gene_barcodes_csv_file))
print("condition_barcodes_csv_file: %s" % (condition_barcodes_csv_file))
print("fastq_file: %s" % (fastq_file))
print("gct_file: %s" % (gct_file))
print("vector_marker_sequence: %s" % (vector_marker_sequence))
print("condition_barcode_len: %i" % (condition_barcode_len))
print("gene_barcode_match_len: %i" % (gene_barcode_match_len))
print("max_mismatches: %i" % (max_mismatches))
print("max_num_reads_to_process: %i" % (max_num_reads_to_process)) 

fastqToGct(gene_barcodes_csv_file, condition_barcodes_csv_file, fastq_file, gct_file, vector_marker_sequence, condition_barcode_len, gene_barcode_match_len, max_mismatches, max_num_reads_to_process)
