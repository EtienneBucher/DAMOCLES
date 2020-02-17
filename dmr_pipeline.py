import sys
sys.dont_write_bytecode = True

import os
from shutil import copyfile

# Modules
import sliding_window_computing
#import merge_intervals
import diff_methylation_test_2
import time
import resource

def launch_dmr_pipeline(genome_fasta, window_size, window_overlap, interval_merging_limit, sample_conf_file, output_directory, chr_arg, gff_file):
    print("Starting DMR analysis, copying files...")
    #if not os.path.exists(output_directory):
    #    os.makedirs(output_directory)
    #    os.makedirs(output_directory + "/data")
    #    os.makedirs(output_directory + "/results")

    copyfile(genome_fasta, output_directory + "/data/genome.fasta")
    copyfile(sample_conf_file, output_directory + "/data/sample.conf")


    print("Sliding windows computing...")
    sliding_window_computing.sliding_windows_init(output_directory + "/data/genome.fasta", window_size, window_overlap, output_directory + "/data/regions_to_test.conf", chr_arg)

    start_time = time.time()

    print("Testing for DMRs...")
    diff_methylation_test_2.diff_methylation_test(sample_conf_file, output_directory + "/data/regions_to_test.conf", output_directory + "/data/genome.fasta", output_directory + "/results/dmrs_not_merged", chr_arg, "CG,CHG,CHH")

    print(("--- %s seconds ---" % (time.time() - start_time)))
    print("Peak memory usage:", "{0:0.1f}".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*0.0000009536743164),"GB")

    #print("Merging intervals...")
    #merge_intervals.merge_intervals(output_directory + "/results/dmrs_not_merged.sign.gff3", output_directory + "/results/dmrs_merged.sign.gff3", interval_merging_limit)


launch_dmr_pipeline(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
