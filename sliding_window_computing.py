import sys
from Bio import SeqIO


def sliding_windows_init(genome_fasta, window_size, window_overlap, output_file, chr_arg):

	window_size = int(window_size)
	window_overlap = int(window_overlap)

	print("Parsing genome...")
	genome_sizes = {record.id: len(str(record.seq)) for record in SeqIO.parse(genome_fasta, "fasta")}
	sliding_windows = {}

	print("Computing sliding windows...")
	for chr in genome_sizes:
		if chr == chr_arg:
			sliding_windows[chr] = compute_sliding_windows(genome_sizes[chr], window_size, window_overlap)
			print(("Chromosome " + chr + " done"))
	output_handler = open(output_file, "w")


	print("Writing results...")
	for chr in sliding_windows:

		window_id_counter = 0
		for sliding_window in sliding_windows[chr]:

			window_id = chr + "-w" + str(window_id_counter)
			to_write = [chr,window_id, str(sliding_window[0]), str(sliding_window[1])]

			output_handler.write("\t".join(to_write) + "\n")

			window_id_counter += 1




def compute_sliding_windows(chr_size, window_size, window_overlap):


    sliding_windows_list = []
    for i in range(0,chr_size,window_overlap):
        sliding_windows_list.append([i, i+window_size])

    sliding_windows_list[0][0] = 1
    sliding_windows_list[-1][1] = chr_size

    return sliding_windows_list

#sliding_windows_init(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
