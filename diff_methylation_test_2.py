# -*- coding: utf-8 -*-
import sys
from Bio import SeqIO
import statsmodels
from statsmodels.sandbox.stats.runs import mcnemar
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy
import scipy
import psutil
import resource
import warnings
from joblib import Parallel, delayed
import os


# Context to test :
# CG CHG CHH CGxCHG CGxCHH CHGxCHH CGxCHGxCHH all custom
def diff_methylation_test(sample_conf_file, regions_file, genome_fasta, output_file, chr_arg, contexts = "all"):

    print("Init genome dictionnary...")
    genome_dict = init_genome(genome_fasta, chr_arg)

    print("Parsing conf file...")
    (sample_conf_data, samples_names) = parse_sample_conf_file(sample_conf_file)

    print("Parsing regions_to_test...")
    regions_to_test = parse_regions_file(regions_file, contexts, chr_arg)

    print("Parsing methylation files...")
    genome_dict = fill_genome_dict(genome_dict, sample_conf_data, chr_arg)

    print("Performing statistical tests...")
    regions_stats_list = perform_statistical_tests(genome_dict, regions_to_test, samples_names)

    print("Writing results...")
    write_results(regions_stats_list, output_file)
    print("Max memory when finished:", "{0:0.1f}".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*0.0000009536743164),"GB")


def write_results(regions_stats_list, output_file):

    output_handler = open(output_file + ".all.gff3", "w")
    output_handler_sign = open(output_file + ".sign.gff3", "w")

    header = ["region_id", "context", "std_replicates", "start", "end", "p-value", ".", "methyl_diff", "info"]
    output_handler.write("\t".join(header) + "\n")
    output_handler_sign.write("\t".join(header) + "\n")

    for r in regions_stats_list:
        #to_write = [region_stats["region_id"],region_stats["context"], str(region_stats["pvalue"]), str(region_stats["non_consistent_percentage"]), str(region_stats["number_tested_cytosines"]), region_stats["over"], str(region_stats["contingency_table"])]

        info=["ID=" + r["region_id"], "context=" + r["context"], "avg_of_std=" + str(r["non_consistent_percentage"]), "number_tested_cytosines=" +  str(r["number_tested_cytosines"]), "overmethylated_in=" + r["over"], "cov_sd=" + str(r["cov_sd"]), "values=" + str(r["contingency_table"]), "std_of_avg=" + str(r["std_of_avg"]), "coverages=" + str(r["indiv_coverages"])]
        to_write = [r["chr"],r["context"], str(r["non_consistent_percentage"]) + "!!!" + str(r["std_of_avg"]), str(r["start"]), str(r["end"]), str(r["pvalue"]), "+", str(r["methyl_diff"]), ";".join(info)]
        #to_write2 = [r["chr"],r["context"],r["region_id"],str(r["start"]), str(r["end"]),str(r["methyl_diff"])] TO BE DONE FOR BETTER OUTPUT
        output_handler.write("\t".join(to_write) + "\n")


        if r["pvalue"] <= 0.05:
            output_handler_sign.write("\t".join(to_write) + "\n")

def average_replicates(cytosines_to_test):

    for cytosine in cytosines_to_test:
        for sample in cytosine["ratios"]:
            if type(cytosine["ratios"][sample]) is list:

                nominateur = 0
                denominateur = 0

                sd_ratio_calc = []
                sd_cov_calc = []


                for tup in cytosine["ratios"][sample]:

                    nominateur += tup[0]
                    denominateur += tup[1]
                    sd_ratio_calc.append(float(tup[0]) / float(tup[1]))
                    sd_cov_calc.append(int(tup[1]))

                ratio = float(nominateur) / float(denominateur)
                cytosine["ratios"][sample] = ratio

                ratio_sd = numpy.std(sd_ratio_calc)
                cov_sd = numpy.std(sd_cov_calc)

                if "sd" not in cytosine:
                    cytosine["sd"] = {}
                cytosine["sd"][sample] = ratio_sd

                if "cov_sd" not in cytosine:
                    cytosine["cov_sd"] = {}
                cytosine["cov_sd"][sample] = cov_sd

    sd = []
    sd_dict = {}
    for cytosine in cytosines_to_test:
        for sample in cytosine["sd"]:
            if sample not in sd_dict:
                sd_dict[sample] = []
            sd_dict[sample].append(cytosine["sd"][sample])

    for sample in sd_dict:
        sd_dict[sample] = numpy.mean(sd_dict[sample])
        sd.append(sample + ":" + str(sd_dict[sample]))

    sd = "$".join(sd)
    if sd == "":
        sd = "NA"

    #-------------------------------------------
    cov_sd = []
    cov_sd_dict = {}
    for cytosine in cytosines_to_test:
        for sample in cytosine["cov_sd"]:
            if sample not in cov_sd_dict:
                cov_sd_dict[sample] = []
            cov_sd_dict[sample].append(cytosine["cov_sd"][sample])

    for sample in cov_sd_dict:
        cov_sd_dict[sample] = numpy.mean(cov_sd_dict[sample])
        cov_sd.append(sample + ":" + str(cov_sd_dict[sample]))

    cov_sd = "$".join(cov_sd)
    if cov_sd == "":
        cov_sd = "NA"

    return (cytosines_to_test, sd, cov_sd)

def average_replicates_simple(cytosines_to_test):

    for cytosine in cytosines_to_test:
        for sample in cytosine["ratios"]:
            if type(cytosine["ratios"][sample]) is list:

                nominateur = 0
                denominateur = 0

                for tup in cytosine["ratios"][sample]:

                    nominateur += tup[0]
                    denominateur += tup[1]

                ratio = float(nominateur) / float(denominateur)
                cytosine["ratios"][sample] = ratio

    return cytosines_to_test


def build_lists_to_compare(cytosines_to_test, samples_names, context):

    lists_to_compare = [[],[]]

    cytosines_to_test = [cytosine for cytosine in cytosines_to_test if cytosine["context"] in context]


    sample_1_name = samples_names[0]
    sample_2_name = samples_names[1]

    number_tested_cytosines = 0

    for cytosine in cytosines_to_test:
        if len(cytosine["ratios"]) == 2:
            lists_to_compare[0].append(cytosine["ratios"][sample_1_name])
            lists_to_compare[1].append(cytosine["ratios"][sample_2_name])
            number_tested_cytosines += 1

    if numpy.mean(lists_to_compare[0]) > numpy.mean(lists_to_compare[1]):
        overmethylated_sample = sample_1_name
    elif numpy.mean(lists_to_compare[0]) < numpy.mean(lists_to_compare[1]):
        overmethylated_sample = sample_2_name
    else:
        overmethylated_sample = "equal_methylation_in_both_samples"


    return (lists_to_compare, number_tested_cytosines, overmethylated_sample)




def perform_one_test(cytosines_to_test, context, samples_names):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        (contingency_table, number_tested_cytosines, overmethylated_sample) = build_lists_to_compare(cytosines_to_test, samples_names, context)

        try:
            #stat_test = statsmodels.sandbox.stats.runs.mcnemar(contingency_table, y=None, exact=True, correction=True)
            stat_test = scipy.stats.wilcoxon(x=contingency_table[0], y=contingency_table[1], zero_method='wilcox', correction=False)
            pvalue = float(stat_test[1])
            contingency_table[0] = numpy.mean(contingency_table[0])
            contingency_table[1] = numpy.mean(contingency_table[1])

            methyl_diff = abs(contingency_table[0] - contingency_table[1])
        except:
            pvalue = 1

            methyl_diff = 0

    #region_id = chr + "_" + str(start) + ":" + str(end)
        test_object = {"contingency_table": contingency_table, "methyl_diff" : methyl_diff, "over": overmethylated_sample, "context": context, "pvalue":pvalue, "number_tested_cytosines": number_tested_cytosines}

    return test_object

import copy
def test_one_region(region_to_test, samples_names, cytosines_to_test, chr):

    contexts_to_test = region_to_test["contexts_to_test"]
    #print(region_to_test)

    start = region_to_test["start"]
    end = region_to_test["end"]

    cytosines_to_test = copy.deepcopy(cytosines_to_test)
    cytosines_for_std_of_avg = copy.deepcopy(cytosines_to_test)

    cytosines_to_test = [cytosine for cytosine in cytosines_to_test if cytosine != -1]
    cytosines_for_std_of_avg = [cytosine for cytosine in cytosines_for_std_of_avg if cytosine != -1]

    #for c in cytosines_to_test:
    #    print(c["ratios"])


    (cytosines_to_test, non_consistent_percentage, cov_sd) = average_replicates(cytosines_to_test)
    #(cytosines_to_test, non_consistent_percentage, cov_sd) = (average_replicates_simple(cytosines_to_test), "NA", "NA")

    test_objects = []
    for context in contexts_to_test:

        test_object = perform_one_test(cytosines_to_test, context, samples_names)
        test_object["region_id"] = region_to_test["name"] + "-" + context
        test_object["non_consistent_percentage"] = non_consistent_percentage
        test_object["cov_sd"] = cov_sd
        test_object["chr"] = chr
        test_object["start"] = start
        test_object["end"] = end

	#print(chr)
        cytosines_for_other = [cytosine for cytosine in cytosines_for_std_of_avg if cytosine["context"] in context]

        res = compute_std_of_avg(cytosines_for_other)
        test_object["std_of_avg"] = res[0]
        test_object["indiv_coverages"] = res[1]

        test_objects.append(test_object)

    return test_objects



def compute_std_of_avg(cytosines_for_std_of_avg):

    replicates_values = {}
    replicates_cov = {}
    #print(cytosines_for_std_of_avg)

    for cytosine in cytosines_for_std_of_avg:



        for sample in cytosine["ratios"]:

            if type(cytosine["ratios"][sample]) is list:

                if sample not in replicates_values:
                    replicates_values[sample] = [[] for elt in cytosine["ratios"][sample]]
                    replicates_cov[sample] = [[] for elt in cytosine["ratios"][sample]]
                else:
                    if len(cytosine["ratios"][sample]) > len(replicates_values[sample]):
                        diff = len(cytosine["ratios"][sample]) - len(replicates_values[sample])
                        for i in range(diff):
                            replicates_values[sample].append([])
                            replicates_cov[sample].append([])




                i = 0
                for tup in cytosine["ratios"][sample]:


                    ratio = float(tup[0]) / float(tup[1])
                    replicates_values[sample][i].append(ratio)
                    replicates_cov[sample][i].append(tup[1])
                    i += 1



    to_compute = {}

    for sample in replicates_values:

        if sample not in to_compute:
            to_compute[sample] = []

        for replicate in replicates_values[sample]:
            replicate = numpy.mean(replicate)
            to_compute[sample].append(replicate)


    std_of_avg = []
    for sample in to_compute:

        std_sample = numpy.std(to_compute[sample])
        std_of_avg.append(sample + ":" + str(std_sample))

    std_of_avg = "$".join(std_of_avg)
    if std_of_avg == "":
        std_of_avg = "NA"

    replicates_cov_write = {}
    for sample in replicates_cov:
        if sample not in replicates_cov_write:
            replicates_cov_write[sample] = []
        for replicate in replicates_cov[sample]:
            replicate = numpy.mean(replicate)
            replicates_cov_write[sample].append(replicate)


    return (std_of_avg, replicates_cov_write)

def perform_statistical_tests(genome_dict, regions_to_test, samples_names):

    regions_stats_list = []

    regions_count = sum([len(regions_to_test[chr]) for chr in regions_to_test])

    done_count = 0
    for chr in regions_to_test:
        for r in regions_to_test[chr]:
            if genome_dict[chr][r["start"]-1:r["end"]] != []:
                regions_stats_list.extend(test_one_region(r, samples_names, genome_dict[chr][r["start"]-1:r["end"]], chr))


    print("regions_stats_list : " + str(len(regions_stats_list)))
        #    done_count += 1
        #    done_pc = float(done_count) / float(regions_count) * 100.0
        #    if done_count % 1000 == 0:
        #        print(str(done_pc) + " % done")

        #regions_stats_list = Parallel(n_jobs=1 ,verbose=100, backend="threading")(delayed(test_one_region)(region_to_test, samples_names, genome_dict, chr) for region_to_test in regions_to_test[chr])
        #regions_stats_list = Parallel(n_jobs=5 ,verbose=0, backend="threading")(delayed(test_one_region)(r, samples_names, genome_dict[chr][r["start"]-1:r["end"]], chr) for r in regions_to_test[chr])


    # test lists to dict to to multiple testing correction on subs
    regions_stats_dict = {}
    for test in regions_stats_list:
        number_cytosines = "cyto" + str(test["number_tested_cytosines"])
        test_context = test["context"]

        if number_cytosines not in regions_stats_dict:
            regions_stats_dict[number_cytosines] = {}

        if test_context not in regions_stats_dict[number_cytosines]:
            regions_stats_dict[number_cytosines][test_context] = []

        regions_stats_dict[number_cytosines][test_context].append(test)


    for number_cytosines in regions_stats_dict:
        for test_context in regions_stats_dict[number_cytosines]:
            pvalues = [elt["pvalue"] for elt in regions_stats_dict[number_cytosines][test_context]]
            pvalues = [1 if numpy.isnan(elt) else elt for elt in pvalues]
            corrected_pvalues = statsmodels.sandbox.stats.multicomp.multipletests(pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]

            i = 0
            for elt in corrected_pvalues:
                regions_stats_dict[number_cytosines][test_context][i]["pvalue"] = elt
                i += 1



    #dict_to_list
    list_to_return = []
    for number_cytosines in regions_stats_dict:
        for test_context in regions_stats_dict[number_cytosines]:
            for test in regions_stats_dict[number_cytosines][test_context]:
                list_to_return.append(test)
    #pvalues = [elt["pvalue"] for elt in regions_stats_list]
    #pvalues = [1 if numpy.isnan(elt) else elt for elt in pvalues]

    #corrected_pvalues = statsmodels.sandbox.stats.multicomp.multipletests(pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]

    #i = 0
    #for elt in corrected_pvalues:
    #    regions_stats_list[i]["pvalue"] = elt
    #    i += 1

    print("final list : " + str(len(list_to_return)))
    return list_to_return
    #return regions_stats_list

# Parse gff3 format
def parse_regions_file(regions_file, contexts, chr_arg):

    regions_to_test = {}



    for line in open(regions_file, "r").readlines():
        sl = line.rstrip().split()

        chr = sl[0]
        name = sl[1]
        start = int(sl[2])
        end = int(sl[3])

        if contexts == "custom":
            contexts_to_test = sl[5].split(",")
        else:
            contexts_to_test = contexts.split(",")

        if contexts_to_test == ["all"]:
            contexts_to_test = ["CG", "CHG", "CHH", "CGxCHG", "CGxCHH", "CHGxCG", "CHGxCHH", "CGxCHGxCHH"]

        if chr == chr_arg:
            if chr not in regions_to_test:
                regions_to_test[chr] = []

            region_to_test_object = {"name": name, "start": start, "end" : end, "contexts_to_test": contexts_to_test}
            regions_to_test[chr].append(region_to_test_object)

    for chr in regions_to_test:
        print(chr + " : " + str(len(regions_to_test[chr])))

    return regions_to_test


def fill_genome_dict(genome_dict, sample_conf_data, chr_arg):

    for sample in sample_conf_data:

        for replicate in sample_conf_data[sample]:

            #print("Processing" + replicate + " ...")


            line_counter = 0
            #print(replicate)
            for line in open(replicate, "r").readlines()[1:]:
                sl = line.split("\t")
                chr = sl[0]

                if chr == chr_arg:
                    pos = int(sl[1])
                    context = sl[3]

                    #ratio = float(sl[4])

                    c_count = int(sl[6])
                    ct_count = int(sl[7])
                    #ratio = binarize_ratio(ratio)

                    if genome_dict[chr][pos-1] == -1:
                        genome_dict[chr][pos-1] = {"context": context, "ratios":{}}

                    if sample not in genome_dict[chr][pos-1]["ratios"]:
                        genome_dict[chr][pos-1]["ratios"][sample] = []

                    #genome_dict[chr][pos-1]["ratios"][sample].append(ratio)
                    genome_dict[chr][pos-1]["ratios"][sample].append([c_count, ct_count])


                    line_counter += 1
                    if line_counter % 1000000 == 0:
                        #print(str(line_counter) + " cytosines processed")
                        #print(genome_dict[chr][pos-1])
                        print(str(line_counter) + " cytosines processed, task memory usage:", "{0:0.1f}".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*0.0000009536743164),"GB        ",end="\r", flush=True)
            print(replicate + " processed")


    return genome_dict


def parse_sample_conf_file(sample_conf_file):

    sample_conf_data = {}
    samples_names = []
    for line in open(sample_conf_file, "r").readlines():

            sl = line.rstrip().split("=")

            sample_conf_data[sl[0]] = sl[1].split(",")
            samples_names.append(sl[0])

    #print(sample_conf_data)
    return (sample_conf_data, samples_names)


def init_genome(genome_fasta, chr_arg):

    genome_dict = {}

    for record in SeqIO.parse(genome_fasta, "fasta"):

        if record.id == chr_arg:
            len_seq = len(str(record.seq))
            list_init = [-1] * (len_seq + 10000000)

            genome_dict[record.id] = list_init

            print(record.id + " intialization done")

    for chr in genome_dict:
        print(chr + " : " + str(len(genome_dict[chr])))

    return genome_dict

