from Bio.SeqIO import FastaIO

biosamples = [line.strip() for line in open("real_samples.txt", "r")]

arg_specific_cdd = dict()
arg_superfamily_cdd = dict()
for part in range(1, 26):
    last_line_fields = list()
    for line in open("../CDD_features/Part{}_hitdata.txt".format(part), "r"):
        fields = line.strip().split('\t')
        if fields[0] == "Query":
            continue
        fields[0] = fields[0].split('>')[1]
        if fields[0] in arg_specific_cdd:
            if (last_line_fields[3] >= fields[3]) & (last_line_fields[4] <= fields[4]):
                arg_specific_cdd[fields[0]] = fields[7]
                if fields[1] == "specific":
                    arg_superfamily_cdd[fields[0]] = fields[10]
                else:
                    arg_superfamily_cdd[fields[0]] = fields[7]
            else:
                arg_specific_cdd[fields[0]] = arg_specific_cdd[fields[0]] + '&' + fields[7]
                if fields[1] == "specific":
                    arg_superfamily_cdd[fields[0]] = arg_superfamily_cdd[fields[0]] + '&' + fields[10]
                else:
                    arg_superfamily_cdd[fields[0]] = arg_superfamily_cdd[fields[0]] + '&' + fields[7]
        else:
            arg_specific_cdd.update({fields[0]:fields[7]})
            if fields[1] == "specific":
                arg_superfamily_cdd.update({fields[0]:fields[10]})
            else:
                arg_superfamily_cdd.update({fields[0]:fields[7]})
        last_line_fields = fields

specific_cdd_amr_ref_count = dict()
superfamily_cdd_amr_ref_count = dict()
arg_amr_ref_count = dict()
with open("../data/database/v2/features.fasta") as features:
    for values in FastaIO.SimpleFastaParser(features):
        arg = values[0].split('|')[-1].upper()
        amr = values[0].split('|')[-2]
        if values[0] in arg_specific_cdd:
            cdd = arg_specific_cdd[values[0]]
            cdd_sup = arg_superfamily_cdd[values[0]]
            if (cdd, amr) in specific_cdd_amr_ref_count:
                specific_cdd_amr_ref_count[(cdd, amr)] += 1
            else:
                specific_cdd_amr_ref_count.update({(cdd, amr):1})
            if (cdd_sup, amr) in superfamily_cdd_amr_ref_count:
                superfamily_cdd_amr_ref_count[(cdd_sup, amr)] += 1
            else:
                superfamily_cdd_amr_ref_count.update({(cdd_sup, amr):1})
        else:
            if (arg, amr) in specific_cdd_amr_ref_count:
                specific_cdd_amr_ref_count[(arg, amr)] += 1
                superfamily_cdd_amr_ref_count[(arg, amr)] += 1
            else:
                specific_cdd_amr_ref_count.update({(arg, amr):1})
                superfamily_cdd_amr_ref_count.update({(arg, amr):1})
        if (arg, amr) in arg_amr_ref_count:
            arg_amr_ref_count[(arg, amr)] += 1
        else:
            arg_amr_ref_count.update({(arg, amr):1})

for prefix in ["spades/deeparg_results", "deeparg_results"]:
    for perc in [30, 50, 80]:
        output_file = "_".join(prefix.split('/')) + "_{}.tsv".format(perc)
        arg_map_diff_freq = dict()
        specific_map_diff_freq = dict()
        superfamily_map_diff_freq = dict()
        for bio in biosamples:
            input = "{}/{}/arg_alignment_identity_{}/AMR_category_difference.tsv".format(
                bio, prefix, perc)
            for line in open(input, "r"):
                fields = line.strip().split('\t')
                if fields[0] == "line#" or len(fields) == 1:
                    continue
                best_cdd = arg_specific_cdd[fields[1]]
                best_cdd_sup = arg_superfamily_cdd[fields[1]]
                best_arg = fields[1].split('|')[-1].upper()
                best_amr = fields[1].split('|')[-2]
                pred_cdd = arg_specific_cdd[fields[7]]
                pred_cdd_sup = arg_superfamily_cdd[fields[7]]
                pred_arg = fields[2]
                pred_amr = fields[6]
                if ((best_cdd_sup, best_cdd, best_arg, best_amr), 
                    (pred_cdd_sup, pred_cdd, pred_arg, pred_amr)) in arg_map_diff_freq:
                    arg_map_diff_freq[
                        ((best_cdd_sup, best_cdd, best_arg, best_amr), 
                         (pred_cdd_sup, pred_cdd, pred_arg, pred_amr))] += 1
                else:
                    arg_map_diff_freq.update({
                        ((best_cdd_sup, best_cdd, best_arg, best_amr), 
                         (pred_cdd_sup, pred_cdd, pred_arg, pred_amr)):1})
                if ((best_cdd_sup, best_cdd, best_amr), 
                    (pred_cdd_sup, pred_cdd, pred_amr)) in specific_map_diff_freq:
                    specific_map_diff_freq[
                        ((best_cdd_sup, best_cdd, best_amr), 
                         (pred_cdd_sup, pred_cdd, pred_amr))] += 1
                else:
                    specific_map_diff_freq.update({
                        ((best_cdd_sup, best_cdd, best_amr), 
                         (pred_cdd_sup, pred_cdd, pred_amr)):1})
                if ((best_cdd_sup, best_amr), 
                    (pred_cdd_sup, pred_amr)) in superfamily_map_diff_freq:
                    superfamily_map_diff_freq[
                        ((best_cdd_sup, best_amr), 
                         (pred_cdd_sup, pred_amr))] += 1
                else:
                    superfamily_map_diff_freq.update({
                        ((best_cdd_sup, best_amr), 
                         (pred_cdd_sup, pred_amr)):1})

        with open(output_file, "w") as output:
            output.write("Superfamily best-hit\tCDD best-hit\tARG best-hit\tAMR best-hit\t")
            output.write("Superfamily-AMR count best-hit\tCDD-AMR count best-hit\tARG-AMR count best-hit\t")
            output.write("Superfamily pred\tCDD pred\tARG pred\tAMR pred\t")
            output.write("Superfamily-AMR count pred\tCDD-AMR count pred\tARG-AMR count pred\t")
            output.write("Superfamily Mapping Frequency\tCDD Mapping Frequency\tARG Mapping Frequency\t")
            output.write("Superfamily Vice-Versa Frequency\tCDD Vice-Versa Frequency\tARG Vice-Versa Frequency\n")

            for key, arg_value in sorted(arg_map_diff_freq.items()):
                
                best_hit_super_count = superfamily_cdd_amr_ref_count[(key[0][0], key[0][3])]
                best_hit_cdd_count = specific_cdd_amr_ref_count[(key[0][1], key[0][3])]
                best_hit_arg_count = arg_amr_ref_count[(key[0][2], key[0][3])]
                pred_super_count = superfamily_cdd_amr_ref_count[(key[1][0], key[1][3])]
                pred_cdd_count = specific_cdd_amr_ref_count[(key[1][1], key[1][3])]
                pred_arg_count = arg_amr_ref_count[(key[1][2], key[1][3])]

                if (key[1], key[0]) in arg_map_diff_freq:
                    arg_vice_versa = arg_map_diff_freq[(key[1], key[0])]
                else:
                    arg_vice_versa = 0

                cdd_key = (
                    (key[0][0], key[0][1], key[0][3]),
                    (key[1][0], key[1][1], key[1][3]))
                cdd_value = specific_map_diff_freq[cdd_key]
                if (cdd_key[1], cdd_key[0]) in specific_map_diff_freq:
                    cdd_vice_versa = specific_map_diff_freq[(cdd_key[1], cdd_key[0])]
                else:
                    cdd_vice_versa = 0

                super_key = (
                    (key[0][0], key[0][3]),
                    (key[1][0], key[1][3]))
                super_value = superfamily_map_diff_freq[super_key]
                if (super_key[1], super_key[0]) in superfamily_map_diff_freq:
                    super_vice_versa = superfamily_map_diff_freq[(super_key[1], super_key[0])]
                else:
                    super_vice_versa = 0

                output.write('\t'.join(key[0]) + '\t')
                output.write("{}\t{}\t{}\t".format(
                    best_hit_super_count, best_hit_cdd_count, best_hit_arg_count))
                output.write('\t'.join(key[1]) + '\t')
                output.write("{}\t{}\t{}\t".format(
                    pred_super_count, pred_cdd_count, pred_arg_count))
                output.write(str(super_value) + '\t' + str(cdd_value) + '\t' + str(arg_value) + '\t')
                output.write(str(super_vice_versa) + '\t' + str(cdd_vice_versa) + '\t' + str(arg_vice_versa) + '\n')