from Bio.SeqIO import FastaIO

biosamples = [line.strip() for line in open("real_samples.txt", "r")]

arg_amr_ref_count = dict()

with open("../data/database/v2/features.fasta") as features:
    for values in FastaIO.SimpleFastaParser(features):
        arg = values[0].split('|')[-1].upper()
        amr = values[0].split('|')[-2]
        if (arg, amr) in arg_amr_ref_count:
            arg_amr_ref_count[(arg, amr)] += 1
        else:
            arg_amr_ref_count.update({(arg, amr):1})

for prefix in ["spades/deeparg_results", "deeparg_results"]:
    for perc in [30, 50, 80]:
        output_file = "_".join(prefix.split('/')) + "_{}.tsv".format(perc)
        map_diff_freq = dict()
        for bio in biosamples:
            input = "{}/{}/arg_alignment_identity_{}/difference_summary.tsv".format(
                bio, prefix, perc)
            long_seq = prefix == "spades/deeparg_results"
            for line in open(input, "r"):
                fields = line.strip().split('\t')
                if fields[0] == "ARG best-hit" or fields[0] == "ORF ID" or len(fields) == 1:
                    continue
                if long_seq:
                    fields = fields[2:]
                    fields.append(1)
                if ((fields[0], fields[1]), (fields[2], fields[3])) in map_diff_freq:
                    map_diff_freq[
                        ((fields[0], fields[1]), (fields[2], fields[3]))
                        ] += int(fields[4])
                else:
                    map_diff_freq.update({
                        ((fields[0], fields[1]), (fields[2], fields[3]))
                        : int(fields[4])})
        with open(output_file, "w") as output:
            output.write("ARG best-hit\tAMR best-hit\tDB count best-hit\tARG pred\tAMR pred\tDB count pred\tFrequency\tVice-Versa Frequency\n")
            for key, value in map_diff_freq.items():
                best_hit_ref_count = arg_amr_ref_count[key[0]]
                pred_ref_count = arg_amr_ref_count[key[1]]
                if (key[1], key[0]) in map_diff_freq:
                    vice_versa = map_diff_freq[(key[1], key[0])]
                else:
                    vice_versa = 0
                output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    key[0][0], key[0][1], best_hit_ref_count, key[1][0], 
                    key[1][1], pred_ref_count, value, vice_versa))