amr_mapping_frequencies = dict()
arg_mapping_frequencies = dict()

for i, line in enumerate(open('AMR_category_difference.tsv')):
    fields = line.strip().split('\t')
    if fields[0] == 'line#' or len(fields) == 1:
        continue
    pred_amr = fields[6]
    best_amr = fields[1].split('|')[-2]

    if (best_amr, pred_amr) in amr_mapping_frequencies:
        amr_mapping_frequencies[(best_amr, pred_amr)] +=1
    else:
        amr_mapping_frequencies.update({(best_amr, pred_amr):1})

with open('AMR_category_difference_frequency.txt', 'w') as amr_output:
    for keys, value in amr_mapping_frequencies.items():
        amr_output.write("({}, {}): {}\n".format(keys[0], keys[1], value))


for i, line in enumerate(open('ARG_category_difference.tsv')):
    fields = line.strip().split('\t')
    if fields[0] == 'line#' or len(fields) == 1:
        continue
    pred_arg = fields[2]
    best_arg = fields[1].split('|')[-1].upper()

    if (best_arg, pred_arg) in arg_mapping_frequencies:
        arg_mapping_frequencies[(best_arg, pred_arg)] +=1
    else:
        arg_mapping_frequencies.update({(best_arg, pred_arg):1})

with open('ARG_category_difference_frequency.txt', 'w') as arg_output:
    for keys, value in arg_mapping_frequencies.items():
        arg_output.write("({}, {}): {}\n".format(keys[0], keys[1], value))