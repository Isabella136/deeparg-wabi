from collections import Counter
mapping_file = './X.mapping.ARG'
alignment_file = './X.align.daa.tsv'

any_label_best_hit = dict()

last_read = ''
for line in open(alignment_file):
    fields = line.strip().split('\t')
    if fields[0] == last_read:
        read = fields[0]
        best = fields[1]
        if any_label_best_hit[read][1] < float(fields[-1]):
            any_label_best_hit[read] = (best, float(fields[-1]))
        continue
    
    read = fields[0]
    best = fields[1]
    
    any_label_best_hit.update({read: (best, float(fields[-1]))})

    last_read = read


amr_misses = 0
arg_misses = 0
total = 0
pred_c = Counter()
best_c = Counter()

with open('AMR_category_difference.tsv', 'w') as amr_output:
    with open('ARG_category_difference.tsv', 'w') as arg_output:
        for i, line in enumerate(open(mapping_file)):
            fields = line.strip().split('\t')
            if fields[0] == '#ARG':
                amr_output.write('line#\tbest hit\t' + line.strip() + '\n')
                arg_output.write('line#\tbest hit\t' + line.strip() + '\n')
                continue
            total += 1

            pred_arg = fields[0]
            pred_amr = fields[4]
            read = fields[3]

            best_arg = any_label_best_hit[read][0].split('|')[-1].upper()
            best_amr = any_label_best_hit[read][0].split('|')[-2]

            pred_c[pred_amr] += 1
            best_c[best_amr] += 1

            if pred_amr != best_amr: 
                amr_output.write(str(i) + '\t' + any_label_best_hit[read][0] + '\t' + line.strip() + '\n')
                amr_misses += 1
            if pred_arg != best_arg: 
                arg_output.write(str(i) + '\t' + any_label_best_hit[read][0] + '\t' + line.strip() + '\n')
                arg_misses += 1
        amr_output.write("misses: {}/{}".format(amr_misses, total) + '\n')
        arg_output.write("misses: {}/{}".format(arg_misses, total))
        amr_output.write("Predicted: {}".format(pred_c) + '\n')
        amr_output.write("Best-hit: {}".format(best_c))