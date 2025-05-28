from collections import Counter
fn = './X.mapping.ARG'

misses = 0
total = 0
c = Counter()
for i, line in enumerate(open(fn)):
    fields = line.strip().split('\t')
    if fields[0] == '#ARG':
        print('line#\t' + line.strip())
        continue
    total += 1

    pred = fields[4]
    best = fields[5]
    best = best.split('|')[3]

    c[best] += 1

    if best != pred: 
        print(str(i) + '\t' + line.strip())
        misses += 1

print("misses: {}/{}".format(misses, total))
print(c)