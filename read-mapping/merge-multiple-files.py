import sys

#Makes an abundance matrix for multiple files - rows=contigs, columns=samples

fp = open(sys.argv[1], 'w')

d = {}
files = sys.argv[2:]
l_data = []

for file in files:
    l_file = []
    fi = open(file, 'r')
    for line in fi:
        l_file.append(line) 
    l_data.append(l_file)

for x in files:
    fp.write('\t%s' % x)

fp.write('\n')

d_transformed = {}

for n, x in enumerate(l_data):
    for m, y in enumerate(x):
        y = y.split(' ')
        contig = y[0]
        abundance = float(y[1])
        if d_transformed.has_key(contig):
            d_transformed[contig].append(abundance)
        else:
            d_transformed[contig] = [abundance]

for key in d_transformed:
    fp.write('%s\t' % key)
    for x in d_transformed[key]:
        fp.write('%f\t' % x)
    fp.write('\n')

