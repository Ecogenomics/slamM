import os, gzip


the_dir = snakemake.input[:-5]
with open('barcoded_reads/summary.txt', 'w') as summary:
    summary.write('barcode\tread_no\tread_length\ttotal_length\n')
    for i in os.listdir(the_dir):
        if i.endswith('.fastq'):
            count = 0
            tot_length = 0
            read_file = os.path.join(the_dir, i)
            with gzip.open(read_file) as f, gzip.open(read_file[:-5] + 'clean.fastq.gz', 'wt') as o:
                while l1 != '':
                    l1 = f.readline()
                    l2 = f.readline()
                    l3 = f.readline()
                    l4 = f.readline()
                    if l2.rstrip() != '':
                        o.write(l1 + l2 + l3 + l4)
                        count += 1
                        tot_length += len(l2.rstrip())
            os.remove(read_file)
            summary.write(i + '\t' + str(count) + '\t' + str(tot_length/count) + '\t' + str(tot_length) + '\n')