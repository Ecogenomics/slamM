import os, gzip


the_dir = snakemake.input.done[:-5]
for i in os.listdir(the_dir):
    if i.endswith('.fastq.gz'):
        read_file = os.path.join(the_dir, i)
        with gzip.open(read_file, 'rb') as f, gzip.open(read_file[:-8] + 'clean.fastq.gz', 'wb') as o:
            while l1 != b'':
                l1 = f.readline()
                l2 = f.readline()
                l3 = f.readline()
                l4 = f.readline()
                if l2 != b'\n':
                    o.write(l1 + l2 + l3 + l4)
