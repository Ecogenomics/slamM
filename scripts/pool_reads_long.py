import pysam
import os
import random

samfile = pysam.AlignmentFile(snakemake.input.bam, 'rb')

contig_bins = {}
outlength = {}
for bins in os.listdir(snakemake.input.metabat_done[:-4]):
    if not bins.startswith('binned_contigs') or bins == "binned_contigs.unbinned":
        continue
    bin = bins.split('.')[1]
    outlength[bin] = 0
    with open(os.path.join(snakemake.input.metabat_done[:-4], bins)) as f:
        for line in f:
            contig_bins[line.rstrip()] = bin
            outlength[bin] += samfile.get_reference_length(line.rstrip())



cutoff = 0.05
outreads = {}
outbases = {}
outreads500 = {}
outbases500 = {}
for read in samfile.fetch(until_eof=True):
    start = True
    if not read.cigartuples is None:
        clipped_start = 0
        clipped_end = 0
        for i in read.cigartuples:
            if i[0] == 4 or i[0] == 5:
                if start:
                    clipped_start += i[1]
                else:
                    clipped_end += i[1]
            else:
                start = False
        if read.reference_name in contig_bins:
            bin = contig_bins[read.reference_name]
            length = read.infer_read_length()
            if not bin in outreads:
                outreads[bin] = set()
                outbases[bin] = 0
                outreads500[bin] = set()
                outbases500[bin] = 0
            if (clipped_start/length <= cutoff and clipped_end/length <= cutoff) or \
                (clipped_start/length <= cutoff and read.reference_end > read.reference_length - 100) or \
                (clipped_end/length <= cutoff and read.reference_start < 100):
                outreads[bin].add(read.query_name)
                outbases[bin] += length
                if length >= 500:
                    outreads500[bin].add(read.query_name)
                    outbases500[bin] += length

try:
    os.makedirs("data/binned_reads")
except FileExistsError:
    pass

min_cov_500 = 200
top_cov = 300


with open(snakemake.output.list, 'w') as o:
    for i in outreads:
        with open("data/binned_reads/r" + i + '.list', 'w') as read_list:
            coverage500 =  outbases500[i] / outlength[i]
            if coverage500 >= top_cov:
                fraction_reads_to_get = top_cov/ coverage500
                reads = random.sample(outreads500[i], int(len(outreads500[i]) * fraction_reads_to_get))
                bases = int(fraction_reads_to_get * outbases500[i])
            elif coverage500 >= min_cov_500:
                reads = outreads500[i]
                bases = outbases500[i]
            else:
                reads = outreads[i]
                bases = outbases[i]
            for j in outreads[i]:
                read_list.write(j + '\n')
        o.write('\t'.join([i, "data/binned_reads/r" + i + '.list', str(outlength[i]), str(bases)]) + '\n')
