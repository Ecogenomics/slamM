import pysam
import os

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
            if (clipped_start/length <= cutoff and clipped_end/length <= cutoff) or \
                (clipped_start/length <= cutoff and read.reference_end > read.reference_length - 100) or \
                (clipped_end/length <= cutoff and read.reference_start < 100):
                outreads[bin].add(read.query_name)
                outbases[bin] += length

try:
    os.makedirs("data/binned_reads")
except FileExistsError:
    pass

with open(snakemake.output.list, 'w') as o:
    for i in outreads:
        with open("data/binned_reads/r" + i + '.list', 'w') as read_list:
            for j in outreads[i]:
                read_list.write(j + '\n')
        o.write('\t'.join([i, "data/binned_reads/r" + i + '.list', str(outlength[i]), str(outbases[i])]) + '\n')
