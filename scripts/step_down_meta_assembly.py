import subprocess
import os
import pysam

def contig_coverage(bamfile, outfile):
    with open(outfile, 'w') as o:
        samfile = pysam.AlignmentFile(bamfile, 'rb')
        for i, length in zip(samfile.references, samfile.lengths):
            bases = 0
            for j in samfile.pileup(i):
                bases += j.nsegments
            o.write(i + '\t' + str(bases/length) + '\n')


def create_read_list(bamfile, outfile, contigs_to_filter):
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    cutoff = 0.05
    mapped_full = set()
    all_reads = set()
    for read in samfile.fetch():
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
            length = read.infer_read_length()
            if clipped_start/length <= cutoff and clipped_end/length <= cutoff and read.reference_name in contigs_to_filter:
                mapped_full.add(read.query_name)
            all_reads.add(read.query_name)
    out_set = all_reads - mapped_full
    with open(outfile, 'w') as o:
        for i in out_set:
            o.write(i + '\n')

overide = False

with open(snakemake.input.cutoffs) as f, open(snakemake.output[0], 'w') as o:
    try:
        os.makedirs(os.path.join('data', 'wtdbg2assembly'))
    except FileExistsError:
        pass
    curr_reads = snakemake.input.fastq
    contig_count = 1
    for line in f:
        read_length = line.rstrip()
        prefix = os.path.join('data', 'wtdbg2assembly', 'w.' + read_length)

        if not os.path.exists(prefix + '.ctg.lay.gz') or overide:
            if read_length == '500':
                process = subprocess.Popen('wtdbg2 -t %s -p 0 -k 15 -AS 2 -s 0.05 --edge-min 2 --rescue-low-cov- edges -i %s -fo %s -L %s' % (snakemake.threads, curr_reads, prefix, read_length),
                                           shell=True, stderr=subprocess.PIPE)
            else:
                process = subprocess.Popen('wtdbg2 -t %s -i %s -fo %s -L %s' % (snakemake.threads, curr_reads, prefix, read_length),
                                           shell=True, stderr=subprocess.PIPE)
            output = process.communicate()[0]
            if process.returncode != 0:
                try:
                    os.remove(prefix + '.ctg.lay.gz')
                except FileNotFoundError:
                    pass
                raise subprocess.CalledProcessError(process.returncode, 'wtdbg2', output=output)

        if not os.path.exists(prefix + '.fna') or overide:
            process = subprocess.Popen('wtdbg-cns -t %s -i %s.ctg.lay.gz -fo %s.fna' % (snakemake.threads, prefix, prefix),
                                       shell=True, stderr=subprocess.PIPE)
            output = process.communicate()[0]
            if process.returncode != 0:
                try:
                    os.remove(prefix + '.fna')
                except FileNotFoundError:
                    pass
                raise subprocess.CalledProcessError(process.returncode, 'wtdbg-cns', output=output)

        if not os.path.exists(prefix + '.bam') or overide:
            process = subprocess.Popen('minimap2 -ax map-ont -L -t %s %s.fna %s  | samtools view -b  > %s.bam' % (snakemake.threads, prefix, curr_reads, prefix),
                                       shell=True, stderr=subprocess.PIPE)
            output = process.communicate()[0]
            if process.returncode != 0:
                try:
                    os.remove(prefix+ '.bam')
                except FileNotFoundError:
                    pass
                raise subprocess.CalledProcessError(process.returncode, 'minimap2', output=output)

        if not os.path.exists(prefix + '.sort.bam.bai') or overide:
            process = subprocess.Popen('samtools sort -@ %s -o %s.sort.bam %s.bam && samtools index %s.sort.bam' % (snakemake.threads, prefix, prefix, prefix),
                                       shell=True, stderr=subprocess.PIPE)
            output = process.communicate()[0]
            if process.returncode != 0:
                try:
                    os.remove(prefix + '.sort.bam.bai')
                except FileNotFoundError:
                    pass
                raise subprocess.CalledProcessError(process.returncode, 'samtools', output=output)

        if not os.path.exists(prefix + '.coverage'):
            contig_coverage(prefix + '.sort.bam', prefix + '.coverage')

        high_cov_contigs = set()
        with open(prefix + '.coverage') as cov_file:
            for line in cov_file:
                contig, cov = line.split()
                cov = float(cov)
                if cov >= float(snakemake.params.minimum_coverage) or (read_length == '500' and snakemake.params.long_only):
                    high_cov_contigs.add(contig)
        if not os.path.exists(prefix + '.reads.list') or overide:
            create_read_list(prefix + '.sort.bam', prefix + '.reads.list', high_cov_contigs)

        if not os.path.exists(prefix + '.fastq.gz') or overide:
            process = subprocess.Popen('seqtk subseq %s %s.reads.list | gzip > %s.fastq.gz' % (curr_reads, prefix, prefix), shell=True, stderr=subprocess.PIPE)
            output = process.communicate()[0]
            if process.returncode != 0:
                try:
                    os.remove(prefix + '.fastq.gz')
                except FileNotFoundError:
                    pass
                raise subprocess.CalledProcessError(process.returncode, 'seqtk', output=output)

        curr_reads = prefix + '.fastq.gz'

        if not os.path.exists(prefix + '.polished.fna') or overide:
            process = subprocess.Popen('wtpoa-cns -t %s -d %s.fna -i %s.sort.bam -fo %s.polished.fna' % (snakemake.threads, prefix, prefix, prefix),
                                       shell=True, stderr=subprocess.PIPE)
            output = process.communicate()[0]
            if process.returncode != 0:
                try:
                    os.remove(prefix + '.polished.fna')
                except FileNotFoundError:
                    pass
                raise subprocess.CalledProcessError(process.returncode, 'wtpoa-cns', output=output)
        with open(prefix + '.polished.fna') as fasta:
            for line in fasta:
                if line.startswith('>') and line.split()[0][1:] in high_cov_contigs:
                    o.write('>ctg.l.' + str(contig_count) + '\n')
                    contig_count += 1
                    getseq = True
                elif line.startswith('>'):
                    getseq = False
                elif getseq:
                    o.write(line)