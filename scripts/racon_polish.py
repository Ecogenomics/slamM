import sys
import os
import subprocess
import random
import shutil


out = "data/racon_polishing"
reads = snakemake.input.fastq
threads = 24
max_cov = snakemake.params.maxcov

if not os.path.exists(reads):
    sys.exit("No read file found")

try:
    os.makedirs(out)
except FileExistsError:
    pass

random.seed(89)
reference = snakemake.input.fasta

for rounds in range(snakemake.params.rounds):
    paf = os.path.join(out, 'alignment.%s.%d.paf' % (snakemake.params.prefix, rounds)
    subprocess.Popen("minimap2 -t %d -x map-ont %s %s > %s" % (snakemake.threads, reference, reads, paf), shell=True).wait()
    cov_dict = {}
    with open(paf) as f:
        for line in f:
            qname, qlen, qstart, qstop, strand, ref, rlen, rstart, rstop = line.split()[:9]
            qlen, qstart, qstop, rlen, rstart, rstop = map(int, [qlen, qstart, qstop, rlen, rstart, rstop])
            if ref in cov_dict:
                cov_dict[ref] += (rstop - rstart) / rlen
            else:
                cov_dict[ref] = (rstop - rstart) / rlen

    high_cov = set()
    low_cov = set()
    for i in cov_dict:
        if cov_dict[i] >= max_cov:
            high_cov.add(i)
        else:
            low_cov.add(i)

    no_cov = set()
    with open(reference) as ref_file, open(os.path.join(out, "filtered.%d.fa" % rounds), 'w') as o:
        for line in ref_file:
            if line.startswith('>'):
                name = line.split()[0][1:]
                if name in low_cov or name in high_cov:
                    o.write(line)
                    getseq = True
                else:
                    no_cov.add(name)
                    getseq = False
            elif getseq:
                o.write(line)

    included_reads = set()
    with open(paf) as f, open(os.path.join(out, "filtered.%d.paf" % rounds), "w") as paf_file:
        for line in f:
            qname, qlen, qstart, qstop, strand, ref, rlen, rstart, rstop = line.split()[:9]
            qlen, qstart, qstop, rlen, rstart, rstop = map(int, [qlen, qstart, qstop, rlen, rstart, rstop])
            if ref in low_cov:
                paf_file.write(line)
                included_reads.add(qname)
            elif ref in high_cov:
                sample_rate = max_cov / cov_dict[ref]
                if random.random() < sample_rate:
                    included_reads.add(qname)
                    paf_file.write(line)


    with open(os.path.join(out, "reads.%d.lst" % rounds), "w") as o:
        for i in included_reads:
            o.write(i + '\n')
    subprocess.Popen("seqtk subseq %s %s/reads.%d.lst | gzip > %s/reads.%d.fastq.gz" % (reads, out, rounds, out, rounds), shell=True).wait()
    subprocess.Popen("racon -m 8 -x -6 -g -8 -w 500 -t %d -u %s/reads.%d.fastq.gz %s/filtered.%d.paf %s/filtered.%d.fa"
                     " > %s/filtered.%d.pol.fa" % (threads, out, rounds, out, rounds, out, rounds, out, rounds), shell=True).wait()

    with open(os.path.join(out, "combined.%d.pol.fa" % rounds), "w") as o:
        with open(os.path.join(out, "filtered.%d.pol.fa" % rounds)) as f:
            gotten_set = set()
            for line in f:
                if line.startswith('>'):
                    gotten_set.add(line.split()[0][1:])
                o.write(line)
        with open(reference) as f:
            for line in f:
                if line.startswith('>'):
                    name = line.split()[0][1:]
                    if name in gotten_set:
                        get_line = False
                    else:
                        get_line = True
                if get_line:
                    o.write(line)
    reference = os.path.join(out, "combined.%d.pol.fa" % rounds)

shutil.copy2(reference, snakemake.output.fasta)