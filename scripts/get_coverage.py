import subprocess
import os
import sys

try:
    os.makedirs('data/binning_bams')
except OSError:
    pass
already_done = set()


setnum = 0
bam_list = []
already_done.add(snakemake.config['long_reads'])
try:
    os.link('data/final_long.sort.bam', 'data/binning_bams/p0.sort.bam')
    os.link('data/final_long.sort.bam.bai', 'data/binning_bams/p0.sort.bam.bai')
except FileExistsError:
    pass
bam_list.append('data/final_long.sort.bam')
set_num = 1
if snakemake.config['short_reads_2'] != 'none':
    already_done.add((snakemake.config['short_reads_1'], snakemake.config['short_reads_2']))
    try:
        os.link('data/final_short.sort.bam', 'data/binning_bams/p1.sort.bam')
        os.link('data/final_short.sort.bam.bai', 'data/binning_bams/p1.sort.bam.bai')
    except FileExistsError:
        pass
    bam_list.append('data/final_short.sort.bam')
    set_num += 1
elif snakemake.params.short_reads_1 != 'none':
    already_done.add(snakemake.config['short_reads_1'])
    try:
        os.link('data/final_short.sort.bam', 'data/binning_bams/p1.sort.bam')
        os.link('data/final_short.sort.bam.bai', 'data/binning_bams/p1.sort.bam.bai')
    except FileExistsError:
        pass
    bam_list.append('data/final_short.sort.bam')
    set_num += 1

if not snakemake.config['profile_read_list'] == 'none':
    with open(snakemake.config['profile_read_list']) as f:
        for line in f:
            if line.split()[0] == 'nanopore':
                r = line.split()[1]
                if not os.path.exists(r):
                    sys.exit('Read file %s does not exits' % r)
                elif r in already_done:
                    pass
                else:
                    already_done.add(r)
                    subprocess.Popen("minimap2 -t %d -ax map-ont -a %s %s |  samtools view -b | samtools sort -o" \
                                     " data/binning_bams/p%d.sort.bam - && samtools index data/binning_bams/p%d.sort.bam" %
                                     (snakemake.threads, snakemake.input.fasta, r, set_num, set_num), shell=True).wait()
                    bam_list.append("data/binning_bams/p%d.sort.bam" % set_num)
                    set_num += 1
            elif line.split()[0] == 'illumina':
                if len(line.split()) == 3:
                    r1 = line.split()[1]
                    r2 = line.split()[2]
                    if not os.path.exists(r1):
                        sys.exit('Read file %s does not exits' % r1)
                    if not os.path.exists(r2):
                        sys.exit('Read file %s does not exits' % r2)
                    if (r1, r2) in already_done:
                        continue
                    already_done.add((r1, r2))
                    r = r1 + ' ' + r2
                else:
                    r = line.split()[1]
                    if not os.path.exists(r):
                        sys.exit('Read file %s does not exits' % r)
                    elif r in already_done:
                        continue
                    already_done.add(r)
                subprocess.Popen("minimap2 -t %d -ax sr -a %s %s |  samtools view -b | samtools sort -o" \
                                 " data/binning_bams/p%d.sort.bam - && samtools index data/binning_bams/p%d.sort.bam" % \
                                 (snakemake.threads, snakemake.input.fasta, r, set_num, set_num), shell=True).wait()
                bam_list.append("data/binning_bams/p%d.sort.bam" % set_num)
                set_num += 1

subprocess.Popen("jgi_summarize_bam_contig_depths --percentIdentity 70 --outputDepth data/metabat.cov " + " ".join(bam_list), shell=True).wait()
try:
    os.makedirs("data/maxbin_cov/")
except OSError:
    pass
with open("data/metabat.cov") as f, open("data/maxbin.cov.list", "w") as o:
    cov_list = []
    contig_list = []
    f.readline()
    for line in f:
        contig = line.split()[0]
        contig_list.append(contig)
        depths = line.split()[3:]
        cov_list.append([])
        for i in range(0, len(depths), 2):
            cov_list[-1].append(depths[i])
    for i in range(len(cov_list[0])):
        with open("data/maxbin_cov/p%d.cov" % i, 'w') as oo:
            for j in range(len(contig_list)):
                oo.write(contig_list[j] + '\t' + cov_list[j][i] + '\n')
        o.write("data/maxbin_cov/p%d.cov\n" % i)
with open("data/binning_bams/done", 'w') as o:
    o.write('done')




