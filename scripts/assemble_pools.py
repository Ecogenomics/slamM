import os
import subprocess
import sys


try:
    os.makedirs("data/final_assemblies")
except FileExistsError:
    pass


out_assemblies = []
with open(snakemake.input.list) as f:
    for line in f:
        mb_bin, long_reads, length, bases_nano, short_reads, bases_ill = line.split()
        long_reads = long_reads[:-5] + '.fastq.gz'
        long_reads = os.path.abspath(long_reads)
        short_reads_1 = short_reads[:-5] + '.1.fastq.gz'
        short_reads_2 = short_reads[:-5] + '.2.fastq.gz'
        short_reads_1 = os.path.abspath(short_reads_1)
        short_reads_2 = os.path.abspath(short_reads_2)
        length, bases_nano, bases_ill = float(length), float(bases_nano), float(bases_ill)
        nano_cov = bases_nano / length
        ill_cov = bases_ill / length
        out_assemblies.append('data/final_assemblies/%s_unicyc/assembly.fasta' % mb_bin)
        if not os.path.exists('data/final_assemblies/%s_unicyc/assembly.fasta' % mb_bin):
            subprocess.Popen("unicycler -t %d -1 %s -2 %s -l %s -o data/final_assemblies/%s_unicyc" % (snakemake.threads, short_reads_1, short_reads_2, long_reads, mb_bin), shell=True).wait()



with open(snakemake.fasta, 'w') as o:
    o.write('assembly\tmax_contig\tcontigs\n')
    count = 0
    for i in out_assemblies:
        if not os.path.exists(i):
            with open(i[:-14] + 'unicycler.log') as f:
                lastline = f.readlines()[-1]
                if lastline.startswith("Error: SPAdes failed to produce assemblies. See spades_assembly/assembly/spades.log for more info."):
                    continue
        with open(i) as assembly:
            for line in assembly:
                if line.startswith('>'):
                    count += 1
                    o.write('>unicycler_' + str(count) + '\n')
                else:
                    o.write(line)