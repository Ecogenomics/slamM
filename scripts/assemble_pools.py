import os
import subprocess


bin_bases_ill_dict = {}
cov_cutoff_nano = 60
cov_cutoff_ill = 20

try:
    os.makedirs("data/final_assemblies")
except FileExistsError:
    pass


out_assemblies = []
with open(snakemake.input.list) as f:
    for line in f:
        mb_bin, long_reads, length, bases_nano, short_reads, bases_ill = line.split()
        long_reads += '.fastq.gz'
        short_reads_1 = short_reads + '.short.1.fastq.gz'
        short_reads_2 = short_reads + '.short.2.fastq.gz'
        length, bases_nano, bases_ill = float(length), float(bases_nano), float(bases_ill)
        nano_cov = bases_nano / length
        ill_cov = bases_ill / length
        if nano_cov > cov_cutoff_nano or ill_cov < cov_cutoff_ill:
            out_assemblies.append('data/final_assemblies/%s_canu/meta.contigs.fasta' % mb_bin)
            subprocess.Popen("canu -d data/final_assemblies/%s_canu -p meta maxThreads=%d"\
            " useGrid=false genomeSize=%d -nanopore-raw %s" % (mb_bin, snakemake.threads, length, long_reads), shell=True).wait()
        else:
            out_assemblies.append('data/final_assemblies/%s_unicyc/contigs.fasta' % mb_bin)
            subprocess.Popen("unicycler -t %d -1 %s -2 %s -l %s -o data/final_assemblies/%s_unicyc" % (snakemake.threads, short_reads_1, short_reads_2, long_reads, mb_bin), shell=True).wait()



with open(snakemake.output.summary, 'w') as o:
    o.write('assembly\tmax_contig\tcontigs\n')
    for i in out_assemblies:
        with open(i) as assembly:
            length_list = []
            for line in assembly:
                if line.startswith('>'):
                    length_list.append(0)
                else:
                    length_list[-1] += len(line.rstrip())
        o.write('\t'.join([i.split('/')[2], str(max(length_list)), str(len(length_list)), str(sum(length_list)), i]) + '\n')