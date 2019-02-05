import os
import subprocess


bin_bases_ill_dict = {}
cov_cutoff_nano = 60
cov_cutoff_ill = 10

try:
    os.makedirs("data/final_assemblies")
except FileExistsError:
    pass


out_assemblies = []
with open(snakemake.input.list) as f:
    for line in f:
        mb_bin, reads, length, bases_nano = line.split()
        reads += '.fastq.gz'
        length, bases_nano = float(length), float(bases_nano)
        if mb_bin in bin_bases_ill_dict:
            bases_ill = bin_bases_ill_dict[mb_bin]
        else:
            bases_ill = 0
        nano_cov = bases_nano / length
        ill_cov = bases_ill / length
        if nano_cov > cov_cutoff_nano or ill_cov < cov_cutoff_ill:
            out_assemblies.append('data/final_assemblies/%s_canu/meta.contigs.fasta' % mb_bin)
            subprocess.Popen("canu -d data/final_assemblies/%s_canu -p meta maxThreads=%d"\
            " useGrid=false genomeSize=%d -nanopore-raw %s" % (mb_bin, snakemake.threads, length, reads), shell=True).wait()


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
        o.write('\t'.join([i.split('/')[2], str(max(length_list)), str(len(length_list))]) + '\n')

