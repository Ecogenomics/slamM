configfile: "config.yaml"

workdir: config["workdir"]


ruleorder: get_reads_list_ref > copy_reads
ruleorder: filter_illumina_ref > ill_copy_reads
ruleorder: process_combination_assembly > process_long_only

rule map_reads_ref:
    input:
        fastq = config["long_reads"],
        reference_filter= config["reference_filter"]
    output:
        "data/raw_mapped_ref.bam"
    conda:
        "envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        "minimap2 -ax map-ont -t {threads} {input.reference_filter} {input.fastq}  | samtools view -b  > {output}"


rule get_umapped_reads_ref:
    input:
        "data/raw_mapped_ref.bam"
    output:
        "data/unmapped_to_ref.list"
    params:
        "no_full"
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/filter_read_list.py"

rule get_reads_list_ref:
    input:
        fastq = config["long_reads"],
        list = "data/unmapped_to_ref.list"
    output:
        "data/long_reads.fastq.gz"
    conda:
        "envs/seqtk.yaml"
    shell:
        "seqtk subseq {input.fastq} {input.list} | gzip > {output}"



rule copy_reads:
    input:
        fastq = config["long_reads"],
    output:
        "data/long_reads.fastq.gz"
    run:
        if input.fastq[-3:] == ".gz":
            shell("ln -s {input.fastq} {output}")
        else:
            shell("cat {input.fastq} | gzip > {output}")
        
rule get_read_lengths:
    input:
        "data/long_reads.fastq.gz"
    output:
        "data/long_read_stats.txt"
    conda:
        "envs/seqtk.yaml"
    shell:
        "seqtk fqchk {input} > {output}"


rule get_read_cutoffs:
    input:
        "data/long_read_stats.txt"
    output:
        "data/cutoffs.txt"
    params:
        "10"
    script:
        "scripts/get_read_cutoffs.py"



rule step_down_meta_assembly:
    input:
        fastq = config["long_reads"],
        cutoffs = "data/cutoffs.txt"
    params:
        minimum_coverage = 60,
        long_only = "short_reads" in config
    output:
        "data/combined_assembly_long.fasta"
    conda:
        "envs/sdmass.yaml"
    threads:
         config["max_threads"]
    script:
        "scripts/step_down_meta_assembly.py"


rule filter_illumina_ref:
    input:
        fastq = config["short_reads"],
        reference_filter= config["reference_filter"]
    output:
        bam = "data/short_unmapped_ref.bam",
        fastq = "data/short_reads.fastq.gz"
    conda:
        "envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        "minimap2 -ax sr -t {threads} {input.reference_filter} {input.fastq}  | samtools view -b -f 12 > {output.bam} &&"\
        "samtools bam2fq {output.bam} | gzip > {output.fastq}"

rule ill_copy_reads:
    input:
        fastq = config["short_reads"],
    output:
        "data/short_reads.fastq.gz"
    run:
        if input.fastq[-3:] == ".gz":
            shell("ln -s {input.fastq} {output}")
        else:
            shell("cat {input.fastq} | gzip > {output}")

rule filter_illumina_wtdbg2:
    input:
        fastq = "data/short_reads.fastq.gz",
        reference = "data/combined_assembly_long.fasta"
    output:
        bam = "data/short_unmapped_wtdg2.bam",
        fastq = "data/short_reads.filt.fastq.gz"
    conda:
        "envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        "minimap2 -ax sr -t {threads} {input.reference} {input.fastq}  | samtools view -b -f 12 > {output.bam} &&"\
        "samtools bam2fq {output.bam} | gzip > {output.fastq}"



rule megahit_assembly:
    input:
        fastq = "data/short_reads.filt.fastq.gz"
    output:
        directory = directory("data/mega_assembly"),
        fasta = "data/mega_assembly/final.contigs.fa"
    threads:
         config["max_threads"]
    conda:
        "envs/megahit.yaml"
    shell:
        "megahit -t {threads} --12 {input.fastq} -o {output.directory}"





rule process_combination_assembly:
    input:
        long_assembly = "data/combined_assembly_long.fasta",
        short_assembly = "data/mega_assembly.fasta",
        illumina_reads = "data/illumina.filt.fastq.gz",
        ill_vs_wtdbg2_bam = "data/ill_vs_wtdbg2.bam"
    output:
        fasta = "data/merged_assembly.fasta",
        coverage = "data/merged_contigs.sort.bam"
    threads:
        config["max_threads"]
    shell:
        "minimap2 -ax sr -t {threads} -a {input.short_assembly} {input.illumina_reads} | samtools view -b > data/mega_assembly.bam &&" \
        "samtools sort -o data/merged_contigs.sort.bam data/long_combined.bam && samtools index data/merged_contigs.sort.bam"



rule process_long_only:
    input:
        fasta = "data/combined_assembly_long.fasta",
        reads = "data/long_reads.fastq.gz"
    output:
        fasta = "data/merged_assembly.fasta",
        bam = "data/merged_contigs.sort.bam"
    conda:
        "envs/minimap2.yaml"
    threads:
        config["max_threads"]
    shell:
        "minimap2 -t {threads} -x map-ont -a {input.fasta} {input.reads} | samtools view -b > data/long_combined.bam &&" \
        "samtools sort -o {output.bam} data/long_combined.bam && samtools index {output.bam} &&" \
        "ln {input.fasta} {output.fasta}"



checkpoint metabat_binning:
    input:
         bam = "data/merged_contigs.sort.bam",
         fasta = "data/merged_assembly.fasta"
    output:
         metabat_done = "data/metabat_bins/done"
    conda:
         "envs/metabat2.yaml"
    shell:
         "jgi_summarize_bam_contig_depths --percentIdentity 75 --outputDepth data/merged_contigs.cov data/merged_contigs.sort.bam && " \
         "mkdir -p data/metabat_bins && " \
         "metabat --seed 89 -l --unbinned -i {input.fasta} -a data/merged_contigs.cov -o data/metabat_bins/binned_contigs && " \
         "touch data/metabat_bins/done"



rule pool_reads_long:
    input:
        bam = "data/merged_contigs.sort.bam",
        metabat_done = "data/metabat_bins/done"
    output:
        list = "data/list_of_lists.txt"
    conda:
        "envs/pysam.yaml"
    params:
        long_only = "short_reads" in config
    script:
        "scripts/pool_reads_long.py"


rule get_fastq_pool:
    input:
        list = "data/list_of_lists.txt",
        reads = "data/long_reads.fastq.gz"
    output:
        fastq_list = "data/binned_reads/done"
    conda:
        "envs/seqtk.yaml"
    shell:
        "gawk '{{print $2}}' {input.list} | while read list; do seqtk subseq {input.reads} $list | gzip > $list.fastq.gz; done &&"\
        "touch data/binned_reads/done"


rule assemble_pools_nano_only:
    input:
        fastq = "data/binned_reads/done",
        list = "data/list_of_lists.txt"
    params:
        illumina = None
    threads:
        config["max_threads"]
    output:
        summary = "data/canu_asembly_summary.txt"
    conda:
        "envs/canu.yaml"
    script:
        "scripts/assemble_pools.py"




