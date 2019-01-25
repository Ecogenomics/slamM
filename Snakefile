

rule all:
   output:
       "web/index.html"


configfile: "config.yaml"

workdir: config["workdir"]


ruleorder: get_reads_list_ref > copy_reads

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
            shell("ln {input.fastq} {output}")
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




rule process_combination_assembly:
    input:
        long_assembly = "data/combined_assembly_long.fasta",
        short_assembly = "data/mega_assembly.fasta",
        illumina_reads = "data/illumina.filt.fastq.gz",
        ill_vs_wtdbg2_bam = "data/ill_vs_wtdbg2.bam"
    output:
        fasta = "data/merged_assembly.fasta",
        coverage = "data/merged_contigs.sort.bam"
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



rule metabat_binning:
    input:
         bam = "data/merged_contigs.sort.bam",
         fasta = "data/merged_assembly.fasta"
    output:
         binned_contigs = "data/metabat_bins.unbinned.fa"
    conda:
         "envs/metabat2.yaml"
    shell:
         "jgi_summarize_bam_contig_depths --outputDepth data/merged_contigs.cov &&" \
         "metabat --seed 89 -l --unbinned -i {input.fasta} -a data/merged_contigs.cov -o data/metabat_bins"





rule pool_reads_long:
    input:
        bam = "data/merged_contigs.sort.bam",
        metabat = "data/metabat_bins.{group}"
    output:
        list = "data/binned_reads.{group}.{genome_size}.{assembler}.list"
    wildcard_constraints:
        group = "\d+"
    conda:
        "envs/pysam.yaml"
    params:
        long_only = "short_reads" in config
    script:
        "scripts/pool_reads_long.py"


rule get_fastq_pool:
    input:
        list = "data/binned_reads.{group}.{genome_size}.{assembler}.list",
        reads = "data/long_reads.fastq.gz"
    output:
        fastq = "data/binned_reads.{group}.{genome_size}.{assembler}fastq"
    conda:
        "envs/seqtk.yaml"
    shell:
        "seqtk subset {input.reads} {input.list} | gzip > {output.fastq}"






rule assemble_reads:
    input:
        fastq = "data/binned_reads.{group}.{genome_size}.c.fastq"
    output:
        directory = "data/canu.{group}/",
        fasta = "data/canu.{group}/meta.contigs.fasta"
    conda:
        "envs/canu.yaml"
    shell:
        "canu -d {output.directory} -p meta gnuplotTested=true useGrid=false genomeSize={wildcards.genome_size} -nanopore_raw {input.fastq}"



#
#
# rule final_assembly: