configfile: "config.yaml"

workdir: config["workdir"]



ruleorder: map_reads_ref > copy_reads
ruleorder: filter_illumina_ref > filter_illumina_ref_interleaved > ill_copy_reads > ill_copy_reads_interleaved
ruleorder: process_combination_assembly > process_long_only
ruleorder: metabat_binning_combined > metabat_binning_long
ruleorder: pool_reads_combo > pool_reads_long
ruleorder: get_fastq_pool_combo > get_fastq_pool_long



# Map nanopore reads to reference genome you want to filter
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


# get list of reads that don't map to genome you want to filter
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


# create new read file with filtered reads
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


# if you don't want to filter the reads using a genome just copy them into the folder
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


# use seqtk to get read lengths
rule get_read_lengths:
    input:
        "data/long_reads.fastq.gz"
    output:
        "data/long_read_stats.txt"
    conda:
        "envs/seqtk.yaml"
    shell:
        "seqtk fqchk {input} > {output}"

# from the seqtk file of read lengths, decide on cutoffs to use for the step down metagenome assembly
rule get_read_cutoffs:
    input:
        "data/long_read_stats.txt"
    output:
        "data/cutoffs.txt"
    params:
        "10"
    script:
        "scripts/get_read_cutoffs_2.py"


# peform the step down metagenome assembly
rule step_down_meta_assembly:
    input:
        fastq = "data/long_reads.fastq.gz",
        cutoffs = "data/cutoffs.txt"
    params:
        minimum_coverage = 30,
        long_only = config["short_reads_1"] == "none"
    output:
        "data/combined_assembly_long.fasta"
    conda:
        "envs/sdmass.yaml"
    threads:
         config["max_threads"]
    script:
        "scripts/step_down_meta_assembly.py"

rule polish_sdmass:
    input:
        fastq = "data/long_reads.fastq.gz",
        fasta = "data/combined_assembly_long.fasta"
    output:
        fasta = "data/combined_assembly_long.pol.fasta"
    params:
        rounds = 2,
        maxcov = 300
    conda:
         "envs/racon.yaml"
    threads:
         config["max_threads"]
    script:
         "scripts/racon_polish.py"


# filter illumina reads that don't map to a user selected reference
rule filter_illumina_ref:
    input:
        fastq_1 = config["short_reads_1"],
        fastq_2 = config["short_reads_2"],
        reference_filter = config["reference_filter"]
    output:
        bam = "data/short_unmapped_ref.bam",
        fastq = "data/short_reads.fastq.gz"
    conda:
        "envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        "minimap2 -ax sr -t {threads} {input.reference_filter} {input.fastq_1} {input.fastq_1}  | samtools view -b -f 12 > {output.bam} &&"\
        "samtools bam2fq {output.bam} | gzip > {output.fastq}"

rule filter_illumina_ref_interleaved:
    input:
        fastq_1 = config["short_reads_1"],
        reference_filter = config["reference_filter"]
    output:
        bam = "data/short_unmapped_ref.bam",
        fastq = "data/short_reads.fastq.gz"
    conda:
        "envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        "minimap2 -ax sr -t {threads} {input.reference_filter} {input.fastq_1} | samtools view -b -f 12 > {output.bam} &&"\
        "samtools bam2fq {output.bam} | gzip > {output.fastq}"


# if no reference provided merge the short reads and copy ot working directory
rule ill_copy_reads:
    input:
        fastq_1 = config["short_reads_1"],
        fastq_2 = config["short_reads_2"]
    output:
        "data/short_reads.fastq.gz"
    conda:
        "envs/seqtk.yaml"
    shell:
        "seqtk mergepe {input.fastq_1} {input.fastq_2} | gzip > {output}"

rule ill_copy_reads_interleaved:
    input:
        fastq = config["short_reads_1"]
    output:
        "data/short_reads.fastq.gz"
    run:
        if input.fastq[-3:] == ".gz":
            shell("ln -s {input.fastq} {output}")
        else:
            shell("cat {input.fastq} | gzip > {output}")


# filter illumina reada against the nanopore wtdbg2 assemblies
rule filter_illumina_wtdbg2:
    input:
        fastq = "data/short_reads.fastq.gz",
        reference = "data/combined_assembly_long.pol.fasta"
    output:
        bam = "data/short_vs_wtdg2.sort.bam",
        fastq = "data/short_reads.filt.fastq.gz"
    conda:
        "envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        "minimap2 -ax sr -t {threads} {input.reference} {input.fastq} | gzip > data/short_vs_wtdbg2.sam.gz && "\
        "zcat data/short_vs_wtdbg2.sam.gz | samtools view -b > data/short_vs_wtdbg2.bam &&"\
        "samtools sort -o {output.bam} data/short_vs_wtdbg2.bam &&"\
        "samtools bam2fq -f 12 {output.bam} | gzip > {output.fastq}"




# assemble filtered illumina reads with megahit
rule megahit_assembly:
    input:
        fastq = "data/short_reads.filt.fastq.gz"
    output:
        fasta = "data/mega_assembly.fasta"
    threads:
         config["max_threads"]
    conda:
        "envs/megahit.yaml"
    shell:
        "megahit -t {threads} --12 {input.fastq} -o data/mega_assembly && ln data/mega_assembly/final.contigs.fa data/mega_assembly.fasta"





rule process_combination_assembly:
    input:
        long_assembly = "data/combined_assembly_long.pol.fasta",
        short_assembly = "data/mega_assembly.fasta",
        illumina_reads = "data/short_reads.filt.fastq.gz",
        ill_vs_wtdbg2_bam = "data/short_vs_wtdg2.sort.bam"
    output:
        fasta = "data/merged_assembly.fasta",
        bam = "data/short_vs_merged_assembly.sort.bam"
    threads:
        config["max_threads"]
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -ax sr -t {threads} -a {input.short_assembly} {input.illumina_reads} | samtools view -b > data/mega_assembly.bam &&" \
        "samtools sort -o data/mega_assembly.sort.bam data/mega_assembly.bam && " \
        "pilon -Xmx16000m --genome {input.short_assembly} --frags data/mega_assembly.bam --threads {threads} --output data/short_pilon_1.fa --fix bases &&"\
        "pilon -Xmx16000m --genome {input.long_assembly} --frags data/short_vs_merged_assembly.sort.bam --threads {threads} --output data/long_pilon_1.fa --fix bases &&"\
        "samtools merge {output.bam} {input.ill_vs_wtdbg2_bam} data/mega_assembly.sort.bam && "\
        "cat data/long_pilon_1.fa data/short_pilon_1.fa > {output.fasta}"



rule process_long_only:
    input:
        fasta = "data/combined_assembly_long.pol.fasta",
        reads = "data/long_reads.fastq.gz"
    output:
        bam = "data/long_reads_vs_comb_long.sort.bam"
    conda:
        "envs/minimap2.yaml"
    threads:
        config["max_threads"]
    shell:
        "minimap2 -t {threads} -ax map-ont -a {input.fasta} {input.reads} |  samtools view -b | " \
        "samtools sort -o {output.bam} - && samtools index {output.bam}"



rule metabat_binning_long:
    input:
         bam = "data/long_reads_vs_comb_long.sort.bam",
         fasta = "data/combined_assembly_long.pol.fasta"
    output:
         metabat_done = "data/metabat_bins/done"
    conda:
         "envs/metabat2.yaml"
    shell:
         "jgi_summarize_bam_contig_depths --percentIdentity 75 --outputDepth data/merged_contigs.cov {input.bam} && " \
         "mkdir -p data/metabat_bins && " \
         "metabat --seed 89 -l -i {input.fasta} -a data/merged_contigs.cov -o data/metabat_bins/binned_contigs && " \
         "touch data/metabat_bins/done"



rule metabat_binning_combined:
    input:
         bam = "data/short_vs_merged_assembly.sort.bam",
         fasta = "data/merged_assembly.fasta"
    output:
         metabat_done = "data/metabat_bins/done"
    conda:
         "envs/metabat2.yaml"
    shell:
         "jgi_summarize_bam_contig_depths --outputDepth data/merged_contigs.cov {input.bam} && " \
         "mkdir -p data/metabat_bins && " \
         "metabat --seed 89 -l -i {input.fasta} -a data/merged_contigs.cov -o data/metabat_bins/binned_contigs && " \
         "touch data/metabat_bins/done"


rule pool_reads_long:
    input:
        long_bam = "data/long_reads_vs_comb_long.sort.bam",
        metabat_done = "data/metabat_bins/done"
    output:
        list = "data/list_of_lists.txt"
    conda:
        "envs/pysam.yaml"
    params:
        long_only = True
    script:
        "scripts/pool_reads.py"



rule pool_reads_combo:
    input:
        long_bam = "data/long_reads_vs_comb_long.sort.bam",
        short_bam = "data/short_vs_merged_assembly.sort.bam",
        metabat_done = "data/metabat_bins/done"
    output:
        list = "data/list_of_lists.txt"
    conda:
        "envs/pysam.yaml"
    params:
        long_only = False
    script:
        "scripts/pool_reads.py"


rule get_fastq_pool_combo:
    input:
        list = "data/list_of_lists.txt",
        long_reads = "data/long_reads.fastq.gz",
        short_reads = "data/short_reads.fastq.gz"
    output:
        fastq_list = "data/binned_reads/done"
    conda:
        "envs/seqtk.yaml"
    shell:
        "gawk '{{print $2}}' {input.list} | while read list; do seqtk subseq {input.long_reads} $list | gzip > $list.fastq.gz; done &&"\
        "gawk '{{print $5}}' {input.list} | while read list; do seqtk seq -1 {input.short_reads} | seqtk subseq - $list | gzip > $list.1.fastq.gz; done &&"\
        "gawk '{{print $6}}' {input.list} | while read list; do seqtk seq -2 {input.short_reads} | seqtk subseq - $list | gzip > $list.2.fastq.gz; done &&"\
        "touch data/binned_reads/done"


rule get_fastq_pool_long:
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


rule assemble_pools:
    input:
        fastq = "data/binned_reads/done",
        list = "data/list_of_lists.txt"
    threads:
        config["max_threads"]
    output:
        summary = "data/canu_assembly_summary.txt"
    conda:
        "envs/final_assembly.yaml"
    script:
        "scripts/assemble_pools.py"

rule cat_pools:
    input:
        summary = "data/canu_assembly_summary.txt"
    output:
        fasta = "data/combined_final_assemblies.fasta"
    script:
        "scripts/concat_canu_unicyc.py"

rule racon_polish_final:
    input:
        fastq = "data/long_reads.fastq.gz",
        fasta = "data/combined_final_assemblies.fasta"
    output:
        fasta = "data/combined_final_assemblies.pol.fasta"
    params:
        rounds = 2
    conda:
        "envs/racon.yaml"
    threads:
        config["max_threads"]
    script:
        "scripts/racon_polish.py"


rule final_cov_long:
    input:
        reads = "data/long_reads.fastq.gz",
        fasta = "data/combined_final_assemblies.pol.fasta"
    output:
        "data/final_cov.sort.bam"
    conda:
        "envs/minimap2.yaml"
    threads:
        config["max_threads"]
    shell:
        "minimap2 -t {threads} -ax map-ont -a {input.fasta} {input.reads} |  samtools view -b | " \
        "samtools sort -o {output.bam} - && samtools index {output.bam}"




rule metabat_binning_2:
    input:
        bam = "data/final_cov.sort.bam",
        fasta = "data/combined_final_assemblies.pol.fasta"
    output:
        metabat_done = "data/metabat_bins_2/done"
    conda:
        "envs/metabat2.yaml"
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth data/combined_final.cov data/{input.bam} && " \
        "mkdir -p data/metabat_bins_2 && " \
        "metabat --seed 89 -l -i {input.fasta} -a data/combined_final.cov -o data/metabat_bins_2/binned_contigs && " \
        "touch data/metabat_bins_2/done"

rule checkm:
    input:
        "data/metabat_bins_2/done"
    output:
        "data/checkm.out"
    conda:
        "envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        "checkm lineage_wf -t {threads} -x fa data/metabat_bins_2 data/checkm > testout.txt"



rule create_webpage:
    input:
        checkm_file = "data/checkm.out",
        metabat_bins = "data/metabat_bins_2/done"
    output:
        "www/index.html"
    threads:
        config["max_threads"]
    script:
        "create_sdmass_webpage.py"
