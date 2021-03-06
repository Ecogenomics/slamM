configfile: "config.yaml"

workdir: config["workdir"]

ruleorder: skip_long_assembly > get_reads_list_ref > copy_reads > short_only
ruleorder: filter_illumina_assembly > short_only
ruleorder: filter_illumina_ref > filter_illumina_ref_interleaved > ill_copy_reads > ill_copy_reads_interleaved
ruleorder: fastqc > fastqc_long
ruleorder: polish_isolate_racon_ill > skip_illumina_polish
ruleorder: combine_assemblies > combine_long_only
ruleorder: instrain > instrain_long
ruleorder: skip_long_assembly > get_high_cov_contigs > short_only
ruleorder: skip_long_assembly > filter_illumina_assembly


onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

onstart:
    long_reads = config["long_reads"]
    short_reads_1 = config["short_reads_1"]
    short_reads_2 = config["short_reads_2"]
    unassembled_long = config["unassembled_long"]
    reference_filter = config["reference_filter"]
    meta_genome_size = config["meta_genome_size"]
    gtdbtk_folder = config["gtdbtk_folder"]
    busco_folder = config["busco_folder"]
    profile_read_list = config["profile_read_list"]
    import os
    if long_reads == "none" and short_reads_1 == "none":
        sys.exit("Need at least one of long_reads, short_reads_1")
    if long_reads != "none" and not os.path.exists(long_reads):
        sys.exit("long_reads does not point to a file")
    if short_reads_1 != "none" and not os.path.exists(short_reads_1):
        sys.exit("short_reads_1 does not point to a file")
    if short_reads_2 != "none" and not os.path.exists(short_reads_2):
        sys.exit("short_reads_2 does not point to a file")
    if unassembled_long != "none" and not os.path.exists(unassembled_long):
        sys.exit("unassembled_long does not point to a file")
    if gtdbtk_folder != "none" and not os.path.exists(gtdbtk_folder):
        sys.stderr.write("gtdbtk_folder does not point to a folder\n")
    if busco_folder != "none" and not os.path.exists(busco_folder):
        sys.stderr.write("busco_folder does not point to a folder\n")


# Filter reads against a reference, i.e. for removing host contamination from the metagenome
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
        "minimap2 -ax map-ont -t {threads} {input.reference_filter} {input.fastq} | samtools view -b  > {output}"


# Get a list of reads that don't map to genome you want to filter
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


# Create new read file with filtered reads
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


# If you don't want to filter the reads using a genome just copy them into the folder
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

# Assembly long reads with metaflye
rule flye_assembly:
    input:
        fastq = "data/long_reads.fastq.gz"
    output:
        fasta = "data/flye/assembly.fasta",
        graph = "data/flye/assembly_graph.gfa",
        info = "data/flye/assembly_info.txt"
    params:
        genome_size = config["meta_genome_size"]
    conda:
        "envs/flye.yaml"
    threads:
        config["max_threads"]
    shell:
        "flye --nano-raw {input.fastq} --meta -o data/flye -t {threads} -g {params.genome_size}"

# Polish the long reads assembly with Racon
rule polish_metagenome_racon:
    input:
        fastq = "data/long_reads.fastq.gz",
        fasta = "data/flye/assembly.fasta"
    conda:
        "envs/racon.yaml"
    threads:
        config["max_threads"]
    params:
        prefix = "racon",
        maxcov = 200,
        rounds = 3,
        illumina = False
    output:
        fasta = "data/assembly.pol.rac.fasta"
    script:
        "scripts/racon_polish.py"


### Steps if illumina data exists
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
        """
        minimap2 -ax sr -t {threads} {input.reference_filter} {input.fastq_1} {input.fastq_1}  |
        samtools view -b -f 12 > {output.bam} && \
        samtools bam2fq {output.bam} | gzip > {output.fastq}
        """


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
        """
        minimap2 -ax sr -t {threads} {input.reference_filter} {input.fastq_1} |
        samtools view -b -f 12 > {output.bam} && \
        samtools bam2fq {output.bam} | gzip > {output.fastq}
        """


# If no reference provided merge the short reads and copy to working directory
rule ill_copy_reads:
    input:
        fastq_1 = config["short_reads_1"],
        fastq_2 = config["short_reads_2"]
    output:
        "data/short_reads.fastq.gz"
    conda:
        "envs/seqtk.yaml"
    shell:
        "rename.sh prefix=SLAM in={input.fastq_1} in2={input.fastq_2} out={output} addpairnum=f"


rule ill_copy_reads_interleaved:
    input:
        fastq = config["short_reads_1"]
    output:
        "data/short_reads.fastq.gz"
    conda:
        "envs/seqtk.yaml"
    shell:
        "rename.sh prefix=SLAM in={input.fastq_1} out={output} addpairnum=f int=t"

# The racon polished long read assembly is polished again with the short reads using Pilon
rule polish_meta_pilon:
    input:
        reads = "data/short_reads.fastq.gz",
        fasta = "data/assembly.pol.rac.fasta"
    output:
        fasta = "data/assembly.pol.pil.fasta",
        bam = "data/pilon.sort.bam"
    threads:
        config["max_threads"]
    params:
        pilon_memory = config["pilon_memory"]
    conda:
        "envs/pilon.yaml"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.fasta} {input.reads} | samtools view -b | 
        samtools sort -o {output.bam} - && \
        samtools index {output.bam} && \
        pilon -Xmx{params.pilon_memory}000m --genome {input.fasta} --frags data/pilon.sort.bam \
        --threads {threads} --output data/assembly.pol.pil --fix bases
        """

# The assembly polished with Racon and Pilon is polished again with the short reads using Racon
rule polish_meta_racon_ill:
    input:
        fastq = "data/short_reads.fastq.gz",
        fasta = "data/assembly.pol.pil.fasta"
    output:
        fasta = "data/assembly.pol.fin.fasta",
        paf = "data/racon_polishing/alignment.racon_ill.0.paf"
    threads:
        config["max_threads"]
    conda:
        "envs/racon.yaml"
    params:
        prefix = "racon_ill",
        maxcov = 200,
        rounds = 1,
        illumina = True
    script:
        "scripts/racon_polish.py"

# High coverage contigs are identified
rule get_high_cov_contigs:
    input:
        info = "data/flye/assembly_info.txt",
        fasta = "data/assembly.pol.fin.fasta",
        graph = "data/flye/assembly_graph.gfa",
        paf = "data/racon_polishing/alignment.racon_ill.0.paf"
    output:
        fasta = "data/flye_high_cov.fasta"
    params:
        min_cov_long = 20.0,
        min_cov_short = 10.0,
        short_contig_size = 200000,
        long_contig_size = 500000
    run:
        ill_cov_dict = {}
        with open(input.paf) as paf:
            for line in paf:
                query, qlen, qstart, qend, strand, ref, rlen, rstart, rend = line.split()[:9]
                ref = ref[:-6]
                if not ref in ill_cov_dict:
                    ill_cov_dict[ref] = 0.0
                ill_cov_dict[ref] += (int(rend) - int(rstart)) / int(rlen)
        count = 0
        with open(input.info) as f:
            f.readline()
            high_cov_set = set()
            short_edges = {}
            for line in f:
                if int(line.split()[1]) > params.long_contig_size:
                    high_cov_set.add(line.split()[0])
                elif int(line.split()[1]) < params.short_contig_size:
                    se1 = line.split()[6].split(',')[0]
                    if se1.startswith('-'):
                        se1 = ("edge_" + se1[1:], True)
                    else:
                        se1 = ("edge_" + se1, False)
                    se2 = line.split()[6].split(',')[-1]
                    if se2.startswith('-'):
                        se2 = ("edge_" + se2[1:], False)
                    else:
                        se2 = ("edge_" + se2, True)
                    if not se1 in short_edges:
                        short_edges[se1] = []
                    short_edges[se1].append(line.split()[0])
                    if not se2 in short_edges:
                        short_edges[se2] = []
                    short_edges[se2].append(line.split()[0])
                if float(line.split()[2]) >= params.min_cov_long or not line.split()[0] in ill_cov_dict or ill_cov_dict[line.split()[0]] <= params.min_cov_short:
                    high_cov_set.add(line.split()[0])
        with open(input.graph) as f:
            filtered_contigs = set()
            for line in f:
                if line.startswith("L"):
                    if (line.split()[1], line.split()[2] == '+') in short_edges and (line.split()[3], line.split()[4] == '-') in short_edges and not line.split()[1] == line.split()[3]:
                        for i in short_edges[(line.split()[1], line.split()[2] == '+')]:
                            filtered_contigs.add(i)
                        for i in short_edges[(line.split()[3], line.split()[4] == '-')]:
                            filtered_contigs.add(i)
        for i in filtered_contigs:
            try:
                high_cov_set.remove(i)
            except KeyError:
                pass
        with open(input.fasta) as f, open(output.fasta, 'w') as o:
            write_line = False
            for line in f:
                if line.startswith('>') and line.split()[0][1:-6] in high_cov_set:
                    write_line = True
                elif line.startswith('>'):
                    write_line = False
                if write_line:
                    o.write(line)


# Illumina reads are filtered against the nanopore assembly.
# Specifically, short reads that do not map to the high coverage long contigs are collected
rule filter_illumina_assembly:
    input:
        reference = "data/flye_high_cov.fasta",
        fastq = "data/short_reads.fastq.gz"
    output:
        bam = "data/sr_vs_long.sort.bam",
        fastq = "data/short_reads.filt.fastq.gz"
    conda:
        "envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        """
        minimap2 -ax sr -t {threads} {input.reference} {input.fastq} |  samtools view -b |
        samtools sort -o {output.bam} - && \
        samtools index {output.bam} && \
        samtools bam2fq -f 12 {output.bam} | gzip > {output.fastq}
        """

# If unassembled long reads are provided, skip the long read assembly
rule skip_long_assembly:
    input:
        fastq = "data/short_reads.fastq.gz",
        unassembled_long = config["unassembled_long"]
    output:
        fastq = "data/short_reads.filt.fastq.gz",
        fasta = "data/flye_high_cov.fasta",
        long_reads = "data/long_reads.fastq.gz"
    shell:
        """
        ln {input.fastq} {output.fastq} && \
        touch {output.fasta} && \
        ln {input.unassembled_long} {output.long_reads}
        """

# If only short reads are provided
rule short_only:
    input:
        fastq = "data/short_reads.fastq.gz"
    output:
        fastq = "data/short_reads.filt.fastq.gz",
        fasta = "data/flye_high_cov.fasta",
        long_reads = "data/long_reads.fastq.gz"
    shell:
        """
        ln {input.fastq} {output.fastq} && \
        touch {output.fasta} && \
        touch {output.long_reads}
        """

# Short reads that did not map to the long read assembly are hybrid assembled with metaspades
# If no long reads were provided, long_reads.fastq.gz will be empty
rule spades_assembly:
    input:
        fastq = "data/short_reads.filt.fastq.gz",
        long_reads = "data/long_reads.fastq.gz"
    output:
        fasta = "data/spades_assembly.fasta"
    threads:
         config["max_threads"]
    params:
         max_memory = config["max_memory"]
    conda:
        "envs/spades.yaml"
    shell:
        """
        minimumsize=500000 && \
        actualsize=$(stat -c%s data/short_reads.filt.fastq.gz) && \
        if [ $actualsize -ge $minimumsize ]
        then
            spades.py --memory {params.max_memory} --meta --nanopore {input.long_reads} --12 {input.fastq} \
            -o data/spades_assembly -t {threads} -k 21,33,55,81,99,127 && \
            ln data/spades_assembly/scaffolds.fasta data/spades_assembly.fasta
        else
            touch {output.fasta}
        fi 
        """

# Short reads are mapped to the spades assembly and jgi_summarize_bam_contig_depths from metabat
# used to calculate the mean coverage of the contigs. The coverage and assembly are then used to bin with metabat.
rule metabat_binning_short:
    input:
         fastq = "data/short_reads.filt.fastq.gz",
         fasta = "data/spades_assembly.fasta"
    output:
         metabat_done = "data/metabat_bins/done",
         bam = "data/short_vs_mega.bam"
    conda:
         "envs/metabat2.yaml"
    threads:
         config["max_threads"]
    shell:
         """
         minimap2 -ax sr -t {threads} {input.fasta} {input.fastq} |  samtools view -b | 
         samtools sort -o data/short_vs_mega.bam - && \
         samtools index {output.bam} && \
         jgi_summarize_bam_contig_depths --outputDepth data/megahit.cov data/short_vs_mega.bam && \
         mkdir -p data/metabat_bins && \
         metabat --seed 89 --unbinned -m 1500 -l -i {input.fasta} -a data/megahit.cov \
         -o data/metabat_bins/binned_contigs && \
         touch data/metabat_bins/done
         """

# Long reads are mapped to the spades assembly
rule map_long_mega:
    input:
        fastq = "data/long_reads.fastq.gz",
        fasta = "data/spades_assembly.fasta"
    output:
        bam = "data/long_vs_mega.bam"
    threads:
        config["max_threads"]
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        minimap2 -t {threads} -ax map-ont -a {input.fasta} {input.fastq} |  samtools view -b |
        samtools sort -o {output.bam} - && \
        samtools index {output.bam}
        """

# Long and short reads that mapped to the spades assembly are pooled (binned) together
rule pool_reads:
    input:
        long_bam = "data/long_vs_mega.bam",
        short_bam = "data/short_vs_mega.bam",
        metabat_done = "data/metabat_bins/done"
    output:
        list = "data/list_of_lists.txt"
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/pool_reads.py"

# Binned read lists are processed to extract the reads associated with each bin
rule get_read_pools:
    input:
        long_reads = "data/long_reads.fastq.gz",
        short_reads = "data/short_reads.fastq.gz",
        list = "data/list_of_lists.txt"
    output:
        "data/binned_reads/done"
    conda:
         "envs/mfqe.yaml"
    shell:
         'eval $(printf "zcat {input.long_reads} | mfqe --fastq-read-name-lists "; for file in data/binned_reads/*.long.list; do printf "$file "; done;'\
         ' printf " --output-fastq-files "; for file in data/binned_reads/*.long.list; do printf "${{file:0:-5}}.fastq.gz "; done; printf "\n") && ' \
         'eval $(printf "seqtk seq -1 {input.short_reads} | mfqe --fastq-read-name-lists "; for file in data/binned_reads/*.short.list; do printf "$file "; done;'\
         ' printf " --output-fastq-files "; for file in data/binned_reads/*.short.list; do printf "${{file:0:-5}}.1.fastq.gz "; done; printf "\n") && ' \
         'eval $(printf "seqtk seq -2 {input.short_reads} | mfqe --fastq-read-name-lists "; for file in data/binned_reads/*.short.list; do printf "$file "; done;'\
         ' printf " --output-fastq-files "; for file in data/binned_reads/*.short.list; do printf "${{file:0:-5}}.2.fastq.gz "; done; printf "\n") && touch {output} || touch {output}'

# Short and long reads for each bin are hybrid assembled with Unicycler
rule assemble_pools:
    input:
        fastq = "data/binned_reads/done",
        list = "data/list_of_lists.txt",
        fasta = "data/spades_assembly.fasta",
        metabat_done = "data/metabat_bins/done"
    threads:
        config["max_threads"]
    output:
        fasta = "data/unicycler_combined.fa"
    conda:
        "envs/final_assembly.yaml"
    script:
        "scripts/assemble_pools.py"

# The long read high coverage assembly from flye and hybrid assembly from unicycler are combined.
# Long and short reads are mapped to this combined assembly.
rule combine_assemblies:
    input:
        short_reads = "data/short_reads.fastq.gz",
        long_reads = "data/long_reads.fastq.gz",
        unicyc_fasta = "data/unicycler_combined.fa",
        flye_fasta = "data/flye_high_cov.fasta"
    output:
        short_bam = "data/final_short.sort.bam",
        fasta = "data/final_contigs.fasta",
        long_bam = "data/final_long.sort.bam"
    conda:
        "envs/metabat2.yaml"
    priority: 1
    threads:
        config["max_threads"]
    shell:
        """
        cat {input.flye_fasta} {input.unicyc_fasta} > {output.fasta} && \
        minimap2 -t {threads} -ax map-ont -a {output.fasta} {input.long_reads} |  samtools view -b | 
        samtools sort -o {output.long_bam} - && samtools index {output.long_bam} && \
        minimap2 -ax sr -t {threads} {output.fasta} {input.short_reads} |  samtools view -b |
        samtools sort -o {output.short_bam} - && \
        samtools index {output.short_bam}
        """


rule combine_long_only:
    input:
        long_reads = "data/long_reads.fastq.gz",
        fasta = "data/assembly.pol.rac.fasta"
    output:
        fasta = "data/final_contigs.fasta",
        bam = "data/final_long.sort.bam"
    priority: 1
    conda:
        "envs/metabat2.yaml"
    threads:
        config["max_threads"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont -a {input.fasta} {input.long_reads} |  samtools view -b | 
        samtools sort -o {output.bam} - && \
        samtools index {output.bam} && \
        ln {input.fasta} {output.fasta}
        """

# jgi_summarize_bam_contig_depths from metabat is used to calculate the coverage of each final contig.
# Only reads with at least an identity of 70% are considered
rule prepare_binning_files:
    input:
        long_reads = "data/long_reads.fastq.gz",
        fasta = "data/final_contigs.fasta",
    output:
        maxbin_coverage = "data/maxbin.cov.list",
        bams = "data/binning_bams/done",
        metabat_coverage = "data/metabat.cov"
    conda:
        "envs/metabat2.yaml"
    threads:
        config["max_threads"]
    script:
        "scripts/get_coverage.py"

# Bin contigs with maxbin
rule maxbin_binning:
    input:
        fasta = "data/final_contigs.fasta",
        maxbin_cov = "data/maxbin.cov.list"
    output:
        "data/maxbin2_bins/done"
    conda:
        "envs/maxbin2.yaml"
    shell:
        """
        mkdir -p data/maxbin2_bins && \
        run_MaxBin.pl -contig {input.fasta} -abund_list {input.maxbin_cov} -out data/maxbin2_bins/maxbin && \
        touch data/maxbin2_bins/done
        """

# Bin contigs with CONCOCT
rule concoct_binning:
    input:
        fasta = "data/final_contigs.fasta",
        bam_done = "data/binning_bams/done"
    output:
        "data/concoct_bins/done"
    conda:
        "envs/concoct.yaml"
    threads:
        config["max_threads"]
    shell:
        """
        mkdir -p data/concoct_working && \
        cut_up_fasta.py {input.fasta} -c 10000 -o 0 --merge_last \
        -b data/concoct_working/contigs_10K.bed > data/concoct_working/contigs_10K.fa && \
        concoct_coverage_table.py data/concoct_working/contigs_10K.bed data/binning_bams/*.sort.bam \
        > data/concoct_working/coverage_table.tsv && \
        concoct --threads {threads}  --composition_file data/concoct_working/contigs_10K.fa \
        --coverage_file data/concoct_working/coverage_table.tsv -b data/concoct_working/ && \
        merge_cutup_clustering.py data/concoct_working/clustering_gt1000.csv \
        > data/concoct_working/clustering_merged.csv && \
        mkdir -p data/concoct_bins && \
        extract_fasta_bins.py {input.fasta} data/concoct_working/clustering_merged.csv \
        --output_path data/concoct_bins/ && \
        touch data/concoct_bins/done
        """

# Bin contigs with metabat in various modes
rule metabat_binning_2:
    input:
        coverage = "data/metabat.cov",
        fasta = "data/final_contigs.fasta"
    output:
        metabat_done = "data/metabat_bins_2/done",
        metabat_sspec = "data/metabat_bins_sspec/done",
        metabat_vspec = "data/metabat_bins_ssens/done",
        metabat_sense = "data/metabat_bins_sens/done"
    conda:
        "envs/metabat2.yaml"
    shell:
        """
        mkdir -p data/metabat_bins_2 && \
        metabat --seed 89 -i {input.fasta} -a {input.coverage} -o data/metabat_bins_2/binned_contigs && \
        touch data/metabat_bins_2/done && \
        metabat1 --seed 89 --sensitive -i {input.fasta} -a {input.coverage} \
        -o data/metabat_bins_sens/binned_contigs && \
        touch data/metabat_bins_sens/done && \
        metabat1 --seed 89 --supersensitive -i {input.fasta} -a {input.coverage} \
        -o data/metabat_bins_ssens/binned_contigs && \
        touch data/metabat_bins_ssens/done && \
        metabat1 --seed 89 --superspecific -i {input.fasta} -a {input.coverage} \
        -o data/metabat_bins_sspec/binned_contigs && \
        touch data/metabat_bins_sspec/done
        """

# DASTool is used to select an optimal, non-redundant set of bins
rule das_tool:
    input:
        fasta = "data/final_contigs.fasta",
        metabat2_done = "data/metabat_bins_2/done",
        concoct_done = "data/concoct_bins/done",
        maxbin_done = "data/maxbin2_bins/done",
        metabat_sspec = "data/metabat_bins_sspec/done",
        metabat_ssens = "data/metabat_bins_ssens/done",
        metabat_sense = "data/metabat_bins_sens/done"
    output:
        das_tool_done = "data/das_tool_bins/done"
    conda:
        "envs/das_tool.yaml"
    threads:
        config["max_threads"]
    shell:
        """
        if [[ ! -f data/maxbin_bins.tsv ]]
        then
            Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2 -e fa > data/metabat_bins_2.tsv && \
            Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/metabat_bins_sspec.tsv && \
            Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/metabat_bins_ssens.tsv && \
            Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/metabat_bins_sens.tsv && \
            Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/concoct_bins.tsv && \
            Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/maxbin_bins.tsv
        fi
        DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads} \
        -i data/metabat_bins_2.tsv,data/metabat_bins_sspec.tsv,data/metabat_bins_ssens.tsv,data/metabat_bins_sens.tsv,data/concoct_bins.tsv,data/maxbin_bins.tsv \
        -c {input.fasta} -o data/das_tool_bins/das_tool && \
        touch data/das_tool_bins/done
        """

# CheckM is used to evaluate the quality of the bins
rule checkm:
    input:
        "data/das_tool_bins/done"
    output:
        "data/checkm.out"
    priority: 2
    conda:
        "envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} -x fa data/das_tool_bins/das_tool_DASTool_bins data/checkm > data/checkm.out'


rule fastqc:
    input:
        "data/short_reads.fastq.gz"
    output:
        "www/short_reads_fastqc.html"
    conda:
        "envs/fastqc.yaml"
    threads:
        config["max_threads"]
    shell:
        "fastqc -o www {input}"

rule fastqc_long:
    output:
        "www/short_reads_fastqc.html"
    shell:
        'echo "no short reads" > {output}'


rule nanoplot:
    input:
        "data/long_reads.fastq.gz"
    output:
        "www/nanoplot/longReadsNanoPlot-report.html"
    conda:
        "envs/nanoplot.yaml"
    threads:
        config["max_threads"]
    shell:
        "NanoPlot -o www/nanoplot -p longReads --fastq {input}"

rule prodigal:
    input:
        fasta = "data/final_contigs.fasta"
    output:
        "data/genes.gff"
    conda:
        "envs/prodigal.yaml"
    shell:
        "prodigal -i {input.fasta} -f gff -o {output} -p meta"

# GTDB-Tk is used to assign a taxonomy to each bin
rule gtdbtk:
    input:
        "data/das_tool_bins/done"
    output:
        done = "data/gtdbtk/done"
    priority: 2
    params:
        gtdbtk_folder = config['gtdbtk_folder']
    conda:
        "envs/gtdbtk.yaml"
    threads:
        config["max_threads"]
    shell:
        """
        export GTDBTK_DATA_PATH={params.gtdbtk_folder} && \
        gtdbtk classify_wf --cpus {threads} --extension fa \
        --genome_dir data/das_tool_bins/das_tool_DASTool_bins \
        --out_dir data/gtdbtk && \
        touch data/gtdbtk/done
        """

# Run BUSCO on the bins
rule busco:
    input:
        "data/das_tool_bins/done"
    output:
        done = "data/busco/done"
    params:
        busco_folder = config["busco_folder"]
    conda:
        "envs/busco.yaml"
    threads:
        config["max_threads"]
    shell:
        """
        mkdir -p data/busco && \
        cd data/busco && \
        minimumsize=500000 && \
        for file in ../das_tool_bins/das_tool_DASTool_bins/*.fa;do 
            actualsize=$(wc -c <\"$file\"); 
            if [ $actualsize -ge $minimumsize ]; then 
                if [ ! -d bacteria_odb10.${{file:39:-3}} ]; then
                    busco -q -c {threads} -i $file -o bacteria_odb10.${{file:39:-3}} \
                    -l {params.busco_folder}/bacteria_odb10 -m geno;
                fi
                if [ ! -d eukaryota_odb10.${{file:39:-3}} ]; then
                busco -q -c {threads} -i $file -o eukaryota_odb10.${{file:39:-3}} \
                -l {params.busco_folder}/eukaryota_odb10 -m geno; 
                fi
                if [ ! -d embryophyta_odb10.${{file:39:-3}} ]; then
                busco -q -c {threads} -i $file -o embryophyta_odb10.${{file:39:-3}} \
                -l {params.busco_folder}/embryophyta_odb10 -m geno; 
                fi
                if [ ! -d fungi_odb10.${{file:39:-3}} ]; then
                busco -q -c {threads} -i $file -o fungi_odb10.${{file:39:-3}} \
                -l {params.busco_folder}/fungi_odb10 -m geno; 
                fi 
                # busco -q -c {threads} -i $file -o metazoa_odb10.${{file:39:-3}} \
                -l {params.busco_folder}/metazoa_odb10 -m geno; 
                # busco -q -c {threads} -i $file -o protists_ensembl.${{file:39:-3}} \
                -l {params.busco_folder}/protists_ensembl -m geno; 
            fi
        done && \
        cd ../../ && \
        touch data/busco/done
        """


rule instrain_long:
    input:
        bam = "data/final_long.sort.bam",
        fasta = "data/final_contigs.fasta"
    output:
        profile = "data/instrain/output/instrain_scaffold_info.tsv"
    params:
        instrain_params = config["instrain_params"]
    threads:
        config["max_threads"]
    conda:
        "envs/instrain.yaml"
    shell:
        """
        inStrain profile {params.instrain_params} --processes {threads} -o data/instrain {input.bam} {input.fasta}
        """


rule instrain:
    input:
        bam = "data/final_short.sort.bam",
        fasta = "data/final_contigs.fasta"
    output:
        profile = "data/instrain/output/instrain_scaffold_info.tsv"
    params:
        instrain_params = config["instrain_params"]
    threads:
        config["max_threads"]
    conda:
        "envs/instrain.yaml"
    shell:
        """
        inStrain profile {params.instrain_params} --processes {threads} {input.bam} {input.fasta} -o data/instrain
        """

# Create summary webpage
rule create_webpage:
    input:
        checkm_file = "data/checkm.out",
        metabat_bins = "data/das_tool_bins/done",
        busco_done = "data/busco/done",
        fasta = "data/final_contigs.fasta",
        long_reads_qc_html = "www/nanoplot/longReadsNanoPlot-report.html",
        short_reads_qc_html = "www/short_reads_fastqc.html",
        genes_gff = "data/genes.gff",
        gtdbtk_done = "data/gtdbtk/done"
        # strain_profile = "data/instrain/output/instrain_scaffold_info.tsv"
    output:
        "www/index.html"
    threads:
        config["max_threads"]
    conda:
        "envs/webpage.yaml"
    script:
        "scripts/create_slamm_webpage.py"


#############################
# processing barcoded reads #
#############################

rule process_reads:
    input:
         html = "QC/read_qc.html",
         image = "QC/velocity.png"
    output:
        "barcoded_reads/done"
    threads:
        config["max_threads"]
    params:
        fastq_pass = config["fastq_pass_dir"]
    shell:
        """
        for file in {params.fastq_pass}/*
        do
        cat $file/* | gzip > barcoded_reads/${{file##*/}}.fastq.gz
        done && \
        touch {output}
        """


rule read_qc:
    input:
         sequence_summary = config["sequence_summary"]
    output:
         "QC/read_qc.html"
    conda:
        "envs/pycoQC.yaml"
    shell:
         "pycoQC -f {input.sequence_summary} -o {output}"

rule read_velocity:
    input:
        sequence_summary = config["sequence_summary"]
    output:
        image = "QC/velocity.png"
    conda:
        "envs/matplotlib.yaml"
    script:
        "scripts/read_stats_long.py"


####################
# isolate assembly #
####################


rule assemble_reads_flye:
    input:
        reads = config["long_reads"]
    output:
        contigs = "isolate/flye/assembly.fasta"
    conda:
        "envs/flye.yaml"
    params:
        genome_size = config["genome_size"]
    threads:
        config["max_threads"]
    shell:
        "flye --nano-raw {input.reads} --threads {threads} -o isolate/flye -g {params.genome_size} --asm-coverage 100"


rule polish_isolate_racon:
    input:
        fastq = config["long_reads"],
        fasta = "isolate/flye/assembly.fasta"
    conda:
        "envs/racon.yaml"
    threads:
        config["max_threads"]
    params:
        prefix = "second",
        maxcov = 1000,
        rounds = 4,
        illumina = False
    output:
        fasta = "isolate/isolate.pol.rac.fasta"
    script:
        "scripts/racon_polish.py"


rule polish_isolate_medaka:
    input:
        reads = config["long_reads"],
        contigs = "isolate/isolate.pol.rac.fasta"
    conda:
        "envs/medaka.yaml"
    threads:
        config["max_threads"]
    params:
        model = config["guppy_model"]
    output:
        fasta = "isolate/isolate.pol.med.fasta"
    shell:
        """
        medaka_consensus -i {input.reads} -d {input.contigs} -o isolate/medaka/ -t {threads} -m {params.model} && \
        cp isolate/medaka/consensus.fasta {output.fasta}
        """


rule polish_isolate_pilon:
    input:
        reads = "data/short_reads.fastq.gz",
        fasta = "isolate/isolate.pol.med.fasta"
    output:
        fasta = "isolate/isolate.pol.pil.fasta",
    threads:
        config["max_threads"]
    conda:
        "envs/pilon.yaml"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.fasta} {input.reads} | samtools view -b | 
        samtools sort -o isolate/pilon.sort.bam - && \
        samtools index isolate/pilon.sort.bam && \
        pilon -Xmx64000m --genome {input.fasta} --frags isolate/pilon.sort.bam --threads {threads} \
        --output isolate/isolate.pol.pil --fix bases
        """


rule polish_isolate_racon_ill:
    input:
        fastq = "data/short_reads.fastq.gz",
        fasta = "isolate/isolate.pol.pil.fasta"
    output:
        fasta = "isolate/isolate.pol.fin.fasta"
    threads:
        config["max_threads"]
    conda:
        "envs/racon.yaml"
    params:
        prefix = "racon_ill",
        maxcov = 1000,
        rounds = 1,
        illumina = True
    script:
        "scripts/racon_polish.py"


rule skip_illumina_polish:
    input:
        fasta = "isolate/isolate.pol.med.fasta"
    output:
        fasta = "isolate/isolate.pol.fin.fasta"
    shell:
        "cp {input.fasta} {output.fasta}"


rule circlator:
    input:
        fasta = "isolate/isolate.pol.fin.fasta",
        reads = config["long_reads"]
    output:
        fasta = "isolate/completed_assembly.fasta"
    threads:
        config["max_threads"]
    conda:
        "envs/circlator.yaml"
    shell:
        """
        circlator all {input.fasta} {input.reads} isolate/circlator && \
        cp isolate/circlator/06.fixstart.fasta {output.fasta}
        """

