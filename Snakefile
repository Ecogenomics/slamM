configfile: "config.yaml"

workdir: config["workdir"]



ruleorder: map_reads_ref > copy_reads
ruleorder: filter_illumina_ref > filter_illumina_ref_interleaved > ill_copy_reads > ill_copy_reads_interleaved
ruleorder: fastqc > fastqc_long
ruleorder: polish_isolate_racon_ill > skip_illumina_polish
ruleorder: final_cov_combo > final_cov_long

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



rule flye_assembly:
    input:
        fastq = "data/long_reads.fastq.gz"
    output:
        fasta = "data/flye/assembly.fasta",
        info = "data/flye/assembly_info.txt"
    conda:
        "envs/flye.yaml"
    threads:
        config["max_threads"]
    shell:
        "flye --nano-raw {input.fastq} --meta -o data/flye -t {threads} -g 200m"


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
        maxcov = 1000,
        rounds = 4
    output:
        fasta = "data/assembly.pol.rac.fasta"
    script:
        "scripts/racon_polish.py"

#
# rule polish_metagenome_medaka:
#     input:
#         reads = config["long_reads"],
#         contigs = "data/assembly.pol.rac.fasta"
#     conda:
#         "envs/medaka.yaml"
#     threads:
#         config["max_threads"]
#     params:
#         model = config["guppy_model"]
#     output:
#         fasta = "data/assembly.pol.med.fasta"
#     shell:
#         "medaka_consensus -i {input.reads} -d {input.contigs} -o data/medaka/ -t {threads} -m {params.model} && " \
#         "cp data/medaka/consensus.fasta {output.fasta}"



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


rule polish_meta_pilon:
    input:
        reads = "data/short_reads.fastq.gz",
        fasta = "data/assembly.pol.rac.fasta"
    output:
        fasta = "data/assembly.pol.pil.fasta",
        bam = "data/pilon.sort.bam"
    threads:
        config["max_threads"]
    conda:
        "envs/pilon.yaml"
    shell:
        "minimap2 -ax sr -t {threads} {input.fasta} {input.reads} | samtools view -b | " \
        "samtools sort -o {output.bam} - && samtools index {output.bam} && " \
        "pilon -Xmx64000m --genome {input.fasta} --frags data/pilon.sort.bam --threads {threads} --output data/assembly.pol.pil --fix bases"


rule polish_meta_racon_ill:
    input:
        reads = "data/short_reads.fastq.gz",
        fasta = "data/assembly.pol.pil.fasta"
    output:
        fasta = "data/assembly.pol.fin.fasta",
        paf = "data/final_assembly.paf.gz"
    threads:
        config["max_threads"]
    conda:
        "envs/racon.yaml"
    shell:
        """zcat {input.reads} | awk '{{if (NR % 8 == 1) {{print $1 "/1"}} else if (NR % 8 == 5) {{print $1 "/2"}} 
        else if (NR % 4 == 3){{print "+"}} else {{print $0}} }}' | gzip > data/short_reads_racon.fastq.gz &&
        minimap2 -x sr -t {threads} {input.fasta} data/short_reads_racon.fastq.gz | gzip > data/final_assembly.paf.gz &&
        racon -t {threads} -u data/short_reads_racon.fastq.gz data/final_assembly.paf.gz {input.fasta} > {output.fasta}"""


rule get_high_cov_contigs:
    input:
        info = "data/flye/assembly_info.txt",
        fasta = "data/assembly.pol.fin.fasta",
        paf = "data/final_assembly.paf.gz"
    output:
        fasta = "data/flye_high_cov.fasta"
    params:
        min_cov_long = 20.0,
        min_cov_short = 20.0
    run:
        import gzip
        ill_cov_dict = {}
        with gzip.open(input.paf) as paf:
            for line in paf:
                query, qlen, qstart, qend, strand, ref, rlen, rstart, rend = line.split()[:9]
                ref = ref.decode("utf-8")[:-6]
                if not ref in ill_cov_dict:
                    ill_cov_dict[ref] = 0.0
                ill_cov_dict[ref] += (int(rend) - int(rstart)) / int(rlen)
        count = 0
        with open(input.info) as f:
            f.readline()
            high_cov_set = set()
            for line in f:
                if line.split()[0] in ill_cov_dict:
                    covme = ill_cov_dict[line.split()[0]]
                else:
                    covme = 0
                if float(line.split()[2]) >= params.min_cov_long or not line.split()[0] in ill_cov_dict or ill_cov_dict[line.split()[0]] <= params.min_cov_short:
                    high_cov_set.add(line.split()[0])
        with open(input.fasta) as f, open(output.fasta, 'w') as o:
            write_line = False
            for line in f:
                if line.startswith('>') and line.split()[0][1:-6] in high_cov_set:
                    write_line = True
                elif line.startswith('>'):
                    write_line = False
                if write_line:
                    o.write(line)



# filter illumina reads against the nanopore assembly
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
        "minimap2 -ax sr -t {threads} {input.reference} {input.fastq} |  samtools view -b | "\
        "samtools sort -o {output.bam} - && samtools index {output.bam} && "\
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
        "if LC_ALL=C gzip -l {input.fastq} | awk 'NR==2 {{exit($2!=0)}}'; then touch {output.fasta}; " \
        "else megahit -t {threads} --12 {input.fastq} -o data/mega_assembly && ln data/mega_assembly/final.contigs.fa data/mega_assembly.fasta; fi"


rule metabat_binning_short:
    input:
         fastq = "data/short_reads.filt.fastq.gz",
         fasta = "data/mega_assembly.fasta"
    output:
         metabat_done = "data/metabat_bins/done",
         bam = "data/short_vs_mega.bam"
    conda:
         "envs/metabat2.yaml"
    shell:
         "minimap2 -ax sr -t {threads} {input.fasta} {input.fastq} |  samtools view -b | "\
         "samtools sort -o data/short_vs_mega.bam - && samtools index {output.bam} && "\
         "jgi_summarize_bam_contig_depths --outputDepth data/megahit.cov data/short_vs_mega.bam && " \
         "mkdir -p data/metabat_bins && " \
         "metabat --seed 89 -m 1500 -l -i {input.fasta} -a data/megahit.cov -o data/metabat_bins/binned_contigs && " \
         "touch data/metabat_bins/done"

rule map_long_mega:
    input:
        fastq = "data/long_reads.fastq.gz",
        fasta = "data/mega_assembly.fasta"
    output:
        bam = "data/long_vs_mega.bam"
    threads:
        config["max_threads"]
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -t {threads} -ax map-ont -a {input.fasta} {input.fastq} |  samtools view -b | " \
        "samtools sort -o {output.bam} - && samtools index {output.bam}"


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
         ' printf " --output-fastq-files "; for file in data/binned_reads/*.short.list; do printf "${{file:0:-5}}.2.fastq.gz "; done; printf "\n") && touch {output}'



rule assemble_pools:
    input:
        fastq = "data/binned_reads/done",
        list = "data/list_of_lists.txt"
    threads:
        config["max_threads"]
    output:
        fasta = "data/unicycler_combined.fa"
    conda:
        "envs/final_assembly.yaml"
    threads:
        config["max_threads"]
    script:
        "scripts/assemble_pools.py"


rule final_cov_combo:
    input:
        short_reads = "data/short_reads.fastq.gz",
        long_reads = "data/long_reads.fastq.gz",
        unicyc_fasta = "data/unicycler_combined.fa",
        flye_fasta = "data/assembly.pol.fin.fasta"
    output:
        short_bam = "data/final_short.sort.bam",
        fasta = "data/final_contigs.fasta",
        long_bam = "data/final_long.sort.bam",
        coverage = "data/final.cov"
    conda:
        "envs/metabat2.yaml"
    threads:
        config["max_threads"]
    shell:
        "cat {input.flye_fasta} {input.unicyc_fasta} > {output.fasta} && " \
        "minimap2 -t {threads} -ax map-ont -a {output.fasta} {input.long_reads} |  samtools view -b | " \
        "samtools sort -o {output.long_bam} - && samtools index {output.long_bam} && " \
        "minimap2 -ax sr -t {threads} {output.fasta} {input.short_reads} |  samtools view -b | "\
        "samtools sort -o {output.short_bam} - && samtools index {output.short_bam} && "\
        "jgi_summarize_bam_contig_depths --outputDepth data/final.cov {output.short_bam}"


rule final_cov_long:
    input:
        long_reads = "data/long_reads.fastq.gz",
        fasta = "data/assembly.pol.rac.fasta"
    output:
        fasta = "data/final_contigs.fasta",
        bam = "data/final_long.sort.bam",
        coverage = "data/final.cov"
    conda:
        "envs/metabat2.yaml"
    threads:
        config["max_threads"]
    shell:
        "minimap2 -t {threads} -ax map-ont -a {input.fasta} {input.long_reads} |  samtools view -b | " \
        "samtools sort -o {output.bam} - && samtools index {output.bam} && " \
        "ln {input.fasta} {output.fasta} && " \
        "jgi_summarize_bam_contig_depths --percentIdentity 80 --outputDepth data/final.cov {output.bam}"



rule metabat_binning_2:
    input:
        coverage = "data/final.cov",
        fasta = "data/final_contigs.fasta"
    output:
        metabat_done = "data/metabat_bins_2/done"
    conda:
        "envs/metabat2.yaml"
    shell:
        "mkdir -p data/metabat_bins_2 && " \
        "metabat --seed 89 -i {input.fasta} -a {input.coverage} -o data/metabat_bins_2/binned_contigs && " \
        "touch data/metabat_bins_2/done"

rule checkm:
    input:
        "data/metabat_bins_2/done"
    output:
        "data/checkm.out"
    conda:
        "envs/checkm.yaml"
    params:
        checkm_folder = config["checkm_folder"]
    threads:
        config["max_threads"]
    shell:
        'var="$(which checkm)" && sed -i "s/\/srv\/whitlam\/bio\/db\/checkm_data\/1.0.0/{params.checkm_folder/g" $\{var:0:-11\}/lib/python2.7/site-packages/checkm/DATA_CONFIG && ' \
        'checkm lineage_wf -t {threads} -x fa data/metabat_bins_2 data/checkm > data/checkm.out'


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


rule gtdbtk:
    input:
        "data/metabat_bins_2/done"
    output:
        summary = "data/gtdbtk/gtdbtk.bac120.summary.tsv"
    params:
        gtdbtk_folder = config['gtdbtk_folder']
    conda:
        "envs/gtdbtk.yaml"
    threads:
        config["max_threads"]
    shell:
        "export GTDBTK_DATA_PATH={params.gtdbtk_folder} && " \
        "gtdbtk classify_wf --cpus {threads} --extension fa --genome_dir data/metabat_bins_2 --out_dir data/gtdbtk"



rule create_webpage:
    input:
        checkm_file = "data/checkm.out",
        metabat_bins = "data/metabat_bins_2/done",
        fasta = "data/final_contigs.fasta",
        long_reads_qc_html = "www/nanoplot/longReadsNanoPlot-report.html",
        short_reads_qc_html = "www/short_reads_fastqc.html",
        genes_gff = "data/genes.gff",
        gtdbtk_summary = "data/gtdbtk/gtdbtk.bac120.summary.tsv"
    output:
        "www/index.html"
    threads:
        config["max_threads"]
    script:
        "scripts/create_sdmass_webpage.py"


#############################
# processing barcoded reads #
#############################


rule process_unsorted_reads:
    input:
         html = "QC/read_qc.html"
    output:
         directory("barcoded_reads")
    threads:
         config["max_threads"]
    conda:
         "envs/qcat.yaml"
    params:
         fastq_pass = config["fastq_pass_dir"]
    shell:
         "cat {params.fastq_pass}/* | qcat --trim -b {output}"



rule read_qc:
    input:
         sequence_summary = config["sequence_summary"]
    output:
         "QC/read_qc.html"
    conda:
        "envs/pycoQC.yaml"
    shell:
         "pycoQC -f {input.sequence_summary} -o {output}"


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
        "flye --nano-raw {input.reads} --threads {threads} -o isolate/flye -g {params.genome_size}"


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
        rounds = 4
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
        "medaka_consensus -i {input.reads} -d {input.contigs} -o isolate/medaka/ -t {threads} -m {params.model} && " \
        "cp isolate/medaka/consensus.fasta {output.fasta}"


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
        "minimap2 -ax sr -t {threads} {input.fasta} {input.reads} | samtools view -b | " \
        "samtools sort -o isolate/pilon.sort.bam - && samtools index isolate/pilon.sort.bam && " \
        "pilon -Xmx64000m --genome {input.fasta} --frags isolate/pilon.sort.bam --threads {threads} --output isolate/isolate.pol.pil --fix bases"


rule polish_isolate_racon_ill:
    input:
        reads = "data/short_reads.fastq.gz",
        fasta = "isolate/isolate.pol.pil.fasta"
    output:
        fasta = "isolate/isolate.pol.fin.fasta"
    threads:
        config["max_threads"]
    conda:
        "envs/racon.yaml"
    shell:
        """zcat {input.reads} | awk '{{if (NR % 8 == 1) {{print $1 "/1"}} else if (NR % 8 == 5) {{print $1 "/2"}} 
        else if (NR % 4 == 3){{print "+"}} else {{print $0}} }}' | gzip > data/short_reads_racon.fastq.gz &&
        minimap2 -x sr -t {threads} {input.fasta} data/short_reads_racon.fastq.gz | gzip > isolate/final_assembly.paf.gz &&
        racon -t {threads} -u data/short_reads_racon.fastq.gz isolate/final_assembly.paf.gz {input.fasta} > {output.fasta}"""

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
        "circlator all {input.fasta} {input.reads} isolate/circlator && cp isolate/circlator/06.fixstart.fasta {output.fasta}"