with open(snakemake.input.summary) as f, open(snakemake.output.fasta, 'w') as o:
    f.readline()
    for line in f:
        name = line.split()[0]
        assembly_file = line.split()[4]
        ctg_no = 1
        with open(assembly_file) as fasta:
            for line in fasta:
                if line.startswith(">"):
                    o.write(">%s.%d\n" % (name, ctg_no))
                    ctg_no += 1
                else:
                    o.write(line)