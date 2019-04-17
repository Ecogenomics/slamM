import os
import sys
with open(snakemake.input.barcode_summary) as f:
    f.readline()
    read_dict = {}
    for line in f:
        read_id, barcode_arrangement, barcode_full_arrangement, barcode_kit, barcode_variant, barcode_score, barcode_front_id, barcode_front_score, barcode_front_refseq, barcode_front_foundseq,\
        barcode_front_foundseq_length, barcode_front_begin_index, barcode_rear_id, barcode_rear_score, barcode_rear_refseq, barcode_rear_foundseq, barcode_rear_foundseq_length, barcode_rear_end_index = line.split()
        if not barcode_arrangement == 'unclassified':
            read_dict[read_id] = (barcode_arrangement, int(barcode_front_foundseq_length) + int(barcode_front_begin_index), - int(barcode_rear_foundseq_length) - int(barcode_rear_end_index))
          
try:
    os.makedirs("processed_reads")
except FileExistsError:
    pass

len_dict, count_dict = {},{}
barcode_files = {}
for i in os.listdir(os.path.join(snakemake.input.nano_dir, "fastq_pass")):
    if i.endswith('.fastq'):
        with open(os.path.join(snakemake.input.reads, i)) as f:
            while True:
                name_line = f.readline()
                if name_line == '':
                    break
                seq = f.readline().rstrip()
                plus_line = f.readline()
                qual = f.readline().rstrip()
                name = name_line.split()[0][1:]
                if name in read_dict:
                    barcode_arrangement, front_index, rear_index = read_dict[name]
                    if barcode_arrangement in barcode_files:
                        o = barcode_files[barcode_arrangement]
                    else:
                        o = open(os.path.join(sys.argv[3], barcode_arrangement + '.fastq'), 'w')
                        barcode_files[barcode_arrangement] = o
                    o.write(name_line)
                    o.write(seq[front_index:rear_index] + '\n')
                    if not barcode_arrangement in len_dict:
                        len_dict[barcode_arrangement] = 0
                        count_dict[barcode_arrangement] = 0
                    len_dict[barcode_arrangement] += len(seq[front_index:rear_index])
                    count_dict[barcode_arrangement] += 1
                    o.write(plus_line)
                    o.write(qual[front_index:rear_index] + '\n')


with open(snakemake.output, 'w') as o:
    for i in sorted(list(len_dict)):
        o.write('\t'.join(map(str, [i, len_dict[i], count_dict[i]])) + '\n')


