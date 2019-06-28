import sys
import matplotlib.pyplot as plt

time_x, time_y = [], []
length_x, length_y = [], []
length_histo = []
bins = 2000
with open(snakemake.input.sequence_summary) as f:
    f.readline()
    for line in f:
        filename_fastq, filename_fast5, read_id, run_id, channel, mux, start_time, duration, num_events, passes_filtering, \
        template_start, num_events_template, template_duration, sequence_length_template, mean_qscore_template, strand_score_template, \
        median_template, mad_template, pore_type, experiment_id, sample_id = line.split()
        bin = int(float(start_time)) // 10800
        while len(time_y) <= bin:
            time_y.append([])
        time_y[bin].append(float(sequence_length_template) / float(duration))
    for i in time_y:
        if len(i) == 0:
            i.append(0)
    plt.violinplot(time_y, showextrema=False, points=100, widths=[0.9 for i in range(len(time_y))])
    labels = []
    for i in range(len(time_y)):
        labels.append(str(i*3) + '-' + str((i+1)*3))
    plt.xticks([i for i in range(1, len(time_y) + 1)], labels, rotation=45)
    plt.xlabel('Interval (hours)')
    plt.ylabel('Speed (nucleotides/second)')

    plt.tight_layout()
    plt.savefig(snakemake.output.image)
