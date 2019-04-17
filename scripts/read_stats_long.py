import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

time_x, time_y = [], []
length_x, length_y = [], []
length_histo = []
bins = 2000
with open(sys.argv[1]) as f:
    f.readline()
    for line in f:
        filename_fastq, filename_fast5, read_id, run_id, channel, mux, start_time, duration, num_events, passes_filtering, \
        template_start, num_events_template, template_duration, sequence_length_template, mean_qscore_template, strand_score_template, \
        median_template, mad_template = line.split()
        bin = int(float(start_time)) // 10800
        while len(time_y) <= bin:
            time_y.append([])
        time_y[bin].append(float(sequence_length_template) / float(duration))
    # columns = []
    # for num, i in enumerate(time_y):
    #     columns.append(str(num * 200))
    # the_data = pd.DataFrame(np.array(time_y[:10]))
    #
    # ax = sns.violinplot(data = the_data, palette = "muted")
    for i in time_y:
        if len(i) == 0:
            i.append(0)
    plt.violinplot(time_y, showextrema=False, points=100, widths=[0.9 for i in range(len(time_y))])
    #plt.set_xticklables('q' for i in range(len(time_y)))
    labels = []
    for i in range(len(time_y)):
        labels.append(str(i*3) + '-' + str((i+1)*3))
    plt.xticks([i for i in range(1, len(time_y) + 1)], labels, rotation=45)
    plt.xlabel('Interval (hours)')
    plt.ylabel('Speed (nucleotides/second)')

    plt.tight_layout()
    plt.savefig(sys.argv[2])

    # heatmap, xedges, yedges = np.histogram2d(time_x, time_y, bins=bins)
    # the_extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    # plt.clf()
    # plt.imshow(heatmap.T, aspect=the_extent[1] / the_extent[3], extent=the_extent)
    # cb = plt.colorbar()
    # cb.set_label('Frequency')
    # plt.title('Nanopore fastq')
    # plt.xlabel('Read length')
    # plt.ylabel('Quality')
    # plt.savefig(sys.argv[2])

