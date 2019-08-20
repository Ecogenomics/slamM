import sys
import matplotlib.pyplot as plt


ys = []
for the_file in sys.argv[1:]:
    len_list = []
    with open(the_file) as f:
        for line in f:
            if line.startswith('>'):
                len_list.append(0)
            else:
                len_list[-1] += len(line.rstrip())
    len_list.sort(reverse=True)
    y = []
    the_sum = 0
    for i in len_list:
        the_sum += i
        y.append(the_sum)
    ys.append(y)

for y, fasta in zip(ys, sys.argv[1:]):
    x = range(1, len(y)+1)
    plt.plot(x, y, label=fasta)
plt.xlabel('Contig #')
plt.xscale('log')
plt.ylabel('Length')
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()