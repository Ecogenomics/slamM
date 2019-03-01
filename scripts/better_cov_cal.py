import sys
import numpy

paf_file = sys.argv[1]


cov_dict = {}
count = 0
with open(paf_file) as f:
    for line in f:
        count += 1
        if count % 10000 == 0:
            print(count)
        qname, qlen, qstart, qstop, strand, ref, rlen, rstart, rstop = line.split()[:9]
        qlen, qstart, qstop, rlen, rstart, rstop = map(int, [qlen, qstart, qstop, rlen, rstart, rstop])
        if not ref in cov_dict:
            cov_dict[ref] = numpy.zeros(rlen)
        for i in range(rstart, rstop):
            cov_dict[ref][i] += 1


for i in cov_dict:
    cov = cov_dict[i]
    average = sum(cov) / len(cov)
    var_tot = 0
    for j in cov:
        var_tot += (j - average) ** 2
    variance = var_tot / (len(cov) - 1)
    print(i, average, variance)