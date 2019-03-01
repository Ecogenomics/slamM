import sys
import os

paf = sys.argv[1]
bin_folder = sys.argv[2]




bin_dict = {}
bin_count = {}
for i in os.listdir(bin_folder):
    bin = i.split('.')[1]
    bin_count[bin] = 0
    with open(os.path.join(bin_folder, i)) as f:
        for line in f:
            if line.startswith('>'):
                bin_count[bin] += 1
                bin_dict[line.split()[0][1:]] = bin



buffer = 50
cov_dict = {}
with open(paf) as aln:
    edge_list = []
    last_q_name = None
    hits = []
    for line in aln:
        qname, qlen, qstart, qstop, strand, ref, rlen, rstart, rstop = line.split()[:9]
        qlen, qstart, qstop, rlen, rstart, rstop = map(int, [qlen, qstart, qstop, rlen, rstart, rstop])
        if not ref in cov_dict:
            cov_dict[ref] = 0
        cov_dict[ref] += (rstop - rstart) / rlen
        if not last_q_name is None and qname != last_q_name:
            for i in hits:
                qname1, qlen1, qstart1, qstop1, strand1, ref1, rlen1, rstart1, rstop1 = i
                mapped_start, mapped_end = False, False
                if qstart1 < buffer and strand == '-' and rstart1 < buffer:
                    mapped_start = True
                if qstart1 < buffer and strand == '+' and rstop1 > rlen1 - buffer:
                    mapped_start = True
                if qstop1 < qlen1 - buffer and strand == '-' and rstop1 > rlen1 - buffer:
                    mapped_end = True
                if qstop1 < qlen1 - buffer and strand == '+' and rstart1 < buffer:
                    mapped_end = True
                for j in hits:
                    if i == j:
                        continue
                    qname2, qlen2, qstart2, qstop2, strand2, ref2, rlen2, rstart2, rstop2 = j
                    if qstart2 < buffer and strand == '-' and rstart2 < buffer and mapped_end:
                        edge_list.append((ref1, ref2, qname1))
                    if qstart2 < buffer and strand == '+' and rstop2 > rlen2 - buffer and mapped_end:
                        edge_list.append((ref1, ref2, qname1))
                    if qstop2 < qlen2 - buffer and strand == '-' and rstop2 > rlen2 - buffer and mapped_start:
                        edge_list.append((ref1, ref2, qname1))
                    if qstop2 < qlen2 - buffer and strand == '+' and rstart2 < buffer:
                        edge_list.append((ref1, ref2, qname1))
            hits = []
        last_q_name = qname
        hits.append((qname, qlen, qstart, qstop, strand, ref, rlen, rstart, rstop))


for i in bin_dict:
    if bin_dict[i] == "3":
        print(i, cov_dict[i])

print('doaa')
for i in bin_dict:
    if bin_dict[i] == "1":
        print(i, cov_dict[i])

sys.exit()


edge_dict = {}

for i in edge_list:
    b_ref, o_ref, rname = i
    if not b_ref in edge_dict:
        edge_dict[b_ref] = {}
    if not o_ref in edge_dict:
        edge_dict[o_ref] = {}
    if not o_ref in edge_dict[b_ref]:
        edge_dict[b_ref][o_ref] = 0
    if not b_ref in edge_dict[o_ref]:
        edge_dict[o_ref][b_ref] = 0
    edge_dict[o_ref][b_ref] += 1
    edge_dict[b_ref][o_ref] += 1

min_edge = 2

edge_list = []
for i in edge_dict:
    for j in edge_dict[i]:
        if edge_dict[i][j] >= min_edge:
            edge_list.append((i, j))


counts = {}
to_bin = {}
for i, j in edge_list:
    if i in counts:
        counts[i] += 1
        if j in bin_dict and i in bin_dict and bin_dict[i] != bin_dict[j]:
            to_bin[i].append(bin_dict[j])
    else:
        counts[i] = 1
        to_bin[i] = []
        if j in bin_dict and i in bin_dict and bin_dict[i] != bin_dict[j]:
            to_bin[i].append(bin_dict[j])




filter_list = set()
for i in counts:
    print(i, counts[i], cov_dict[i], to_bin[i])
    if counts[i] > 4:
        print('ding')
        filter_list.add(i)

new_edge_list = []
for i in edge_list:
    if not i[0] in filter_list and not i[1] in filter_list:
        new_edge_list.append(i)

edge_list = new_edge_list
unbinned_bins = {}
bin_edges = {}
for o, b in edge_list:
    if o in bin_dict:
        obin = bin_dict[o]
    else:
        obin = 'unbinned'
    if b in bin_dict:
        bbin = bin_dict[b]
    else:
        bbin = 'unbinned'
    if obin == bbin:
        pass
    elif obin == 'unbinned':
        if not o in unbinned_bins:
            unbinned_bins[o] = {}
        if not bbin in unbinned_bins[o]:
            unbinned_bins[o][bbin] = 0
        unbinned_bins[o][bbin] += 1
    elif bbin != 'unbinned':
        if not obin in bin_edges:
            bin_edges[obin] = {}
        if not bbin in bin_edges[obin]:
            bin_edges[obin][bbin] = 0
        bin_edges[obin][bbin] += 1
        if obin == '1' and bbin == '3':
            print(o, b, 'ding')
        elif obin == '3' and bbin == '1':
            print(o, b, 'dong')

for i in bin_edges:
    for j in bin_edges[i]:
        print(i, j, bin_edges[i][j], bin_edges[i][j], bin_count[i], bin_count[j])



