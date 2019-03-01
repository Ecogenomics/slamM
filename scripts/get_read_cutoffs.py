with open(snakemake.input[0]) as f:
    read_lengths = []
    f.readline()
    f.readline()
    f.readline()
    last_read_no = int(f.readline().split()[1])
    for line in f:
        read_length, read_no = map(int, line.split()[:2])
        if read_no != last_read_no:
            for x in range(read_no, last_read_no):
                read_lengths.append(read_length - 1)
            last_read_no = read_no
    for x in range(0, last_read_no):
        read_lengths.append(read_length)
    read_lengths.sort(reverse=True)
    bp_sum = sum(read_lengths)
    y = int(snakemake.params[0])
    cutoffs = [x/10000 * bp_sum for x in range(10000//y, 9999, 10000//y)]
    read_cutoffs = []
    for i in read_lengths:
        if cutoffs == []:
            break
        bp_sum -= i
        if bp_sum < cutoffs[-1]:
            cutoffs.pop()
            read_cutoffs.append(i)
    if read_cutoffs[-1] > 500:
        read_cutoffs.append(500)


with open(snakemake.output[0], 'w') as o:
    o.write("bases=" + str(bp_sum))
    for i in read_cutoffs:
        o.write(str(i) + '\n')


