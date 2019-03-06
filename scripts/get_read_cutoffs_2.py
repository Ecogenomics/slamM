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
    y = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0]
    print(bp_sum)
    max_read = 8000
    cutoffs = [i*bp_sum for i in y]
    print(cutoffs)
    read_cutoffs = []
    index = 0
    bp_sum = 0
    bp_sum_max_read = 0
    for i in read_lengths:
        if index == len(cutoffs):
            break
        bp_sum += i
        if i >= max_read:
            bp_sum_max_read  += i
        if bp_sum >= cutoffs[index]:
            index += 1
            print(i, index)
            if i <= 500:
                read_cutoffs.append(500)
                break
            else:
                read_cutoffs.append(i)

print(read_cutoffs)
with open(snakemake.output[0], 'w') as o:
    for num, i in enumerate(read_cutoffs):
        if i > max_read:
            fraction = cutoffs[num] / bp_sum_max_read
            o.write(str(max_read) + '\t' + str(fraction) + '\n')
        else:
            o.write(str(i) + '\t1.00\n')


