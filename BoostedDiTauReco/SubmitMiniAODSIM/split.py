lines_per_file = 10
smallfile = None
with open('QCD_Pt-15to7000.txt') as bigfile:
    for lineno, line in enumerate(bigfile):
        if lineno % lines_per_file == 0:
            if smallfile:
                smallfile.close()
            small_filename = 'QCD_Pt-15to7000_{}.txt'.format((lineno + lines_per_file)/10)
            smallfile = open(small_filename, "w")
        smallfile.write(line)
    if smallfile:
        smallfile.close()
