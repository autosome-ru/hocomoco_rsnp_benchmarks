with open('./motifs_list.txt') as infile, \
        open('selex_benchmark_list.txt', 'w') as outfile:
    for line in infile:
        outfile.write("./selex_pwm_benchmark.py " + f"{line.strip()[:-4]}\n")
