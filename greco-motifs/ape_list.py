with open('tf_list.txt') as infile, open('ape_list.txt', 'w') as outfile:
    for line in infile:
        outfile.write(f"java -cp ./ape.jar ru.autosome.ape.PrecalculateThresholds ./motifs_organized_v5_pwm/{line.strip()} ./motifs_organized_thr --background 0.25,0.25,0.25,0.25 --pvalues 1e-15,1.0,1.01,mul --discretization 10000 --silent\n")
