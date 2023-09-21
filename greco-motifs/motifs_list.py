with open("motifs_names_list.txt") as infile, open("motifs_list.txt", 'w') as outfile:
    for line in infile:
        tf_name = line.strip().split("@")[0]
        pwm_name = line.strip()
        outfile.write(f"./{tf_name}/{pwm_name}\n")
