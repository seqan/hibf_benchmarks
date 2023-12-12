## INPUT:
# 0         1           2       3               4       5               6           7                           8
# KEY	    SUBKEY	    LEVEL	BIT_SIZE	    IBFS	AVG_LOAD_FACTOR	TBS_TOO_BIG	AVG_TBS_TOO_BIG_ELEMENTS	AVG_MAX_ELEMENTS
# alpha	    alpha=0.5	0	    66459366336	    1	    78.60	        1	        6999443	                    212522373
# tmax=1024	U+R	        0	    128499762176	1	    99.76	        1018	    16593073	                15880273

## OUTPUT:
# [Subkeys]
# [Level 0]
# [Level 1]
# [Level 2]
# [Level 3]
# [Avg Load Factor]

import os
import csv

input_path = "/Users/leoxnard/Documents/JupiterNotebooks/data/all.sizes"
output_directory = "/Users/leoxnard/Documents/JupiterNotebooks/sizes"
keys = ['k', 'alpha', 'hash', 'r-fpr', 'none', 'U+R', 'U']
repeat = 3

try:
    os.mkdir(output_directory)
    print("Directory %s created" % output_directory)
except OSError:
    print("files replaced in directory %s" % output_directory)

for key in keys:
    output_path = f'{output_directory}/{key}.csv'
    subkeys = []
    bit_size_levels = [[] for _ in range(4)]
    avg_load_factor_levels = [[] for _ in range(4)]
    with open(input_path, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        expected_level_counter = 0
        for line in reader:
            if line[1] == key or line[0] == key:
                if int(line[2]) != expected_level_counter:
                    for i in range(expected_level_counter, 4):
                        avg_load_factor_levels[i].append(0)
                        bit_size_levels[i].append(0)
                    expected_level_counter = 0
                if line[1] == key and expected_level_counter == 0: subkeys.append(float(line[0].split('=')[1]))
                elif expected_level_counter == 0: subkeys.append(float(line[1].split('=')[1]))
                bit_size_levels[int(line[2])].append(round(float(line[3]) / (8*1024*1024*1025), 2))
                avg_load_factor_levels[int(line[2])].append(round(float(line[5])/100, 2))
                expected_level_counter += 1
                if expected_level_counter == 4: expected_level_counter = 0
        for i in range(expected_level_counter, 4):
            avg_load_factor_levels[i].append(0)
            bit_size_levels[i].append(0)

    subkeys, bit_size_level_0, bit_size_level_1, bit_size_level_2, bit_size_level_3, avg_load_factor_level_0, avg_load_factor_level_1, avg_load_factor_level_2, avg_load_factor_level_3  = zip(*sorted(zip(subkeys, bit_size_levels[0], bit_size_levels[1], bit_size_levels[2], bit_size_levels[3], avg_load_factor_levels[0], avg_load_factor_levels[1], avg_load_factor_levels[2], avg_load_factor_levels[3])))
    all_lists = [subkeys,bit_size_level_0, bit_size_level_1, bit_size_level_2, bit_size_level_3, avg_load_factor_level_0, avg_load_factor_level_1, avg_load_factor_level_2, avg_load_factor_level_3]
    with open(output_path, "w", newline='') as file:
        writer = csv.writer(file)
        for list_ in all_lists:
            writer.writerow(list_)