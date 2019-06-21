import json
import itertools
import csv
import re

# Parse output
def read_output(filename):
    with open(filename, 'r') as json_file:
        data = json.load(json_file)

    n_covered = data['n_covered']
    total_lib_size = data['total_lib_size']
    
    return (n_covered, total_lib_size)
    

# Parse time file
def read_time(filename):
    with open(filename, 'r') as time_file:
        time_lines = time_file.read().split('\n')

    for line in time_lines:
        line = line.strip()
        if line.startswith('Elapsed'):
            time = line.replace('Elapsed (wall clock) time (h:mm:ss or m:ss):', '').strip()

            split_time = time.split(':')

            if len(split_time) == 2:
                h = 0
                m = float(split_time[0])
                s = float(split_time[1])

            elif len(split_time) == 3:
                h = float(split_time[0])
                m = float(split_time[1])
                s = float(split_time[2])

            t = (h * 3600) + (m * 60) + s

        elif line.startswith('Maximum'):
            mem = float(line.replace('Maximum resident set size (kbytes):', '').strip())
            
    return (t, mem)
            

# Parse log file
def read_log(filename):
    with open(filename, 'r') as handle:
        test = handle.readlines()

    data = []

    keeper = False

    for line in test:
        line = line.strip()
        if line.startswith('Expl Unexpl'):
            keeper = True

        if re.search(r'[0-9]+s$', line):

            if keeper:
                splitline = line.split()
                data_keep = []
                data_keep.append(-1 * float(splitline[-5]))
                data_keep.append(-1 * float(splitline[-4]))
                data_keep.append(float(splitline[-3].replace('%', '')) if splitline[-3] != '-' else 'NA')
                data_keep.append(float(splitline[-1].replace('s', '')))

                data.append(data_keep)
                
    return data


if __name__ == '__main__':
    
    # Handle SwiftLib comparison datasets
    dirs = ['ilp', 'sl']
    lib_sizes = ['100000', '1000000', '10000000', '100000000', '1000000000']
    sublibs = ['1', '2']
    
    base_name = '../results/sl_comparison/{}/gfp_239_{}_{}'
    
    sl_results = [['method', 'lib_limit', 'sublibs', 'n_covered', 'total_lib_size', 'time', 'max_mem']]
    
    for i in itertools.product(dirs, lib_sizes, sublibs):
        
        filename = base_name.format(i[0], i[1], i[2])
        
        if i[0] == 'ilp':
            method = 'DeCoDe'
        else:
            method = 'SwiftLib'
        
        out = read_output(filename + '.json')
        time = read_time(filename + '.time')
        
        data_line = [method, i[1], i[2], out[0], out[1], time[0], time[1]]
        
        sl_results.append(data_line)
        
    with open('../results/sl_comparison/results.csv', 'w') as results_file:
        writer = csv.writer(results_file)
        writer.writerows(sl_results)
        
    # Handle multi sublibrary data 
    sublibs = ['1', '2', '3', '4', '8', '12']
    
    base_name = '../results/multi_sublib/gfp_exclude_long_10000000_{}'
    
    ms_results = [['method', 'lib_limit', 'sublibs', 'n_covered', 'total_lib_size', 'time', 'max_mem']]
    ms_logs = [['sublibs', 'solution', 'bound', 'gap', 'time']]
    
    for i in sublibs:
        
        filename = base_name.format(i)
        
        method = 'DeCoDe'
        
        out = read_output(filename + '.json')
        time = read_time(filename + '.time')
        log = read_log(filename + '.log')
        
        for array in log:
            array.insert(0, i)
            ms_logs.append(array)
        
        data_line = [method, 10000000, i, out[0], out[1], time[0], time[1]]
        
        ms_results.append(data_line)
        
    with open('../results/multi_sublib/results.csv', 'w') as results_file:
        writer = csv.writer(results_file)
        writer.writerows(ms_results)
    
    with open('../results/multi_sublib/log.csv', 'w') as results_file:
        writer = csv.writer(results_file)
        writer.writerows(ms_logs)

        