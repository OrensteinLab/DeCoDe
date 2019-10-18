import json

sl_p1_sol2 = [
        'ADEGHIKLMNPQRSTV',
        'ADEGKLMNRSTVW',
        'AILSTV',
        'AS',
        'ALMSTV',
        'AGST',
        'ADEGKNRST',
        'DEHKLMNQRSTW',
        'ADEGHIKLMNPQRSTV'
]

with open('results/sl_comparison/ilp/1xbi_320000000_4.json', 'r') as jsonfile:
    dc_lib_dict = json.load(jsonfile)
    
    
sl_n_covered = 0

for seq in dc_lib_dict['sequences']:
    match_count = 0
    for idx, aa in enumerate(list(seq)):
        if aa in sl_p1_sol2[idx]:
            match_count += 1
            
    if match_count == 9:
        sl_n_covered += 1
        
if __name__ == '__main__':
    print('DeCoDe n covered:\t{}\nSwiftLib n covered:\t{}'.format(dc_lib_dict['n_covered'], sl_n_covered))
            