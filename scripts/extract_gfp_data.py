import requests
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
import os
import time
from tqdm import *
import json

if __name__ == '__main__':
    
    # Setup Firefox to work with selenium for FPbase
    options = Options()
    options.headless = True
    driver = webdriver.Firefox(options=options)
    driver.get('https://www.fpbase.org/protein/avgfp')
    time.sleep(5)
    source = driver.page_source
    driver.quit()
    
    # Parse the HTML
    soup = BeautifulSoup(source, 'html.parser')
    
    # Get all the fluorescent proteins in the tree
    fps = {}

    for node in soup.find_all('g', {'class': 'node'}):
        link = node.find('a')['xlink:href']
        safe_name = os.path.basename(os.path.normpath(link))
        text = node.find('text')
        if text:
            name = text.contents[0]
        else:
            name = 'avGFP'
        fps[safe_name] = {'link': link,
                          'name': name}
    
    # Get the data for each individual protein
    for protein in tqdm(fps):
        resp = requests.get('https://www.fpbase.org/protein/{}'.format(protein))
        soup = BeautifulSoup(resp.content, 'html.parser')

        seq_div = soup.find('div', {'class': 'sequence'})
        sequence = seq_div.find('div', {'class': 'aminosequence'}).contents[0]

        fps[protein]['sequence'] = sequence

        deriv_div = seq_div.find('p', {'class': 'text-muted'})

        if deriv_div:
            parent = os.path.basename(os.path.normpath(deriv_div.find('a')['href']))
            mutations = deriv_div.find('span', {'class': 'mut-rel-root'}).contents[0].split('/')

            fps[protein]['parent'] = parent
            fps[protein]['mutations'] = mutations
            
    # Write out the sequence datasets
    with open('../examples/gfp/gfp.fa', 'w') as handle:
        for fp in fps:
            handle.write('>{}\n'.format(fp))
            handle.write('{}\n'.format(fps[fp]['sequence']))
        
    with open('../examples/gfp/gfp_exclude_long.fa', 'w') as handle:
        for fp in fps:
            seq = fps[fp]['sequence']
            if len(seq) in [238, 239]:
                handle.write('>{}\n'.format(fp))
                handle.write('{}\n'.format(seq))
            
    with open('../examples/gfp/gfp_239.fa', 'w') as handle:
        for fp in fps:
            seq = fps[fp]['sequence']
            if len(seq) == 239:
                handle.write('>{}\n'.format(fp))
                handle.write('{}\n'.format(seq))
        
    # Write all of the data out to a JSON file
    with open('../aux_data/gfp.json', 'w') as handle:
        json.dump(fps, handle)
