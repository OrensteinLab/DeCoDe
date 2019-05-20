from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.keys import Keys
import pandas as pd

def run_swift_lib(input_string, input_type, lib_size_limit, oligo_limit):
    # Setup selenium 
    options = Options()
    options.headless = True
    driver = webdriver.Firefox(options=options)
    
    # Get the SwiftLib website
    driver.get('http://rosettadesign.med.unc.edu/SwiftLib/')
    
    # Get the tabs
    manual_tab = driver.find_element_by_id('ui-id-1')
    csv_tab = driver.find_element_by_id('ui-id-2')
    fasta_tab = driver.find_element_by_id('ui-id-3')
    

    # Get the buttons/clickable elements for fasta input
    fasta_input = driver.find_element_by_id('fasta')
    fasta_update_table_button = driver.find_element_by_id('table_from_fasta')

    # Get the buttons/clickable elements for csv input
    csv_input = driver.find_element_by_id('csvaacounts')
    csv_update_table_button = driver.find_element_by_id('table_from_csv')

    # Get references to other buttons and fields
    multiple_dcs_button = driver.find_element_by_id('allow_mult_dcs')
    libsize_field = driver.find_element_by_id('libsize_upper')
    max_primers_field = driver.find_element_by_id('max_primers_total')
    launch_button = driver.find_element_by_id('launchbutton')

    # Activate multiple codons
    multiple_dcs_button = driver.find_element_by_id('allow_mult_dcs')
    multiple_dcs_button.click()

    # Handle FASTA input
    if input_type == 'fasta':
        fasta_tab.click()
        fasta_input.click()
        fasta_input.send_keys(input_string)
        fasta_update_table_button.click()

    # Handle CSV input
    elif input_type == 'csv':
        csv_tab.click()
        csv_input.click()
        csv_input.send_keys(Keys.CONTROL + 'a')
        csv_input.send_keys(Keys.DELETE)
        csv_input.send_keys(input_string)
        csv_update_table_button.click()

    # Go back to the main tab
    manual_tab.click()
    
    # Enter the max library size 
    libsize_field.click()
    libsize_field.send_keys(lib_size_limit)
    
    # Enter the max oligo count
    max_primers_field.click()
    max_primers_field.send_keys(oligo_limit)
    
    # Run SwiftLib
    launch_button.click()

    # Get the table
    table_output = driver.find_element_by_id('scrollhere')
    html_out = table_output.get_attribute('outerHTML')
    
    # Read the table
    swift_lib_codon_table = pd.read_html(html_out)[0]
    
    return swift_lib_codon_table

