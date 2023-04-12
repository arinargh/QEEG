# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:17:55 2023

@author: arina
"""

import pandas as pd, numpy as np
import os, re

# step1: amalgamate all channels into a master dataset

root_dir = 'C:\\Users\\arina\\OneDrive\\Data_Portfolio\\qeeg\\'

### create master summfile - PSA 
def amalg_psa(root_dir, fillna=None):
    inpath = root_dir + 'raw_output\\'
    outpath = root_dir + 'cleaned_output\\'
    
    # list of PSG macro-data in PSA masterfile
    dup_cols = ['TIB_min', 'TST_min', 'SOL_min', 'ROL_min', 'WASO_min',
           'SE_percent', 'REM_min', 'S1_min', 'S2_min', 'S3_min', 'S4_min',
           'NREM_min', 'SWS_min', 'REM_percent', 'S1_percent', 'S2_percent',
           'S3_percent', 'S4_percent', 'NREM_percent', 'SWS_percent']
    
    file_list = [file for file in os.listdir(inpath) if 'psa_summary_' in file]

    ch_list = list()
    
    print('=====\nPSA summary files found for the following channels:')
    
    for file in file_list: 
        ch_catch = re.search('psa_summary_(\w\w)', file)
        ch_list.append(ch_catch.group(1))
        print("*{}".format(ch_catch.group(1)))
    
    print('=====\n')
    
    df_dict = dict()
    
    for f in file_list: 
        ch_df = pd.read_excel(inpath + f, na_values='.')
        
        df_dict[re.search('psa_summary_(\w\w)', f).group(1)] = ch_df
    
    # create master summary file
    first_ch = list(df_dict.keys())[0]
    master_df = df_dict[first_ch]
    del df_dict[first_ch]
    
    for ch,df in df_dict.items():
        # remove dup columns for each ch df: 
        df_dict[ch].drop(dup_cols, axis=1, inplace=True)
            
        master_df = master_df.merge(df_dict[ch], how='outer', on='PID')
        
    
    # insert divider columns for easier navigation in SPSS
    master_df.insert(21, 'END_QEEG_PSG_DATA', '')
    master_df.insert(1, 'START_QEEG_PSG_DATA', '')
    
    master_df.insert(23, 'START_QEEG_PSA_DATA', '')
    master_df['END_QEEG_PSA_DATA'] = ''
    
    # drop S3 (N3) columns - repeat of SWS data
    # drop S4 (N4) columns - empty 
    drop_cols = [x for x in list(master_df.columns) if 'S3' in x or 'S4' in x]
    
    master_df.drop(columns=drop_cols, inplace=True)
    
    # export master summary file
    if fillna==True: 
        master_df.fillna(value=999, inplace=True)
        print("Missing values replaced with 999\n")
        
    master_df.to_csv(outpath+'\\psa_master_summfile.csv', index=False)
    
    print("Done! Please check the file: 'psa_master_summfile.csv'\n({})\n".format(outpath))

### create master summfile - SPD
def amalg_spd(root_dir, fillna=None):
    inpath = root_dir + 'raw_output\\'
    outpath = root_dir + 'cleaned_output\\'
    
    
    file_list = [file for file in os.listdir(inpath) if 'spd_summary_' in file]

    ch_list = list()
    
    print('=====\nSpindle summary files found for the following channels:')

    for file in file_list: 
        ch_catch = re.search('spd_summary_(\w\w)', file)
        ch_list.append(ch_catch.group(1))
        print("*{}".format(ch_catch.group(1)))

    print('=====\n')

    df_dict = dict()
    
    for f in file_list: 
        ch_df = pd.read_excel(inpath + f, na_values='0')
        ch_df.rename(columns={'Threshold':'{}_Threshold'.format(re.search('spd_summary_(\w\w)', f).group(1))}, inplace=True)
        df_dict[re.search('spd_summary_(\w\w)', f).group(1)] = ch_df
    
    # create master summary file
    first_ch = list(df_dict.keys())[0]
    master_df = df_dict[first_ch]
    del df_dict[first_ch]
    
    for ch,df in df_dict.items():            
        master_df = master_df.merge(df_dict[ch], how='outer', on='PID')
    
    # insert divider columns 
    master_df.insert(1, 'START_QEEG_SPD_DATA', '')
    master_df['END_QEEG_SPD_DATA'] = ''
    
    # export master summary file
    if fillna==True: 
        master_df.fillna(value=999, inplace=True)
        print("Missing values replaced with 999\n")
    
    master_df.to_csv(outpath+'\\spd_master_summfile.csv', index=False)

    print("Done! Please check the file: 'spd_master_summfile.csv'\n({})\n".format(outpath))


###############################################################################

### create annot file 

def export_annot(root_dir):
    annot_df = pd.DataFrame()
    
    # get all ppt rows 
    psa = root_dir + 'cleaned_output\\psa_master_summfile.csv'
    spd = root_dir + 'cleaned_output\\spd_master_summfile.csv'
    
    psa_df = pd.DataFrame()
    spd_df = pd.DataFrame()
    
    try: 
        psa_df = pd.read_csv(psa)
    except:
        print("No PSA master summary file found. Skipping...\n")
    
    try:
        spd_df = pd.read_csv(spd)
    except:
        print("No spindle master summary file found. Skipping...\n")
        
    if not psa_df.empty and not spd_df.empty:
        annot_df['PID'] = pd.concat([psa_df, spd_df], axis=0, join='outer', ignore_index=True)['PID']
        annot_df.drop_duplicates(inplace=True)
    elif psa_df.empty:
        annot_df['PID'] = spd_df['PID']
    elif spd_df.empty:
        annot_df['PID'] = psa_df['PID']
    
    
    # get all ch columns
    file_list = [file for file in os.listdir(root_dir + 'raw_output\\') if 'summary_' in file]
    
    for f in file_list:
        annot_df['{}_Quality'.format(re.search('summary_(\w\w)', f).group(1))] = ''
        
    
    # export annot_file
    annot_df.to_csv(root_dir + 'cleaned_output\\annot_bcr.csv', index=False)
    
    print("Annotations file created. Please mark channels to be excluded in the file: 'annot_bcr.csv \n({})".format(root_dir+'cleaned_output'))
    
    
###############################################################################

### perform BCR

def bcr_psa(root_dir):
    inpath = root_dir + 'cleaned_output\\'
    
    # read BCR and PSA files
    try:
        annot_file = pd.read_csv(inpath + '\\annot_bcr.csv', na_values='')
    except FileNotFoundError:
        print("\nERROR: BCR annotations file missing\n")
            
    try: 
        vqc_df = pd.read_csv(inpath + '\\psa_master_summfile_bcr.csv')
    except FileNotFoundError:
        vqc_df = pd.read_csv(inpath + '\\psa_master_summfile.csv')

    # define channel list and BCR dictionary
    file_list = [file for file in os.listdir(root_dir + 'raw_output\\') if 'psa_summary_' in file]
    
    ch_list = list()
    
    print('=====\nBCR performed on PSA data for the following channels:')

    for file in file_list: 
        ch_catch = re.search('summary_(\w\w)', file)
        ch_list.append(ch_catch.group(1))
        print("*{}".format(ch_catch.group(1)))

    print('=====\n')
    
    bcr_dict = dict()
    
    # add to BCR dictionary
    
    for ch in ch_list: 
        bcr_list = annot_file[ annot_file['{}_Quality'.format(ch)].str.upper()=='BAD' ]['PID'].to_list()
        bcr_dict[ch] = bcr_list
    
    # remove bad channels
    
    for ch, ppt_list in bcr_dict.items():
        ch_cols = list(vqc_df.filter(like=ch).columns)
        for ppt in ppt_list:
            vqc_df.loc[ vqc_df['PID']==ppt, ch_cols ] = 111111
    
    vqc_df.to_csv(inpath + 'psa_master_summfile_bcr.csv', index=False)
    
    print("Bad quality channel values replaced with 111111\n")
    print("Done! Please check the file: 'psa_master_summfile_bcr.csv'\n({})\n".format(inpath))



def bcr_spd(root_dir):
    inpath = root_dir + 'cleaned_output\\'
    
    # read BCR and SPD files
    try:
        annot_file = pd.read_csv(inpath + '\\annot_bcr.csv', na_values='')
    except FileNotFoundError:
        print("\nERROR: BCR annotations file missing\n")
            
    try: 
        vqc_df = pd.read_csv(inpath + '\\spd_master_summfile_bcr.csv')
    except FileNotFoundError:
        vqc_df = pd.read_csv(inpath + '\\spd_master_summfile.csv')

    # define channel list and BCR dictionary
    file_list = [file for file in os.listdir(root_dir + 'raw_output\\') if 'spd_summary_' in file]
    
    ch_list = list()
    
    print('=====\nBCR performed on PSA data for the following channels:')

    for file in file_list: 
        ch_catch = re.search('summary_(\w\w)', file)
        ch_list.append(ch_catch.group(1))
        print("*{}".format(ch_catch.group(1)))

    print('=====\n')
    
    bcr_dict = dict()
    
    # add to BCR dictionary
    
    for ch in ch_list: 
        bcr_list = annot_file[ annot_file['{}_Quality'.format(ch)].str.upper()=='BAD' ]['PID'].to_list()
        bcr_dict[ch] = bcr_list
    
    # remove bad channels
    
    for ch, ppt_list in bcr_dict.items():
        ch_cols = list(vqc_df.filter(like=ch).columns)
        for ppt in ppt_list:
            vqc_df.loc[ vqc_df['PID']==ppt, ch_cols ] = 111111
    
    vqc_df.to_csv(inpath + 'spd_master_summfile_bcr.csv', index=False)
    
    print("Bad quality channel values replaced with 111111\n")
    print("Done! Please check the file: 'spd_master_summfile_bcr.csv'\n({})\n".format(inpath))