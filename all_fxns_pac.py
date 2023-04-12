# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 08:46:57 2022

@author: arina
"""

import os, re
import pandas as pd

root_dir = 'C:\\Users\\arina\\Desktop\\arina_hdeeg_pac\\testfxn\\'
ppt_id = ['MCI03JP', 'MCI02JG']
stage = 'nrem'


# create sub-directories for ppt in results folder 
def make_direct(ppt_id, root_dir, stage):
    stage = stage.lower()
    outpath = root_dir + 'pac_{}\\results\\'.format(stage)
    print('------------------------------------------\n')

    for p in ppt_id:
        print('Creating results folder for {}....\n'.format(p))
        try:
            os.makedirs(outpath+p)
        except FileExistsError:
            pass
                
    print("\nFolders successfully created for all participants.\n")
    

def ppa(ppt_id, root_dir, stage, spindle_type):
    stage = stage.lower()
    spindle_type = spindle_type.lower()
    inpath = root_dir + 'pac_{}\\source\\'.format(stage)
    print('------------------------------------------\n')

    for p in ppt_id:
        # get all files from ppt directory
        files = os.listdir(inpath+'{}\\{}\\'.format(p, spindle_type))
        print('Summarising preferred phase angle output for {}....\n'.format(p))
        
        master_df = pd.DataFrame(columns=['channel', 'radians', 'degrees'], dtype=object)
        
        for f in files: 
            # detect all channels per file
            file_chs = re.findall('[E]\d+', f)
            first_ch = file_chs[0]
            
            # read file as df, clean col names and add channel col 
            df = pd.read_csv(inpath+'{}\\{}\\{}'.format(p, spindle_type, f))
            df['channel'] = file_chs
            df.rename(columns={ first_ch+'radians':'radians', first_ch+'degrees':'degrees' }, inplace=True)
            df = df[['channel', 'radians', 'degrees']]
            
            # concat current df to master df
            master_df = pd.concat([master_df, df])
            
        
        # further clean cols: channel and radians 
        master_df['channel'] = master_df['channel'].str.replace('E', '').astype('int')
        master_df['radians'] = master_df['radians'].str.replace('[', '').str.replace(']', '').astype('float')
        
        master_df.sort_values(by=['channel'], inplace=True)
        master_df = master_df.set_index('channel')
        
        master_df.to_csv(root_dir+'pac_{}\\results\\{}\\{}_preferred_phase_{}spd_{}.csv'.format(stage, p, p, spindle_type, stage))

    print("\nPreferred phase angles successfully summarised for all participants.\n")


