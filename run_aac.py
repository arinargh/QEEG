# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 19:36:32 2022

@author: arina
"""

import csv, os, time
import scipy.stats as stats, numpy as np, pandas as pd
from datetime import timedelta

from all_fxns import *



######## MANUAL INPUT ########

ppt_id = 'MCI02JG'
root_dir = 'C:\\Users\\arina\\Desktop\\arina_hdeeg_aac_compare\\aac_nrem\\'

channels = ['E' + str(i) for i in range(0,3) if i not in [0]]

stage = 'NREM' # NREM N2 SWS

#######################################################################################################

check_dir = root_dir + 'results\\cleaned\\{}\\spindle\\fast'.format(ppt_id)

if not os.path.exists(check_dir):
    make_direct(ppt_id, root_dir)
    

def slowwaves():    
    filetype = 'slowwave'
    merge_files(filetype, ppt_id, channels, root_dir)
    
    df_nophys = remove_outliers_physio(filetype, ppt_id, channels, root_dir)
    df_annot = check_outliers_stats(filetype, ppt_id, channels, root_dir)
    df_nostats = remove_outliers_stats(df_annot, filetype, ppt_id, channels, root_dir)
    
    output_summary(df_nophys, df_nostats, filetype, ppt_id, channels, root_dir)
    
    hypno = load_hypno(ppt_id, root_dir)
    split_events(hypno, filetype, ppt_id, channels, root_dir)
    

def fast_spindles():
    filetype = 'spindle'
    spindle_type = 'fast'
    merge_files(filetype, ppt_id, channels, root_dir, spindle_type)
    
    df_nophys = remove_outliers_physio(filetype, ppt_id, channels, root_dir, spindle_type)
    df_annot = check_outliers_stats(filetype, ppt_id, channels, root_dir, spindle_type)
    df_nostats = remove_outliers_stats(df_annot, filetype, ppt_id, channels, root_dir, spindle_type)
    
    output_summary(df_nophys, df_nostats, filetype, ppt_id, channels, root_dir, spindle_type)
    
    hypno = load_hypno(ppt_id, root_dir)
    split_events(hypno, filetype, ppt_id, channels, root_dir, spindle_type)
    

def slow_spindles():
    filetype = 'spindle'
    spindle_type = 'slow'
    merge_files(filetype, ppt_id, channels, root_dir, spindle_type)
    
    df_nophys = remove_outliers_physio(filetype, ppt_id, channels, root_dir, spindle_type)
    df_annot = check_outliers_stats(filetype, ppt_id, channels, root_dir, spindle_type)
    df_nostats = remove_outliers_stats(df_annot, filetype, ppt_id, channels, root_dir, spindle_type)
    
    output_summary(df_nophys, df_nostats, filetype, ppt_id, channels, root_dir, spindle_type)
    
    hypno = load_hypno(ppt_id, root_dir)
    split_events(hypno, filetype, ppt_id, channels, root_dir, spindle_type)
    


def get_aac_fast_spindles():
    spindle_type = 'fast'
    
    check_summfile(spindle_type, ppt_id, root_dir, channels, stage)
    hypno = load_hypno(ppt_id, root_dir)
    stage_duration = calc_duration(hypno, stage)
    
    aac(stage_duration, spindle_type, ppt_id, root_dir, channels, stage)
    
    
def get_aac_slow_spindles():
    spindle_type = 'slow'
        
    check_summfile(spindle_type, ppt_id, root_dir, channels, stage)
    hypno = load_hypno(ppt_id, root_dir)
    stage_duration = calc_duration(hypno, stage)
    
    aac(stage_duration, spindle_type, ppt_id, root_dir, channels, stage)