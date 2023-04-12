# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 09:16:59 2022

@author: arina
"""

import os, re
import pandas as pd

from all_fxns_pac import *

######## MANUAL INPUT ########

root_dir = 'C:\\Users\\arina\\Desktop\\arina_hdeeg_pac\\testfxn\\'
ppt_id = ['MCI03JP', 'MCI02JG']
stage = 'nrem'

#######################################################################################################

stage = stage.lower()

for p in ppt_id: 
    check_dir = root_dir + 'results\\pac_{}\\results\\{}'.format(stage, p)
    
    if not os.path.exists(check_dir):
        make_direct(ppt_id, root_dir, stage)
        

def get_ppa_fast_spindles():
    spindle_type = 'fast'
    ppa(ppt_id, root_dir, stage, spindle_type)

def get_ppa_slow_spindles():
    spindle_type = 'slow'
    ppa(ppt_id, root_dir, stage, spindle_type)
