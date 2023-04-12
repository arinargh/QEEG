# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 18:41:23 2022

@author: arina
"""

import csv, os, time
import scipy.stats as stats, numpy as np, pandas as pd
from datetime import timedelta


# create sub-directories for ppt in results folder 
def make_direct(ppt_id, root_dir):
    outpath = root_dir + 'results\\'
    direct = ['merged', 'cleaned_nrem', 'cleaned_n2', 'cleaned_sws',
              'co_occur_nrem', 'co_occur_n2', 'co_occur_sws',
              'aac_summary_nrem', 'aac_summary_n2', 'aac_summary_sws',
              'outlier_check']
    spdl_dir = ['slow', 'fast']
    neuro_dir = ['slowwave', 'spindle']
    
    print('------------------------------------------\n')
    print('Creating results folders for {}....\n'.format(ppt_id))
    
    for d in direct: 
        ppt_dir = outpath + "{}\\{}\\".format(d, ppt_id)
        if d == 'merged' or d=='outlier_check' or d.startswith('cleaned'):
            try:
                for n_d in neuro_dir:
                    os.makedirs(ppt_dir + n_d)
            except FileExistsError:
                pass
                
            try:
                for s_d in spdl_dir:
                    os.mkdir(ppt_dir + 'spindle\\{}'.format(s_d))
            except FileExistsError:
                pass
        else:
            try:
                for s_d in spdl_dir:
                    os.makedirs(ppt_dir + s_d)
            except FileExistsError:
                pass
                
    print("\nFolders successfully created.\n")
    
# merge event time files with main event files
def merge_files(filetype, ppt_id, channels, root_dir, spindle_type=None):
    if filetype == 'slowwave':
        inpath = root_dir + 'source\\{}\\{}\\'.format(ppt_id, filetype)
        outpath = root_dir + 'results\\merged\\{}\\{}\\'.format(ppt_id, filetype)
        print('------------------------------------------\n')
        print("Adding slowwave trough times to event files for {}: channels {} to {}...\n".format(ppt_id, channels[0], channels[-1]))
        
    elif filetype == 'spindle':
        inpath = root_dir + 'source\\{}\\{}\\{}\\'.format(ppt_id, filetype, spindle_type)
        outpath = root_dir + 'results\\merged\\{}\\{}\\{}\\'.format(ppt_id, filetype, spindle_type)
        print('------------------------------------------\n')
        print("Adding {} spindle peak times to event files for {}: channels {} to {}...\n".format(spindle_type, ppt_id, channels[0], channels[-1]))

    # if input file exists, merge. else, log error
    errpath = root_dir + 'errors\\'
    errfile = errpath + '{}_error_log'.format(ppt_id)
    field_names = ['channel', 'filetype', 'spindle_type', 'step', 'desc', 'folder']
    
    if not os.path.exists(errfile+'.csv'):
        with open(errfile+'.csv', 'w', newline='') as f:
            dict_writer = csv.DictWriter(f, fieldnames=field_names)
            dict_writer.writeheader()
    
    for ch in channels:
        # defining input and output files
        if filetype == 'slowwave':
            infile = inpath + '{}_{}_{}'.format(ppt_id, ch, filetype)
            infile_times = inpath + '{}_{}_{}_detsw_Staresina2015_0.50-1.25Hz'.format(ppt_id, ch, filetype)
            outfile = outpath + '{}_{}_{}_troughtime'.format(ppt_id, ch, filetype)
            
        elif filetype == 'spindle':
            infile = inpath + '{}_{}_{}'.format(ppt_id, ch, filetype)
            infile_times = inpath + 'arina_{}_{}_{}_times'.format(ppt_id, ch, filetype)
            outfile = outpath + '{}_{}_{}_peaktime'.format(ppt_id, ch, filetype)
        

        if not os.path.exists(infile_times+'.csv'):
            if filetype == 'slowwave':
                print("**\nERROR: No slowwave trough time file in source folder for channel {}. Skipping...**".format(ch))
            elif filetype == 'spindle':
                print("**\nERROR: No spindle peak time file in source folder for channel {}. Skipping**".format(ch))
                
            err_desc = "No event time (peak/trough) file found"

            if filetype=='slowwave':
                err = {'channel':ch, 'filetype':filetype, 'spindle_type':'n/a', 'step':'merge', 'desc':err_desc, 'folder':inpath}
            elif filetype=='spindle':
                err = {'channel':ch, 'filetype':filetype, 'spindle_type':spindle_type, 'step':'merge', 'desc':err_desc, 'folder':inpath}
        
            with open(errfile+'.csv', 'a+', newline='') as f:
                dict_writer = csv.DictWriter(f, fieldnames=field_names)
                dict_writer.writerow(err)


        try:
            with open(infile+'.csv', 'r') as f:
                reader = csv.reader(f)
                data = list(reader)
                data_headers = data[3]
                data_content = data[8:]
                
        except FileNotFoundError:
            print("**\nERROR: No file found in source folder for channel {}**".format(ch))
            err_desc = "No event file found"
            
            if filetype=='slowwave':
                err = {'channel':ch, 'filetype':filetype, 'spindle_type':'n/a', 'step':'merge', 'desc':err_desc, 'folder':inpath}
            elif filetype=='spindle':
                err = {'channel':ch, 'filetype':filetype, 'spindle_type':spindle_type, 'step':'merge', 'desc':err_desc, 'folder':inpath}

            field_names = list(err.keys())
            
            with open(errfile+'.csv', 'a+', newline='') as f:
                dict_writer = csv.DictWriter(f, fieldnames=field_names)
                dict_writer.writerow(err)

        else:
            df_infile = pd.DataFrame(data_content, columns=data_headers)
            df_infile_times = pd.read_csv(infile_times+'.csv')
            
            if filetype=='slowwave':
                df_infile['Trough time (s)'] = df_infile_times['trough_time']
            elif filetype=='spindle':
                df_infile['Peak time (s)'] = df_infile_times['peak_time']
                
            df_infile.to_csv(outfile+'.csv', index=False)
            print("\tFiles successfully merged for {}".format(ch))
            
    print("\nAll files successfully merged.\n")
    
# detect and remove outliers 
def remove_outliers_physio(filetype, ppt_id, channels, root_dir, spindle_type=None):
    print('------------------------------------------\n')
    print("Removing {} physiological outliers for {}...\n".format(filetype, ppt_id))
    
    def check_outliers_physio(df):
        """Criteria for physiological outlier: 
            - Slowwave (a): peak to peak amplitude over 1000 uV
            - Slowwave (b): trough amplitude greater than or equal to -15 uV AND peak amplitude less than or equal to 10 uV
            - Spindle: peak to peak amplitude over 200 uV"""
            
        if filetype=='slowwave':
            if df['p2p_amplitude_uV'] > 1000:
                df['outlier_physio'] = True
            elif (df['trough_amplitude_uV']>=-15) & (df['peak_amplitude_uV']<=10):
                df['outlier_physio'] = True
            else:
                df['outlier_physio'] = False
     
        elif filetype=='spindle':
            if df['p2p_amplitude_uV'] > 200: 
                df['outlier_physio'] = True
            else:
                df['outlier_physio'] = False
            
        return df
        
    
    if filetype=='slowwave':
        inpath = root_dir + 'results\\merged\\{}\\slowwave\\'.format(ppt_id)
        outpath = root_dir + 'results\\outlier_check\\{}\\slowwave\\'.format(ppt_id)
    elif filetype=='spindle':
        inpath = root_dir + 'results\\merged\\{}\\spindle\\{}\\'.format(ppt_id, spindle_type)
        outpath = root_dir + 'results\\outlier_check\\{}\\spindle\\{}\\'.format(ppt_id, spindle_type)
    
    ch_dict = dict()

    for ch in channels:
        # read in input file from 'merge' sub-directory
        if filetype == 'slowwave':
            infile = inpath + '{}_{}_slowwave_troughtime'.format(ppt_id, ch)
        elif filetype == 'spindle':
            infile = inpath + '{}_{}_spindle_peaktime'.format(ppt_id, ch)
        
        infile_df = pd.read_csv(infile+'.csv')
        
        # keep in only cols needed for physiological outlier check. cols renamed for clarity
        outlier_df = infile_df.iloc[:, [0, 9, 10, 11]].copy()
        outlier_df.rename(columns={'Segment index':'segment_index', 'Min. amplitude (uV)':'trough_amplitude_uV', 'Max. amplitude (uV)':'peak_amplitude_uV', 'Peak-to-peak amplitude (uV)':'p2p_amplitude_uV'}, inplace=True)

        # apply physiological outlier check criteria
        df_outliers = outlier_df.apply(check_outliers_physio, axis=1).copy()
        df_outliers = df_outliers[df_outliers.outlier_physio==1]
        df_outliers['channel'] = ch
        
        # rearrange columns
        if filetype == 'slowwave':
            df_outliers = df_outliers[['channel', 'segment_index', 'trough_amplitude_uV', 'peak_amplitude_uV', 'p2p_amplitude_uV']].copy()
        elif filetype=='spindle':
            df_outliers = df_outliers[['channel', 'segment_index', 'p2p_amplitude_uV']].copy()

        # remove all physio outliers from infile
        outlier_segments = df_outliers['segment_index'].to_list()
        
        df_outliers['total_events_original'] = infile_df.shape[0]
        df_outliers['physio_outliers_abs'] = len(outlier_segments)
        df_outliers['physio_outliers_rel'] = round((len(outlier_segments)/infile_df.shape[0])*100,2)
        df_outliers['total_events_nophysio_outliers'] = df_outliers['total_events_original'] - df_outliers['physio_outliers_abs']
        
        df_outliers = df_outliers[['channel', 'total_events_original', 'physio_outliers_abs', 'physio_outliers_rel', 'total_events_nophysio_outliers']]

        
        infile_cleaned_df = infile_df[~(infile_df['Segment index'].isin(outlier_segments))].copy()

        if filetype=='slowwave':
            infile_cleaned_df.to_csv(inpath+'{}_{}_slowwave_troughtime_nophys.csv'.format(ppt_id, ch), index=False)
        elif filetype=='spindle':
            infile_cleaned_df.to_csv(inpath+'{}_{}_spindle_peaktime_nophys.csv'.format(ppt_id, ch), index=False)
            
        print('\tPhysio outliers removed for {}'.format(ch))


        # add df to channel dictionary to eventually create master physio outlier file
        if df_outliers.empty == True:
            no_outliers = {'channel':ch, 'total_events_original':infile_df.shape[0], 'physio_outliers_abs':0, 'physio_outliers_rel':0, 'total_events_nophysio_outliers':infile_df.shape[0]}
            df_no_outliers= pd.DataFrame(data=no_outliers, index=[0])
            
            ch_dict[ch] = df_no_outliers
            
        elif df_outliers.empty == False: 
            ch_dict[ch] = df_outliers
        
        # write masterfile for each ppt     
        outfile_df_dict = {'channel':np.NaN, 'total_events_original':np.NaN, 'physio_outliers_abs':np.NaN, 'physio_outliers_rel':np.NaN, 'total_events_nophysio_outliers':np.NaN}
        outfile_df = pd.DataFrame(data=outfile_df_dict, index=[0])
        
        for ch, df in ch_dict.items():
            outfile_df = pd.concat( [ outfile_df, ch_dict[ch] ] )
            
        outfile_df.drop_duplicates(subset=['channel'], inplace=True)
        
    outfile_df.dropna(how='all', inplace=True)

    
    if filetype=='slowwave':
        outfile_df.to_csv(outpath+'{}_slowwave_outliers_physio.csv'.format(ppt_id), index=False)
    elif filetype=='spindle':
        outfile_df.to_csv(outpath+'{}_spindle_outliers_physio.csv'.format(ppt_id), index=False)


    return outfile_df


def check_outliers_stats(filetype, ppt_id, channels, root_dir, spindle_type=None):
    if filetype=='slowwave':
        inpath = root_dir + 'results\\merged\\{}\\slowwave\\'.format(ppt_id)
        #outpath = root_dir + 'results\\outlier_check\\{}\\slowwave\\'.format(ppt_id)
    elif filetype=='spindle':
        inpath = root_dir + 'results\\merged\\{}\\spindle\\{}\\'.format(ppt_id, spindle_type)
        #outpath = root_dir + 'results\\outlier_check\\{}\\spindle\\{}\\'.format(ppt_id, spindle_type)
        
    ch_dict = dict()
    
    for ch in channels: 
        if filetype=='slowwave':
            infile_clean = inpath + '{}_{}_slowwave_troughtime_nophys'.format(ppt_id, ch)
        elif filetype=='spindle':
            infile_clean = inpath + '{}_{}_spindle_peaktime_nophys'.format(ppt_id, ch)
            
        infile_clean_df = pd.read_csv(infile_clean+'.csv')
        
        # keep only the cols relevant for stats outlier check
        stats_outlier_df = infile_clean_df.iloc[:, [0,11]].copy()
        stats_outlier_df.rename(columns={'Segment index':'segment_index', 'Peak-to-peak amplitude (uV)':'p2p_amplitude_uV'}, inplace=True)
        
 
        # calculate z-scores for the amplitude columns
        stats_outlier_z = stats_outlier_df.iloc[:, 1:].apply(stats.zscore).copy()
        if len(stats_outlier_z) != 0:
            for col in ['p2p_amplitude_uV']:
                stats_outlier_df['z_{}'.format(col)] = stats_outlier_z[col]
    
            # flagging statstical outliers 
            for col in ['p2p_amplitude_uV']:
                stats_outlier_df['outlier_{}'.format(col)] = stats_outlier_df['z_{}'.format(col)].apply(lambda x: True if abs(x)>=2.5 else False)
        
            for col in ['p2p_amplitude_uV']:
                try:
                    stats_outlier_df['outlier_{}_perc'.format(col)] = round((stats_outlier_df['outlier_{}'.format(col)].value_counts()[True]/infile_clean_df.shape[0])*100,2)
                except KeyError:
                    stats_outlier_df['outlier_{}_perc'.format(col)] = 0
    
                
        # rearranging columns
        stats_outlier_df['channel'] = ch
        
        # isolating the outliers 

        outfile_df_ch = stats_outlier_df[stats_outlier_df.outlier_p2p_amplitude_uV==True].copy()
        
        # calculating proportion of file outliers make up 
        outfile_df_ch['total_outliers_perc'] = round((outfile_df_ch.shape[0]/infile_clean_df.shape[0])*100,2)
            
        # add to channel dictionary 
        ch_dict[ch] = outfile_df_ch
    
    # write masterfile for each ppt
    master_ch = list(ch_dict.keys())[0]
    outfile_df = ch_dict[master_ch]
    del ch_dict[master_ch]
    
    for ch,df in ch_dict.items():
        outfile_df = pd.concat([outfile_df, ch_dict[ch]])
    
    outfile_df['remove_yn'] = True
    
    return outfile_df

def remove_outliers_stats(annot_df, filetype, ppt_id, channels, root_dir, spindle_type=None):
    print('\nRemoving {} statstical outliers for {}...\n'.format(filetype, ppt_id))
    
    ch_dict = dict()
    
    inpath = root_dir + 'results\\merged\\{}\\'.format(ppt_id)
    outpath_cleaned = root_dir + 'results\\cleaned_nrem\\{}\\'.format(ppt_id)
    outpath_outliercheck = root_dir + 'results\\outlier_check\\{}\\'.format(ppt_id)
    
    # read event files (no phys outliers)
    for ch in channels:
        if filetype=='slowwave':
            infile = inpath + 'slowwave\\{}_{}_slowwave_troughtime_nophys'.format(ppt_id, ch)
        elif filetype=='spindle':
            infile = inpath + 'spindle\\{}\\{}_{}_spindle_peaktime_nophys'.format(spindle_type, ppt_id, ch)
        df_nophys = pd.read_csv(infile+'.csv')
        
        df_stats_outliers = annot_df[(annot_df.channel==ch) & (annot_df.remove_yn==True)].copy()
        outlier_segments = df_stats_outliers['segment_index'].to_list()
        
        total_events_nophysio = df_nophys.shape[0]
        
        # summarise outlier numbers 
        df_stats_outliers['total_events_nophysio_outliers'] = total_events_nophysio
        df_stats_outliers['stats_outliers_abs'] = len(outlier_segments)
        df_stats_outliers['stats_outliers_rel'] = round(len(outlier_segments)/total_events_nophysio*100,2)
        df_stats_outliers['total_events_nostats_outliers'] = df_stats_outliers['total_events_nophysio_outliers'] - df_stats_outliers['stats_outliers_abs']
        
        df_stats_outliers = df_stats_outliers[ ['channel', 'total_events_nophysio_outliers', 'stats_outliers_abs','stats_outliers_rel', 'total_events_nostats_outliers'] ]
        
        if df_stats_outliers.empty == True:
            no_outliers = {'channel':ch, 'total_events_nophysio_outliers':total_events_nophysio, 'stats_outliers_abs':0,
                           'stats_outliers_rel':0, 'total_events_nostats_outliers':total_events_nophysio}
            df_no_stats_outliers = pd.DataFrame(data=no_outliers, index=[0])
            
            ch_dict[ch] = df_no_stats_outliers
            
        elif df_stats_outliers.empty == False:
            ch_dict[ch] = df_stats_outliers
            
        # export cleaned event files    
        infile_cleaned_full_df = df_nophys[~(df_nophys['Segment index'].isin(outlier_segments))].copy()
        
        if filetype=='slowwave':
            infile_cleaned_full_df.to_csv(outpath_cleaned+'slowwave\\{}_{}_slowwave_troughtime_clean_nrem.csv'.format(ppt_id, ch), index=False)
        elif filetype=='spindle':
            infile_cleaned_full_df.to_csv(outpath_cleaned+'spindle\\{}\\{}_{}_spindle_peaktime_clean_nrem.csv'.format(spindle_type, ppt_id, ch), index=False)
        
        print('\tStats outliers removed for {}'.format(ch))
        
        # write outlier check masterfile for each ppt
        outfile_df_dict = {'channel':np.NaN, 'total_events_nophysio_outliers':np.NaN}
        outfile_df = pd.DataFrame(data=outfile_df_dict, index=[0])
        
        for ch, df in ch_dict.items():
            outfile_df = pd.concat( [ outfile_df, ch_dict[ch] ] )
            
        outfile_df.drop_duplicates(subset=['channel'], inplace=True)
        
    outfile_df.dropna(how='all', inplace=True)

    if filetype=='slowwave':
        outfile_df.to_csv(outpath_outliercheck + 'slowwave\\{}_slowwave_outliers_stats.csv'.format(ppt_id), index=False)
    elif filetype=='spindle':
        outfile_df.to_csv(outpath_outliercheck + 'spindle\\{}\\{}_spindle_outliers_stats.csv'.format(spindle_type, ppt_id), index=False)

    
    return outfile_df


def output_summary(df_nophys, df_nostats, filetype, ppt_id, channels, root_dir, spindle_type=None):
    outpath = root_dir + 'results\\outlier_check\\{}\\'.format(ppt_id)
    
    df_nophys.reset_index(drop=True, inplace=True)
    df_nostats.reset_index(drop=True, inplace=True)
    df_nostats = df_nostats.iloc[:,2:].copy()
    
    outfile = pd.concat([df_nophys, df_nostats], axis=1, join='inner')
    
    if filetype=='slowwave':
        outfile.to_csv(outpath + 'slowwave\\{}_outliers_slowwave_summary_nrem.csv'.format(ppt_id), index=False)
    elif filetype=='spindle':
        outfile.to_csv(outpath + 'spindle\\{}\\{}_outliers_spindle_summary_nrem.csv'.format(spindle_type, ppt_id), index=False)
        
    print("\nAll outliers successfully removed.\n")
    

# load hypnogram and staging file
def load_hypno(ppt_id, root_dir):
    inpath_hypno = root_dir + 'hypno\\'
    
    # import staging file and get first epoch time
    infile_scored = inpath_hypno + '{}_ScoredEvents'.format(ppt_id)
    scored_df = pd.read_table(infile_scored+'.txt', header=None, delimiter=',')
    scored_df[0] = pd.to_timedelta(scored_df[0])
    epoch_time = scored_df.iloc[0,0]

    
    
    infile_hypno = inpath_hypno + '{}_staging'.format(ppt_id)
    hypno_df = pd.read_table(infile_hypno+'.txt', header=None, delimiter=',')
    hypno_df.rename(columns={0:'stage'}, inplace=True)
    
    

    # create list of times 
    time_30s = timedelta(seconds=30)
    start_times_list = list()
    
    for i in range(hypno_df.shape[0]):
        start_times_list.append(epoch_time)
        epoch_time += time_30s
    
    
    # update hypno_df with relevant cols
    hypno_df['epoch'] = hypno_df.index + 1
    hypno_df['epoch_start'] = start_times_list
    hypno_df['epoch_end'] = hypno_df['epoch_start'] + timedelta(seconds=30)
    
    return hypno_df

# section events to N2 or SWS
def split_events(hypno_df, filetype, ppt_id, channels, root_dir, spindle_type=None):
    inpath_events = root_dir + 'results\\cleaned_nrem\\{}\\'.format(ppt_id)

    outpath_n2 = root_dir + 'results\\cleaned_n2\\{}\\'.format(ppt_id)
    outpath_sws = root_dir + 'results\\cleaned_sws\\{}\\'.format(ppt_id)
    

    hypno_df['epoch_start'] = pd.to_timedelta(hypno_df['epoch_start'])
    hypno_df['epoch_end'] = pd.to_timedelta(hypno_df['epoch_end'])
    
    first_epoch = hypno_df.iloc[0,2]
    
    # creating dictionaries 
    hypno_start = hypno_df['epoch_start'].to_list()
    hypno_end = hypno_df['epoch_end'].to_list()
    
    times_both = list()
    
    for i in range(len(hypno_start)):
        times_both.append( [hypno_start[i],hypno_end[i]] )  
        
    epoch_stage_dict = dict( zip(hypno_df['epoch'], hypno_df['stage']) )
    epoch_times_dict = dict( zip(hypno_df['epoch'], times_both) )
    
    
    def convert(seconds):
        return time.strftime("%H:%M:%S", time.gmtime(seconds))
    
    def check_epochs(df_events):
        for epoch,times in epoch_times_dict.items():
            epoch_start = times[0]
            epoch_end = times[1]
        
            if (epoch_start<=df_events['start_time_final']) & (df_events['start_time_final']<=epoch_end):
                df_events['epoch'] = epoch
            
        return df_events

    def check_stage(df_events):
        for epoch,stage in epoch_stage_dict.items():
            if df_events['epoch'] == epoch:
                df_events['stage'] = stage
            
        return df_events
    
    #########################
    print('------------------------------------------\n')
    print("Splitting detected events by stage...\n")
    
    for ch in channels:
        
        ### EVENT FILE 
    
        # import event file
        if filetype=='slowwave':
            events = inpath_events + 'slowwave\\{}_{}_slowwave_troughtime_clean_nrem'.format(ppt_id, ch)
        elif filetype=='spindle':
            events = inpath_events + 'spindle\\{}\\{}_{}_spindle_peaktime_clean_nrem'.format(spindle_type, ppt_id, ch)
        
        df_events = pd.read_csv(events+'.csv')
        
        
        # convert time to timedelta
        df_events['start_time_hours'] = df_events['Start time'].apply(convert)
        df_events['start_time_hours'] = pd.to_timedelta(df_events['start_time_hours'])
        
        
        # create timedelta col to compare to staging file 
        df_events['start_time_final'] = df_events['start_time_hours'] + first_epoch
    
        # apply functions
        df_events_export = df_events.apply(check_epochs, axis=1).copy()
        df_events_export = df_events_export.apply(check_stage, axis=1).copy()
        


        # ### EXPORT EPOCHS FILE
        
        # keep only relevant cols
        cols_of_interest = ['Segment index', 'Start time', 'End time', 'start_time_final', 'epoch', 'stage']
        df_events_export = df_events_export.loc[:,cols_of_interest].copy()
        
        
        ### EXPORT EVENT FILES PER STAGE
        
        # remove irrelevant cols
        df_events.drop(columns=['start_time_hours', 'start_time_final'], inplace=True)
        
        df_epochs_sws = df_events_export[df_events_export.stage=='3'].copy()
        sws_segments_list = df_epochs_sws['Segment index'].to_list()
        
        # splitting df based on N2/SWS epochs 
        df_events_sws = df_events[df_events['Segment index'].isin(sws_segments_list)].copy()
        df_events_n2 = df_events[~(df_events['Segment index'].isin(sws_segments_list))].copy()
        
        print("\tEvents split for {}".format(ch))
        
        # export file
        if filetype=='slowwave':
            export_file_n2 = outpath_n2 + 'slowwave\\{}_{}_slowwave_troughtime_clean_n2'.format(ppt_id, ch)
            export_file_sws = outpath_sws + 'slowwave\\{}_{}_slowwave_troughtime_clean_sws'.format(ppt_id, ch)
            
        elif filetype=='spindle':
            export_file_n2 = outpath_n2 + 'spindle\\{}\\{}_{}_spindle_peaktime_clean_n2'.format(spindle_type, ppt_id, ch)
            export_file_sws = outpath_sws + 'spindle\\{}\\{}_{}_spindle_peaktime_clean_sws'.format(spindle_type, ppt_id, ch)
            
        df_events_sws.to_csv(export_file_sws+'.csv', index=False)
        df_events_n2.to_csv(export_file_n2+'.csv', index=False)

    print("\nAll events successfully split.")
    
# calculate total stage duration for densities
def calc_duration(hypno_df, stage):
    epochs_cat = hypno_df['stage'].value_counts()
    
    if stage=='N2' or stage=='NREM2':
        epochs_total = epochs_cat['2']
    elif stage=='SWS' or stage=='NREM3':
        epochs_total = epochs_cat['3']
    elif stage=='NREM' or stage=='NREM2NREM3':
        epochs_total = epochs_cat['2'] + epochs_cat['3']
        
    stage_dur = (epochs_total*30)/60
    return stage_dur

def check_summfile(spindle_type, ppt_id, root_dir, channels, stage):  
    if stage=='NREM2':
        stage='n2'
    elif stage=='NREM3':
        stage='sws'
    elif stage=='NREM2NREM3':
        stage='nrem'
    else:
        stage=stage.lower()
        
    inpath = root_dir + 'results\\aac_summary_{}\\{}\\{}\\'.format(stage, ppt_id, spindle_type)
    infile = inpath + '{}_aac_summary_{}'.format(ppt_id, stage)
    
    if not os.path.exists(infile+'.csv'):
        with open(infile+'.csv', 'a+', newline='') as f:
            field_names = ['ppt_id', 'channel',
                           'total_slowwave_count', 'total_spindle_count', 'coupled_spindle_count', 'uncoupled_spindle_count',
                           'total_slowwave_density', 'total_spindle_density', 'coupled_spindle_density', 'uncoupled_spindle_density']
            writer = csv.writer(f)
            writer.writerow(field_names)
        

# calculate AAC
def aac(stage_duration, spindle_type, ppt_id, root_dir, channels, stage):
    print('------------------------------------------\n')
    print("Calculating co-occurrence for {}...\n".format(ppt_id))
    
    if stage=='NREM2':
        stage='n2'
    elif stage=='NREM3':
        stage='sws'
    elif stage=='NREM2NREM3':
        stage='nrem'
    else:
        stage=stage.lower()
    
    # define filepaths
    inpath = root_dir + 'results\\cleaned_{}\\{}\\'.format(stage, ppt_id)
    outpath_occur = root_dir + 'results\\co_occur_{}\\{}\\{}\\'.format(stage, ppt_id, spindle_type)
    outpath_summ = root_dir + 'results\\aac_summary_{}\\{}\\{}\\'.format(stage, ppt_id, spindle_type)
    
    errpath = root_dir + 'errors\\'
    errfile = errpath + '{}_error_log'.format(ppt_id)
    field_names = ['channel', 'filetype', 'spindle_type', 'step', 'desc', 'folder']
    
    if not os.path.exists(errfile+'.csv'):
        with open(errfile+'.csv', 'w', newline='') as f:
            dict_writer = csv.DictWriter(f, fieldnames=field_names)
            dict_writer.writeheader()

    # columns to rename
    sw_cols = ['Segment index', 'Start time', 'End time', 'Trough time (s)']
    sw_cols_renamed = {'Segment index':'slowwaves_segment_index',
                         'Start time':'slowwaves_start_time',
                         'End time':'slowwaves_end_time',
                         'Trough time (s)':'slowwaves_trough_time'}

    spd_cols = ['Segment index', 'Start time', 'End time', 'Peak time (s)']
    spd_cols_renamed = {'Segment index':'spindles_segment_index',
                            'Start time':'spindles_start_time',
                            'End time':'spindles_end_time',
                            'Peak time (s)':'spindles_peak_time'}
    
    
    ## functions to check co-occurrence
    def get_cols(df, cols_list, rename_dict):
        df_iso = df[cols_list].copy()
        for i in cols_list: 
            df_iso[i] = df_iso[i].astype('float')
        
        df_iso.rename(columns=rename_dict, inplace=True)
        
        return df_iso
    
    def check_occur(spd_df):
        for (seg, e_time) in sw_dict.items(): 
            
            if spd_df['co_occur'] != True: 
                before_trough = e_time-1.2
                after_trough = e_time+1.2
                
                mask_1 = ( before_trough<=spd_df['spindles_start_time']<=after_trough ) & ( before_trough<=spd_df['spindles_peak_time']<=after_trough )
                mask_2 = ( before_trough<=spd_df['spindles_peak_time']<=after_trough ) & ( before_trough<=spd_df['spindles_end_time']<=after_trough )
                
                if (mask_1) | (mask_2):
                    spd_df['co_occur'] = True
                    spd_df['co_occur_slowwave_segment_index'] = seg
                    spd_df['co_occur_slowwave_trough_time'] = e_time
        
        return spd_df
    
    
    # iterate over channels 
    for ch in channels:
        
        infile_sw = inpath + 'slowwave\\{}_{}_slowwave_troughtime_clean_{}'.format(ppt_id, ch, stage)
        infile_spd = inpath + 'spindle\\{}\\{}_{}_spindle_peaktime_clean_{}'.format(spindle_type, ppt_id, ch, stage)
        outfile_occur = outpath_occur + '{}_{}_co_occur_{}'.format(ppt_id, ch, stage)
        outfile_summ = outpath_summ + '{}_aac_summary_{}'.format(ppt_id, stage)
        
        try:
            df_sw = pd.read_csv(infile_sw+'.csv')
        except FileNotFoundError:
            print("**\nERROR: No slowwave file found in cleaned_{} folder for channel {}**".format(stage, ch))
            err_desc = "No cleaned slowwave event file found"
            err = {'channel':ch, 'filetype':'slowwave', 'spindle_type':'n\a', 'step':'co_occur', 'desc':err_desc, 'folder':inpath}
            field_names = list(err.keys())
            
            with open(errfile+'.csv', 'a+', newline='') as f:
                dict_writer = csv.DictWriter(f, fieldnames=field_names)
                dict_writer.writerow(err)
        else:
            df_sw_iso = get_cols(df_sw, sw_cols, sw_cols_renamed)  
            sw_events = len(df_sw_iso['slowwaves_segment_index'])
            sw_dict = dict(zip(df_sw_iso.slowwaves_segment_index, df_sw_iso.slowwaves_trough_time))

        
        try:
            df_spd = pd.read_csv(infile_spd+'.csv')
        except FileNotFoundError:
            print("**\nERROR: No spindle file found in cleaned_{} folder for channel {}**".format(stage, ch))
            err_desc = "No cleaned spindle event file found"
            err = {'channel':ch, 'filetype':'spindle', 'spindle_type':spindle_type, 'step':'co_occur', 'desc':err_desc, 'folder':inpath}
            field_names = list(err.keys())
            
            with open(errfile+'.csv', 'a+', newline='') as f:
                dict_writer = csv.DictWriter(f, fieldnames=field_names)
                dict_writer.writerow(err)
                
        else:
            df_spd_iso = get_cols(df_spd, spd_cols, spd_cols_renamed)            
            df_spd_iso['co_occur'] = False 
            df_spd_iso['co_occur_slowwave_segment_index'] = np.NaN
            df_spd_iso['co_occur_slowwave_trough_time'] = np.NaN
            

        
        df_co_occur = df_spd_iso.apply(check_occur, axis=1).copy()
            
        df_co_occur = df_co_occur.iloc[:, [0,3,4,5,6]]
        df_co_occur.to_csv(outfile_occur+'.csv', index=False)
        
        
        
            
        ## write to summary file
        absolute_results = df_co_occur.co_occur.value_counts()
        #relative_results = df_co_occur.co_occur.value_counts(normalize=True)
            
        with open(outfile_summ+'.csv', 'a+', newline='') as f:
            field_names = ['ppt_id', 'channel',
                           'total_slowwave_count', 'total_spindle_count', 'coupled_spindle_count', 'uncoupled_spindle_count',
                           'total_slowwave_density', 'total_spindle_density', 'coupled_spindle_density', 'uncoupled_spindle_density']
            
            summ = {'ppt_id':ppt_id, 'channel':ch,
                    'total_slowwave_count':sw_events, 'total_spindle_count':absolute_results.sum(), 'coupled_spindle_count':absolute_results[True], 'uncoupled_spindle_count':absolute_results[False],
                    'total_slowwave_density':sw_events/stage_duration, 'total_spindle_density':absolute_results.sum()/stage_duration, 'coupled_spindle_density':absolute_results[True]/stage_duration, 'uncoupled_spindle_density':absolute_results[False]/stage_duration}
            
            dict_writer = csv.DictWriter(f, fieldnames=field_names)
            dict_writer.writerow(summ)
            
        print("\tCo-occurrence calculated for {}".format(ch))
        
    print("\nCo-occurrence successfully calculated for all files.\n")


    

