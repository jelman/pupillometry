# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 12:00:45 2016

@author: jelman

This script takes as input a data file that has been parsed with 
the parse_oddball_LCIP.py script. It performs interpolation, artifact 
removal, and averaging. It calculates pupil dilation for targets and 
non-targets, then calculates a signal to noise ratio for use in 
further analysis. 
"""

import itertools
import os
import re
import sys
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0,'/home/jelman/netshare/K/code/Pupillometry/LCIP')
import parse_oddball_LCIP as parse_ao

import matlab_wrapper



#def create_ao_template(df):
#    trialnums = get_trial_nums(df)
#    templateDF = pd.DataFrame({"ID": df['Subject ID'],
#                           "Task (AO)": ['AO'+str(x) for x in trialnums],
#                           "Time": df['Time'],
#                           "Trial #": trialnums,
#                           "Min": np.nan,
#                           "Max": np.nan,
#                           "Pupil Profile": df['Pupil Profile']},
#                           columns=['ID','Task (AO)','Time','Trial #',
#                                    'Min','Max','Pupil Profile'])
#    pprofileDF = templateDF.pop('Pupil Profile').str.split('\t', expand=True)
#    pprofileDF = pprofileDF.rename(columns=lambda x: str(x+1) + 'st data pt')
#    templateDF = pd.concat([templateDF,pprofileDF], axis=1)
#    templateDF['Min'] = templateDF.ix[:,'1st data pt':].min(axis=1, skipna=True)
#    templateDF['Max'] = templateDF.ix[:,'1st data pt':].max(axis=1, skipna=True)
#    return templateDF



def get_pupil_profile(df):
    timeprofile = pd.Series(df['Time Profile'].str.split('\t'))
    time_dt = pd.to_datetime(timeprofile.values[0], unit='s')
    time_idx = pd.DatetimeIndex(time_dt)
    pupilprofile = df['Pupil Profile'].str.split('\t').values[0]
    pupilprofile = pd.DataFrame({'PupilProfile':pupilprofile},
                                 index=time_idx)
    pupilprofile.PupilProfile = pupilprofile.PupilProfile.astype('float')
    return pupilprofile
    

def upsample_interp(pupilprofile):
    upsampled = pupilprofile.resample('67ms').mean()
    interpolated = upsampled.interpolate()
    return interpolated


def stublinks(ts):
    matlab = matlab_wrapper.MatlabSession()
    matlab.eval("addpath(genpath('/home/jelman/netshare/K/code/Pupillometry/LCIP'))")
    res = matlab.workspace.stublinks_new(ts.PupilProfile.values,0,0.,[],0.1,15.)
    origdata = res['RescaleData'].tolist()
    blinktimes = pd.Series(res['BlinkTimes'].tolist(), 
                           index=ts.index)
    noblinkdata = pd.Series(res['NoBlinks'].tolist(), 
                            index=ts.index)
    noblinkdata_unsmooth = pd.Series(res['NoBlinksUnsmoothed'].tolist(),
                                     index=ts.index)
    return origdata, blinktimes, noblinkdata, noblinkdata_unsmooth
    

def clean_sess(sessdf, qc_dir, subid, session):
    pupilprofile = get_pupil_profile(sessdf)
    interpolated = upsample_interp(pupilprofile)
    origdata, blinktimes, noblinkdata, noblinkdata_unsmooth = stublinks(interpolated)
    fname = os.path.join(qc_dir, ''.join([subid,'_',str(session),'_stublinks.png']))
    plt.plot(range(len(origdata)), origdata, sns.xkcd_rgb["pale red"], 
             range(len(noblinkdata)), noblinkdata, sns.xkcd_rgb["denim blue"], 
             blinktimes, sns.xkcd_rgb["amber"], lw=1)
    plt.title(''.join(['Subject: ',subid,', Session:',str(session)]), fontsize=20)
    plt.savefig(fname)
    plt.close()
    return noblinkdata, blinktimes
    
    
def clean_all(parsed_df, qc_dir): 
    grp_df = parsed_df.groupby(['Subject ID','Session'])
    noblinkdata_dict = {}
    blinktimes_dict = {}
    for (subid, session), sessdf in grp_df:
        sess_noblinkdata, sess_blinktimes = clean_sess(sessdf, qc_dir, subid, session)
        noblinkdata_dict.update({(subid,session) : sess_noblinkdata})
        blinktimes_dict.update({(subid,session) : sess_blinktimes})
    noblinkdata = pd.DataFrame.from_dict(noblinkdata_dict)
    blinktimes = pd.DataFrame.from_dict(blinktimes_dict)
    noblinkdata.columns = zip(noblinkdata.columns.get_level_values(0), noblinkdata.columns.get_level_values(1))
    blinktimes.columns = zip(blinktimes.columns.get_level_values(0), blinktimes.columns.get_level_values(1))
    return noblinkdata, blinktimes


def get_events(pupilprofile, eprime):
    trg_onsets = eprime.index[eprime['Target_Detected']==1]
    std_onsets = eprime.index[eprime['Correct_Rejection']==1][1:]
    
    # Get event indices
    trg_event_idx = pupilprofile.index.searchsorted(trg_onsets)
    std_event_idx = pupilprofile.index.searchsorted(std_onsets)
    
    trg_events = pupilprofile.index[trg_event_idx]
    std_events = pupilprofile.index[std_event_idx]
    return trg_events, std_events


def calc_max_dilation(events, pupilprofile, blinks, tpre=0.5, tpost=2.5):
    all_dilations = []
    for event in events:
        pre_event = event - pd.tseries.offsets.relativedelta(seconds=tpre)
        baseline = pupilprofile[pre_event:event].mean()
        post_event = event + pd.tseries.offsets.relativedelta(seconds=tpost)
        if blinks[event:post_event].mean() >= .5:
            print 'Blinks during >50% of trial, skipping...'
            continue
        else:
            normed_post_event = pupilprofile[event:post_event] - baseline
        tdelta = pd.tseries.offsets.relativedelta(seconds=tpost)
        max_dilation = normed_post_event[event:event+tdelta].max()
        all_dilations.append(max_dilation)
    return np.mean(all_dilations)
    
    
def calc_sess_snr(noblink_series, blinktimes, ao_eprime):
    blinktime_series = blinktimes[noblink_series.name]
    subid, sess = noblink_series.name
    eprime_sess = ao_eprime[(ao_eprime['Subject_ID']==subid) & 
                                 (ao_eprime['Session']==sess)]
    eprime_sess.index = pd.to_datetime(list(eprime_sess.Tone_Onset), unit='ms')
    trg_events, std_events = get_events(noblink_series, eprime_sess)
    trg_dilation = calc_max_dilation(trg_events, noblink_series, blinktime_series)
    std_dilation = calc_max_dilation(std_events, noblink_series, blinktime_series)
    sess_snr = trg_dilation / std_dilation
    return sess_snr
   
   
def get_blink_pct(blinktimes):
    blinkpct = blinktimes.apply(np.mean)
    blinkpct.name = 'blink_pct'
    blinkpct.index = pd.MultiIndex.from_tuples(blinkpct.index) 
    blinkpct.index.names = ['Subject_ID', 'Session']
    return blinkpct
    
    
def calc_subj_snr(noblinkdata, blinktimes, ao_eprime):
    subj_sess_snr = noblinkdata.apply(calc_sess_snr, args=(blinktimes, ao_eprime))
    subj_sess_snr.name = 'SNR'
    subj_sess_snr.index = pd.MultiIndex.from_tuples(subj_sess_snr.index)    
    return pd.DataFrame(subj_sess_snr)
    
    
def proc_oddball(pupil_fname, behav_fname, outdir):
    parsed_df =  parse_ao.parse_pupil_data(pupil_fname)
    fulltrials = parsed_df['Measurement Duration'].str.replace('sec','').astype('float') > 290
    parsed_df = parsed_df.ix[fulltrials]
    ao_eprime = pd.read_csv(behav_fname)
    sessioninfo = parsed_df[['Subject ID', 'Time', 'Session', 'Device ID', 'Eye Measured']]
    qc_dir = os.path.join(outdir,'QC') 
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    noblinkdata, blinktimes = clean_all(parsed_df, qc_dir=qc_dir)
    blinkpct = get_blink_pct(blinktimes, qc_dir)
    outfile = os.path.join(qc_dir, "Blink_Percentages.csv")
    blinkpct.to_csv(outfile, index=True, header=True)
    subj_sess_snr = calc_subj_snr(noblinkdata, blinktimes, ao_eprime)
    subj_info_snr = sessioninfo.merge(subj_sess_snr, left_on=['Subject ID','Session'], right_index=True)
    # subj_info_snr = subj_info_snr.merge(blinkpct, left_on=['Subject ID','Session'], right_index=True)
    outfile = os.path.join(outdir, 'LCIP_Oddball_SNR.csv')
    subj_info_snr.to_csv(outfile, index=False, header=True)    
    
    
#if __name__ == '__main__':
#
#'/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/raw/OddballP300_LCI_Pilot_007-1_8-29-2016.csv'
#    if len(sys.argv) == 1:
#        print 'USAGE: %s <file> <output directory>' % os.path.basename(sys.argv[0])
#        print 'Takes pupillometer data text file outputs csv files for use'
#        print 'in further analysis.'
#    else:
#        filename = sys.argv[1]
#        outdir = sys.argv[2]
#        parse_pupil_data(filename, outdir)


#########################################################
#                         TESTING                       #
#########################################################
pupil_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/raw/R_20160913_0902_20160915_1246.dat.txt'
behav_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/raw/oddball/OddballP300_LCI_Pilot_AllSubjects.csv'
outdir = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data'
qc_plots=True
# Plot timecourse
sns.set_style('ticks')
plt_pre = 0.5
plt_post = 4
plt_pre_event = event_time - pd.tseries.offsets.relativedelta(seconds=plt_pre)

plt_baseline = pupilprofile[plt_pre_event:event_time].mean()
plt_post_event = event_time + pd.tseries.offsets.relativedelta(seconds=plt_post)
plt_normed_post_event = pupilprofile[event_time:plt_post_event] - plt_baseline
plt_normed_post_event[event_time:plt_post_event].plot()


#################################
#           Misc. code          #
#################################

## Save out data
#    aoCols = ['Subject ID', 'Time', 'Device ID', 'Eye Measured', 'Record ID',
#              'Profile Normal', 'Measurement Duration',
#              'Pupil Profile']  
#    ao_outfile = os.path.join(outdir,'LCIP_AO_Data.csv')
#    save_to_csv(ao_df, ao_outfile, aoCols)

##  To get index into float
# a = noblinkdata.index - noblinkdata.index[0]
# a.components.milliseconds/1000.
# or
# a.microseconds/1000000.