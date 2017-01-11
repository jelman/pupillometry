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

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import parse_oddball_LCIP as parse_ao
from datetime import datetime
import matlab_wrapper


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
    res = matlab.workspace.stublinks_new(ts.PupilProfile.values,0,0.,[],0.5,15.)
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
             blinktimes.values, sns.xkcd_rgb["amber"], lw=1)
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


def calc_trial_dilations(events, pupilprofile, blinks, tpre=.5, tpost_start=1, tpost_end=2):
    mean_dilations = []
    max_dilations = []
    for event in events:
        pre_event = event - pd.tseries.offsets.relativedelta(seconds=tpre)
        baseline = pupilprofile[pre_event:event].mean()
        post_event_start = event + pd.tseries.offsets.relativedelta(seconds=tpost_start)
        post_event_end = event + pd.tseries.offsets.relativedelta(seconds=tpost_end)
        if blinks[event:post_event_end].mean() >= .5:
            print 'Blinks during >50% of trial, skipping...'
            continue
        else:
            normed_post_event = pupilprofile[post_event_start:post_event_end] - baseline
            mean_dil = normed_post_event.mean()
            max_dil = normed_post_event.max()
            mean_dilations.append(mean_dil)
            max_dilations.append(max_dil)
    return np.nanmean(mean_dilations), np.nanmean(max_dilations)
  
    
def calc_sess_snr(noblink_series, blinktimes, ao_eprime):
    blinktime_series = blinktimes[noblink_series.name]
    subid, sess = noblink_series.name
    eprime_sess = ao_eprime[(ao_eprime['Subject_ID']==subid) & 
                                 (ao_eprime['Session']==sess)]
    eprime_sess.index = pd.to_datetime(list(eprime_sess.Tone_Onset), unit='ms')
    trg_events, std_events = get_events(noblink_series, eprime_sess)
    trg_mean_dil, trg_max_dil = calc_trial_dilations(trg_events, noblink_series, blinktime_series)
    std_mean_dil, std_max_dil = calc_trial_dilations(std_events, noblink_series, blinktime_series)
    sess_snr_max = trg_max_dil / std_max_dil
    sess_snr_mean = trg_mean_dil / std_mean_dil 
    resultdict = dict(Trg_max = trg_max_dil, Trg_mean = trg_mean_dil,
                      Std_max = std_max_dil, Std_mean = std_mean_dil,
                      SNR_max=sess_snr_max, SNR_mean=sess_snr_mean)
    return pd.Series(resultdict)
   
   
def get_blink_pct(blinktimes):
    blinkpct = blinktimes.apply(np.mean)
    blinkpct.name = 'blink_pct'
    blinkpct.index = pd.MultiIndex.from_tuples(blinkpct.index) 
    blinkpct.index.names = ['Subject ID', 'Session']
    return blinkpct
    
    
def calc_subj_snr(noblinkdata, blinktimes, ao_eprime):
    subj_sess_snr = noblinkdata.apply(calc_sess_snr, args=(blinktimes, ao_eprime))
    subj_sess_snr = subj_sess_snr.T
    subj_sess_snr.index = pd.MultiIndex.from_tuples(subj_sess_snr.index)    
    return pd.DataFrame(subj_sess_snr)
    
    
def event_waveforms(trg_events, std_events, pupilprofile, blinks, tpre=.5, pltpre=2, pltpost=4):
    subj_events = pd.DataFrame(columns=["Subject","Session","Condition","Time","Dilation"])
    for trg_event in trg_events:
        base_event = trg_event - pd.tseries.offsets.relativedelta(seconds=tpre)
        pltpre_event = trg_event - pd.tseries.offsets.relativedelta(seconds=pltpre)
        baseline = pupilprofile[base_event:trg_event].mean()
        post_event = trg_event + pd.tseries.offsets.relativedelta(seconds=pltpost)
        if blinks[trg_event:post_event].mean() >= .5:
            print 'Blinks during >50% of trial, skipping...'
            continue
        else:
            normed_trg_event = pupilprofile[pltpre_event:post_event] - baseline
            normed_trg_event.index = (normed_trg_event.index - trg_event) / np.timedelta64(1, 's')
            normed_trg_event.index.name = 'Time'
            normed_trg_event = pd.DataFrame(normed_trg_event).unstack().reset_index(name="Dilation")
            normed_trg_event = normed_trg_event.rename(columns={"level_0":"Subject","level_1":"Session"})
            normed_trg_event["Condition"] = "Target"
            subj_events = subj_events.append(normed_trg_event, ignore_index=True)
    for std_event in std_events:
        base_event = std_event - pd.tseries.offsets.relativedelta(seconds=tpre)
        pltpre_event = std_event - pd.tseries.offsets.relativedelta(seconds=pltpre)
        baseline = pupilprofile[base_event:std_event].mean()
        post_event = std_event + pd.tseries.offsets.relativedelta(seconds=pltpost)
        if blinks[std_event:post_event].mean() >= .5:
            print 'Blinks during >50% of trial, skipping...'
            continue
        else:
            normed_std_event = pupilprofile[pltpre_event:post_event] - baseline
            normed_std_event.index = (normed_std_event.index - std_event) / np.timedelta64(1, 's')
            normed_std_event.index.name = 'Time'
            normed_std_event = pd.DataFrame(normed_std_event).unstack().reset_index(name="Dilation")
            normed_std_event = normed_std_event.rename(columns={"level_0":"Subject","level_1":"Session"})
            normed_std_event["Condition"] = "Standard"
            subj_events = subj_events.append(normed_std_event, ignore_index=True)
    mean_subj = subj_events.groupby(['Subject','Condition','Session','Time'])['Dilation'].mean().reset_index()    
    return mean_subj

def plot_dilation(noblinkdata, blinktimes, ao_eprime):
    all_events = pd.DataFrame(columns=["Subject","Condition","Session","Time","Dilation"])
    for colname, col in noblinkdata.iteritems():
        noblink_series = noblinkdata[colname]
        blinktime_series = blinktimes[colname]
        subid, sess = noblink_series.name
        eprime_sess = ao_eprime[(ao_eprime['Subject_ID']==subid) & 
                                     (ao_eprime['Session']==sess)]
        eprime_sess.index = pd.to_datetime(list(eprime_sess.Tone_Onset), unit='ms')
        trg_events, std_events = get_events(noblink_series, eprime_sess)
        mean_subj = event_waveforms(trg_events, std_events, noblink_series, blinktime_series, tpre=.5, pltpre=2, pltpost=4)
        all_events = all_events.append(mean_subj)
        all_events = all_events.groupby(['Subject','Condition','Time'])['Dilation'].mean().reset_index()
    sns.tsplot(data=all_events, time="Time", 
       condition="Condition", unit="Subject", value="Dilation")
        
def proc_oddball(pupil_fname, behav_fname, outdir):
    parsed_df =  parse_ao.parse_pupil_data(pupil_fname, outdir)
    parsed_df[~parsed_df['Subject ID'].str.contains("LCIP99")]
    fulltrials = parsed_df['Measurement Duration'].str.replace('sec','').astype('float') > 290
    parsed_df = parsed_df.ix[fulltrials]
    ao_eprime = pd.read_csv(behav_fname)
    ao_eprime.loc[:,'Tone_Onset'] = ao_eprime['Tone_Onset'] + 3000
    sessioninfo = parsed_df[['Subject ID', 'Time', 'Session', 'Device ID', 'Eye Measured']]
    qc_dir = os.path.join(outdir,'QC') 
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    noblinkdata, blinktimes = clean_all(parsed_df, qc_dir=qc_dir)
    blinkpct = get_blink_pct(blinktimes)
    tstamp = datetime.now().strftime("%Y%m%d")
    blink_fname = "Blink_Percentages_" + tstamp + ".csv"
    outfile = os.path.join(qc_dir, blink_fname)
    blinkpct.to_csv(outfile, index=True, header=True)
    subj_sess_snr = calc_subj_snr(noblinkdata, blinktimes, ao_eprime)
    subj_info_snr = sessioninfo.merge(subj_sess_snr, left_on=['Subject ID','Session'], right_index=True)
    subj_info_snr = subj_info_snr.merge(pd.DataFrame(blinkpct), left_on=['Subject ID','Session'], right_index=True)
    snr_fname = 'LCIP_Oddball_SNR_' + tstamp + '.csv'
    outfile = os.path.join(outdir, snr_fname)
    subj_info_snr.to_csv(outfile, index=False, header=True)    
    plot_dilation(noblinkdata, blinktimes, ao_eprime)
    
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print 'USAGE: %s <raw pupil file> <behavioral file> <output directory>' % os.path.basename(sys.argv[0])
        print 'Takes pupillometer data text file and behavioral data as inputs.'
        print 'Parses pupil data, removes artifacts and calculates signal to noise'
        print 'ratio (target vs. non-targets). Outputs csv files for use in'
        print 'further analysis.'
    else:
        pupil_fname = sys.argv[1]
        behav_fname = sys.argv[2]
        outdir = sys.argv[3]
        proc_oddball(pupil_fname, behav_fname, outdir)


#########################################################
#                         TESTING                       #
#########################################################
pupil_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/raw/R_20161207_1224_20161207_1311.dat.txt'
behav_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/raw/oddball/OddballP300_LCI_Pilot_AllSubjects_12232016.csv'
outdir = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data'
