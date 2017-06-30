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


def calc_trial_dilations(events, pupilprofile, blinks, tpre=.5, tpost=2):
    mean_dilations = []
    max_dilations = []
    sd_dilations = []
    for event in events:
        pre_event = event - pd.tseries.offsets.relativedelta(seconds=tpre)
        baseline = pupilprofile[pre_event:event].mean()
        post_event = event + pd.tseries.offsets.relativedelta(seconds=tpost)
        if blinks[pre_event:post_event].mean() >= .33:
            print 'Blinks during >33% of trial, skipping...'
            continue
        else:
            normed_post_event = pupilprofile[event:post_event] - baseline
            mean_dil = normed_post_event.mean()
            max_dil = normed_post_event.max()
            sd_dil = normed_post_event.std()
            mean_dilations.append(mean_dil)
            max_dilations.append(max_dil)
            sd_dilations.append(sd_dil)
    return np.nanmean(mean_dilations), np.nanmean(max_dilations), np.nanmean(sd_dilations)


def h_pupil(t,n=10.1,t_max=930,f=1./(10**27) ):
  # n+1 = number of laters
  # t_max = response maximum
  # f = scaling factor
  h = f*(t**n)*np.exp(-n*t/t_max)
  h[0] = 0
  return h

   
def calc_prf_fit(pupilprofile, blinks, eprime_sess, tpre=.2, tpost=4):
    trg_events, std_events = get_events(pupilprofile, eprime_sess)
    tuneparlist = []
    for event in trg_events:
        pre_event = event - pd.tseries.offsets.relativedelta(seconds=tpre)
        baseline = pupilprofile[pre_event:event].mean()
        post_event = event + pd.tseries.offsets.relativedelta(seconds=tpost)
        if blinks[pre_event:post_event].mean() >= .33:
            print 'Blinks during >33% of trial, skipping...'
            continue
        else:
            normed_post_event = pupilprofile[event:post_event] - baseline
            t = np.linspace(0,tpost,len(normed_post_event)) * 1000
            prf = h_pupil(t)
            tunepar = np.sum(normed_post_event * prf)
            tuneparlist.append(tunepar)
    return np.median(tuneparlist)
    
    
    
def tune_offset(noblink_series, blinktime_series, eprime_sess):
    offset_stats = {}
    for offset in np.arange(-2500, 2500, 250):
        eprime_sess_copy = eprime_sess.copy()
        eprime_sess_copy.loc[:,'Tone_Onset'] = eprime_sess_copy.Tone_Onset + offset
        eprime_sess_copy.index = pd.to_datetime(list(eprime_sess_copy.Tone_Onset), unit='ms')
        offset_stats[offset] = calc_prf_fit(noblink_series, blinktime_series, eprime_sess_copy)
    tuneSeries = pd.Series(offset_stats)
    maxparam = tuneSeries.argmax()
    return maxparam
    
    
def get_blink_pct(blinktimes):
    blinkpct = blinktimes.apply(np.mean)
    blinkpct.name = 'blink_pct'
    blinkpct.index = pd.MultiIndex.from_tuples(blinkpct.index) 
    blinkpct.index.names = ['Subject ID', 'Session']
    return blinkpct


def calc_sess_stats(noblink_series, blinktimes, ao_eprime):
    """
    DIFF: Target max - Standard max
    CNR1: Target max / Standard SD
    CNR2: (Target max - Standard max) / Standard SD
    CNR3: Target SD / Standard SD
    """
    subid, sess = noblink_series.name
    blinktime_series = blinktimes[noblink_series.name]
    eprime_sess = ao_eprime[(ao_eprime['Subject_ID']==subid) & 
                             (ao_eprime['Session']==sess)] 
    offset = tune_offset(noblink_series, blinktime_series, eprime_sess)
    eprime_sess.loc[:,'Tone_Onset'] = eprime_sess.Tone_Onset + offset
    eprime_sess.index = pd.to_datetime(list(eprime_sess.Tone_Onset), unit='ms')
    trg_events, std_events = get_events(noblink_series, eprime_sess)
    trg_mean_dil, trg_max_dil, trg_sd_dil = calc_trial_dilations(trg_events, noblink_series, blinktime_series)
    std_mean_dil, std_max_dil, std_sd_dil = calc_trial_dilations(std_events, noblink_series, blinktime_series)
    sess_diff = trg_max_dil - std_max_dil
    sess_cnr1 = trg_max_dil / std_sd_dil
    sess_cnr2 = (trg_max_dil - std_max_dil) / std_sd_dil
    sess_cnr3 = trg_sd_dil / std_sd_dil
    resultdict = dict(Trg_max = trg_max_dil, Trg_mean = trg_mean_dil,
                      Std_max = std_max_dil, Std_mean = std_mean_dil,
                      DIFF = sess_diff, CNR1 = sess_cnr1, 
                      CNR2 = sess_cnr2, CNR3 = sess_cnr3, 
                      Offset=offset)
    return pd.Series(resultdict)
  
    
def calc_subj_stats(noblinkdata, blinktimes, ao_eprime):
    subj_sess_snr = noblinkdata.apply(calc_sess_stats, args=(blinktimes, ao_eprime))
    subj_sess_snr = subj_sess_snr.T
    subj_sess_snr.index = pd.MultiIndex.from_tuples(subj_sess_snr.index)    
    return pd.DataFrame(subj_sess_snr)
    

def event_waveform(event, pupilprofile, blinks, condition, tpre=.5, pltpre=1, pltpost=4):
    base_event = event - pd.tseries.offsets.relativedelta(seconds=tpre)
    pltpre_event = event - pd.tseries.offsets.relativedelta(seconds=pltpre)
    baseline = pupilprofile[base_event:event].mean()
    post_event = event + pd.tseries.offsets.relativedelta(seconds=pltpost)
    if blinks[base_event:post_event].mean() >= .33:
        print 'Blinks during >33% of trial, skipping...'
    else:
        normed_event = pupilprofile[pltpre_event:post_event] - baseline
        normed_event.index = (normed_event.index - event) / np.timedelta64(1, 's')
        normed_event.index.name = 'Time'
        normed_event = pd.DataFrame(normed_event).unstack().reset_index(name="Dilation")
        normed_event = normed_event.rename(columns={"level_0":"Subject","level_1":"Session"})
        normed_event["Condition"] = condition
        return normed_event
  
        
def subj_waveforms(trg_events, std_events, pupilprofile, blinks, **kwargs):
    subj_events = pd.DataFrame(columns=["Subject","Session","Condition","Time","Dilation"])
    for trg_event in trg_events:
        normed_event = event_waveform(trg_event, pupilprofile, blinks, "Target")
        subj_events = subj_events.append(normed_event, ignore_index=True)
    for std_event in std_events:
        normed_event = event_waveform(std_event, pupilprofile, blinks, "Standard")
        subj_events = subj_events.append(normed_event, ignore_index=True)
    mean_subj = subj_events.groupby(['Subject','Condition','Session','Time'])['Dilation'].mean().reset_index()    
    return mean_subj

    
def plot_dilation(noblinkdata, blinktimes, ao_eprime, outdir, offset_info=None):
    all_events = pd.DataFrame(columns=["Subject","Condition","Session","Time","Dilation"])
    for colname, col in noblinkdata.iteritems():
        noblink_series = noblinkdata[colname]
        blinktime_series = blinktimes[colname]
        subid, sess = noblink_series.name
        eprime_sess = ao_eprime[(ao_eprime['Subject_ID']==subid) & 
                                     (ao_eprime['Session']==sess)]
        if offset_info is not None:
            offset = offset_info.ix[(subid, sess), "Offset"]
            eprime_sess.loc[:,'Tone_Onset'] = eprime_sess.Tone_Onset + offset         
        eprime_sess.index = pd.to_datetime(list(eprime_sess.Tone_Onset), unit='ms')
        trg_events, std_events = get_events(noblink_series, eprime_sess)
        mean_subj = subj_waveforms(trg_events, std_events, noblink_series, blinktime_series)
        p = sns.tsplot(data=mean_subj, time="Time", 
                       condition="Condition", unit="Session", value="Dilation").figure
        outfile = os.path.join(outdir, str(subid) + '_' + str(sess) + '_PSTC.png')
        p.savefig(outfile)  
        p.clear()
        all_events = all_events.append(mean_subj)
        all_events = all_events.groupby(['Subject','Condition','Time'])['Dilation'].mean().reset_index()
    p = sns.tsplot(data=all_events, time="Time", 
                   condition="Condition", unit="Subject", value="Dilation")
    outfile = os.path.join(outdir, 'AllSubjects_PSTC.png')
    p.figure.savefig(outfile)
     
    
def proc_oddball(pupil_fname, behav_fname, outdir):
    parsed_df =  parse_ao.parse_pupil_data(pupil_fname, outdir)
    parsed_df = parsed_df[~parsed_df['Subject ID'].str.contains("LCIP99")]
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
    subj_sess_stats = calc_subj_stats(noblinkdata, blinktimes, ao_eprime)
    subj_info_stats = sessioninfo.merge(subj_sess_stats, left_on=['Subject ID','Session'], right_index=True)
    subj_info_stats = subj_info_stats.merge(pd.DataFrame(blinkpct), left_on=['Subject ID','Session'], right_index=True)
    stats_fname = 'LCIP_Oddball_Stats_' + tstamp + '.csv'
    outfile = os.path.join(outdir, stats_fname)
    subj_info_stats.to_csv(outfile, index=False, header=True)    
    plot_dilation(noblinkdata, blinktimes, ao_eprime, qc_dir, subj_sess_stats)
    
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
#pupil_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/raw/R_20161207_1224_20161207_1311.dat.txt'
#behav_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/raw/oddball/OddballP300_LCI_Pilot_AllSubjects_12232016.csv'
#outdir = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data'
