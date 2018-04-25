#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 10:47:05 2018

@author: jelman
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
from lmfit import minimize, Parameters, Parameter, report_fit
from datetime import datetime, timedelta
import scipy as sci
import ast
from nitime.analysis import FilterAnalyzer, NormalizationAnalyzer
from nitime.timeseries import TimeSeries
import nitime.timeseries as ts
import nitime.analysis as nta
import nitime.viz as viz
from nipy.modalities.fmri.glm import GeneralLinearModel


def get_offset_params():
    params = Parameters()
    # Offset parameters
    params.add('offset', value=0, min=-2.5, max=2.5)
    return params


def get_events(pupilprofile, eprime):
    trg_onsets = eprime.index[eprime['Tone']==2]
    std_onsets = eprime.index[eprime['Tone']==1][1:]   
    return trg_onsets, std_onsets


def hp_filter(signal, sample_rate):
    F = FilterAnalyzer(signal, ub=7.5, lb=0.01)
    signal.data = F.fir.data
    return signal


def pupil_IRF(x):
    s1 = 1000
    n1 = 10.1
    tmax1 = 0.930
    return s1 * ((x**n1) * (np.e**((-n1*x)/tmax1)))


def add_offset(events, offset):
    events_new = events + timedelta(seconds=offset)
    return events_new



def plot_response_kernel(response, kernel, interval, title, plotdir):
    # Turn interactive plotting off
#    plt.ioff()
    # plot:
    trial_x = np.linspace(0, interval[1], len(response))
    f = plt.figure(figsize = (8,6))
    ax = f.add_subplot(1,1,1) 
    ax.plot(trial_x, response, label='all response')
    ax.plot(trial_x, kernel, label='all fit')
    ax.set_xlabel('Time from event (s)')
    ax.set_ylabel('Pupil size')
    ax.axhline(0,color = 'k', lw = 0.5, alpha = 0.5)
    ax.legend(loc=4)
    ax.set_title(title)
    sns.despine(offset=10)
    fname = title + '_FIRplot.png'
    outfile = os.path.join(plotdir, fname)
    f.savefig(outfile)
    plt.close(f)
 
    
def plot_glm(signal_hp, pred_events, title, plotdir):
    # Turn interactive plotting off
    plt.ioff()
    # plot:
    signal_x = np.arange(len(signal_hp))
    f = plt.figure(figsize = (10,8))
    ax = f.add_subplot(1,1,1) 
    ax.plot(signal_x, signal_hp, label='raw')
    ax.plot(signal_x, pred_events[0], label='targets')
    ax.plot(signal_x, pred_events[1], label='standards')
    ax.plot(signal_x, np.sum(pred_events, axis=0), label='all')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Pupil size')
    ax.axhline(0,color = 'k', lw = 0.5, alpha = 0.5)
    ax.legend()
    ax.set_title(title)
    sns.despine(offset=10)
    fname = title + '_glmplot.png'
    outfile = os.path.join(plotdir, fname)
    f.savefig(outfile)
    plt.close(f)
    
    
def convolve_reg(event_ts, kernel):
    return sci.signal.fftconvolve(event_ts, kernel, 'full')[:-(len(kernel)-1)]


def run_glm(signal_hp, regs):
    intercept = np.ones(len(signal_hp))
    design_matrix = np.matrix(np.vstack((intercept,[reg for reg in regs]))).T
    betas = np.array(((design_matrix.T * design_matrix).I * design_matrix.T) * np.matrix(signal_hp).T).ravel()
    pred_events = [betas[i+1]*regs[i] for i in range(len(betas)-1)]
    return betas, pred_events



def get_event_ts(signal_filt, events):
    event_secs = (events - pd.Timestamp(0)).total_seconds()
    event_idx = signal_filt.time.index_at(event_secs)
    event_reg = np.zeros(len(signal_filt))
    event_reg[event_idx] = 1
    event_ts = ts.TimeSeries(event_reg, sampling_rate=15., time_unit='s')
    return event_ts

    
    
def set_offset_ls(offset_params, trg_onsets, std_onsets, signal_filt):
    offset = offset_params['offset'].value
    trg_new = add_offset(trg_onsets, offset)
    std_new = add_offset(std_onsets, offset)
    trg_ts = get_event_ts(signal_filt, trg_new)    
    std_ts = get_event_ts(signal_filt, std_new)
    kernel_x = np.linspace(0, 4, len_et)
    kernel = pupil_IRF(kernel_x)
    trg_reg = convolve_reg(trg_ts, kernel)
    std_reg = convolve_reg(std_ts, kernel)
    intercept = np.ones_like(signal_filt.data)
    X = np.array(np.vstack((intercept, trg_reg, std_reg)).T)
#    X = np.atleast_2d(trg_reg).T
    Y = np.atleast_2d(signal_filt).T
    model = GeneralLinearModel(X)
    model.fit(Y, model='ar1')    
    return model.get_mse()[0]


def plot_event(signal_filt, trg_ts, std_ts, kernel, title):
    plt.ioff()
    trg_era = nta.EventRelatedAnalyzer(signal_filt, trg_ts, len_et=60, correct_baseline=True)
    std_era = nta.EventRelatedAnalyzer(signal_filt, std_ts, len_et=60, correct_baseline=True)
    f = viz.plot_tseries(
        ts.TimeSeries(data=np.vstack([kernel,trg_era.eta.data, std_era.eta.data]),
                  sampling_rate=trg_era.sampling_rate, time_unit='s'))
    fname = title + '_glmplot.png'
    outfile = os.path.join(plotdir, fname)
    f.savefig(outfile)
    plt.close(f)
   
    
def analyze_ts(trg_onsets, std_onsets, signal_filt, offset, title):
    trg_new = add_offset(trg_onsets, offset)
    std_new = add_offset(std_onsets, offset)
    trg_ts = get_event_ts(signal_filt, trg_new)    
    std_ts = get_event_ts(signal_filt, std_new)
    kernel_x = np.linspace(0, 4, len_et)
    kernel = pupil_IRF(kernel_x)
    trg_reg = convolve_reg(trg_ts, kernel)
    std_reg = convolve_reg(std_ts, kernel)
    plot_event(signal_filt, trg_ts, std_ts, kernel, title)
    intercept = np.ones_like(signal_filt.data)
    X = np.array(np.vstack((intercept, trg_reg, std_reg)).T)
    Y = np.atleast_2d(signal_filt).T
    model = GeneralLinearModel(X)
    model.fit(Y, model='ar1')    
    int_beta, trg_beta, std_beta = model.get_beta().T[0] 
    cval = [0,1,-1]
    con = model.contrast(cval)    
    zval = con.z_score()[0]
    resultdict = {'trg_beta':trg_beta, 'std_beta':std_beta, 'zval':zval}
    return resultdict

def calc_sess_stats(noblink_series, blinktimes, ao_eprime):
    """
    DIFF: Target max - Standard max
    CNR1: Target max / Standard SD
    CNR2: (Target max - Standard max) / Standard SD
    CNR3: Target SD / Standard SD
    CNR4: (Target max - Standard max) / Standard max
    """
    subid, sess = noblink_series.name
    title = "%s_%s" %(subid, sess)
    noblink_series.index = pd.to_datetime(noblink_series.index)
    blinktime_series = blinktimes[noblink_series.name]
    eprime_sess = ao_eprime[(ao_eprime['Subject_ID']==subid) & 
                             (ao_eprime['Session']==sess)] 
    pupilts = ts.TimeSeries(noblink_series.dropna(), sampling_rate=sample_rate)
    blinktime_series = blinktime_series.dropna()
    eprime_sess.index = pd.to_datetime(list(eprime_sess.Tone_Onset), unit='ms')
    signal_filt = hp_filter(pupilts, sample_rate)
    trg_onsets, std_onsets = get_events(noblink_series, eprime_sess)
    # Estimate offset between e-prime and pupillometer timestamps  
    offset_params = get_offset_params()
    offset_result = minimize(set_offset_ls, offset_params, method='powell', args=(trg_onsets, std_onsets, signal_filt))
    offset = float(offset_result.params['offset'])
    resultdict = analyze_ts(trg_onsets, std_onsets, signal_filt, offset, title)
    return pd.Series(resultdict)


def calc_subj_stats(noblinkdata, blinktimes, ao_eprime):
    subj_sess_snr = noblinkdata.apply(calc_sess_stats, args=(blinktimes, ao_eprime))
    subj_sess_snr = subj_sess_snr.T
    subj_sess_snr.index = pd.MultiIndex.from_tuples(subj_sess_snr.index, names=("Subject","Session"))    
    subj_sess_snr['contrast_beta'] = subj_sess_snr['trg_beta'] - subj_sess_snr['std_beta']
    return subj_sess_snr



#----------------------------------------------------------------------------#

# Set filenames
noblinkdata_fname = "/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/LCIP_Oddball_NoBlinkData20180212.csv"
blinktimes_fname = "/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/LCIP_Oddball_BlinkTimes20180212.csv"
behav_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/raw/oddball/OddballP300_LCI_Pilot_AllSubjects_12062017.csv'
plotdir = "/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/GLMplots"
outdir = "/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/stats"
tstamp = datetime.now().strftime("%Y%m%d")
stats_fname = 'LCIP_Oddball_GLM' + tstamp + '.csv'

# Set parameters
interval = [-1, 4]
sample_rate = 15.
len_et=60

############################
### Start running script ###
############################
# Load data
noblinkdata = pd.read_csv(noblinkdata_fname, index_col=0, parse_dates=True)  
blinktimes = pd.read_csv(blinktimes_fname, index_col=0, parse_dates=True)  
ao_eprime = pd.read_csv(behav_fname)
ao_eprime.loc[:,'Tone_Onset'] = ao_eprime['Tone_Onset'] + 3000
# Convert strings to tuples for heirarchical column names
columns = pd.MultiIndex.from_tuples([ast.literal_eval(item) for item in noblinkdata.columns])
noblinkdata.columns = columns
columns = pd.MultiIndex.from_tuples([ast.literal_eval(item) for item in blinktimes.columns])
blinktimes.columns = columns
if not os.path.exists(plotdir):
    os.makedirs(plotdir)    

subj_sess_stats = calc_subj_stats(noblinkdata, blinktimes, ao_eprime)
subj_sess_stats = subj_sess_stats.reset_index()
outfile = os.path.join(outdir, stats_fname)
subj_sess_stats.to_csv(outfile, index=False, header=True)    

