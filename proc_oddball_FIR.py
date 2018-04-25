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
from datetime import datetime
import scipy as sci
import ast
sys.path.append('/home/jelman/netshare/K/code/Pupillometry/FIRDeconvolution/src')
from FIRDeconvolution import FIRDeconvolution


def get_offset_params():
    params = Parameters()
    # Offset parameters
    params.add('offset', value=0, min=-2.5, max=2.5)
    return params


def get_IRF_params():
    params = Parameters()
    # Pupil impulse response parameters
    params.add('s1', value=5000., min=1e-28, max=1e8)
    params.add('n1', value=10.1, min=6, max=15)
    params.add('tmax1', value=0.930, min=0.1, max=3.0)
    return params


def get_events(pupilprofile, eprime):
    trg_onsets = eprime.index[eprime['Tone']==2]
    std_onsets = eprime.index[eprime['Tone']==1][1:]
    # Get event indices
    trg_event_idx = pupilprofile.index.searchsorted(trg_onsets)
    std_event_idx = pupilprofile.index.searchsorted(std_onsets)
    trg_events = pupilprofile.index[trg_event_idx]
    std_events = pupilprofile.index[std_event_idx]
    
    return trg_events, std_events


def hp_filter(signal, sample_rate):
    # High pass:
    hp = 0.01
    hp_cof_sample = hp /  (sample_rate / 2)
    bhp, ahp = sci.signal.butter(3, hp_cof_sample, btype='high')
    signal_hp = sci.signal.filtfilt(bhp, ahp, signal)
    return signal_hp


def single_pupil_IRF(params, x):
    s1 = params['s1']
    n1 = params['n1']
    tmax1 = params['tmax1']
    return s1 * ((x**n1) * (np.e**((-n1*x)/tmax1)))

def single_pupil_IRF_ls(params, x, data):
    s1 = params['s1'].value
    n1 = params['n1'].value
    tmax1 = params['tmax1'].value
    model = s1 * ((x**n1) * (np.e**((-n1*x)/tmax1)))
    return model - data


def add_offset(events, offset):
    events_new = ((events - pd.Timestamp(0)).astype('timedelta64[ms]') / 1000.) + offset
    return events_new


def run_deconvolve(offset, events, event_names, signal_hp, interval):
    """
    Run FIR deconvolution.
    offset (float) : offset in sec between eprime and pupillometer timestamps
    events (list) : list 
    """
    a = FIRDeconvolution(signal=signal_hp, 
                         events=events, event_names=event_names, sample_frequency=15.0, 
                         deconvolution_frequency = 15.,
                         deconvolution_interval=interval)
    a.create_design_matrix()
    a.regress()
    a.betas_for_events()
    a.calculate_rsq()
    return a


def offset_fit_all(offset, trg_events, std_events, signal_hp, interval):
    trg_new = add_offset(trg_events, offset)
    std_new = add_offset(std_events, offset)
    events = [np.sort(np.append(np.array(trg_new), np.array(std_new)))]
    event_names = ["allevents"]
    a = run_deconvolve(offset, events, event_names, signal_hp, interval)
    response = np.array(a.betas_per_event_type[0]).ravel()
    # baseline the kernels:
    response = response - response[0].mean()
    trial_x = np.linspace(0, interval[1], len(response))
    irf_params = get_IRF_params()
    result = minimize(single_pupil_IRF_ls, irf_params, method='powell', args=(trial_x, response))
    kernel = single_pupil_IRF(result.params, trial_x)
    return result, response, kernel


def offset_fit_all_ls(params, trg_events, std_events, signal_hp, interval):
    offset = params['offset'].value
    trg_new = add_offset(trg_events, offset)
    std_new = add_offset(std_events, offset)
    events = [np.sort(np.append(np.array(trg_new), np.array(std_new)))]
    event_names = ["allevents"]
    a = run_deconvolve(offset, events, event_names, signal_hp, interval)
    response = np.array(a.betas_per_event_type[0]).ravel()
    # baseline the kernels:
    response = response - response[0].mean()
    trial_x = np.linspace(0, interval[1], len(response))
    irf_params = get_IRF_params()
    result = minimize(single_pupil_IRF_ls, irf_params, method='powell', args=(trial_x, response))
    kernel = single_pupil_IRF(result.params, trial_x)
    return np.sum((response - kernel)**2)


def offset_std_all(offset, trg_events, std_events, signal_hp, interval):
    trg_new = add_offset(trg_events, offset)
    std_new = add_offset(std_events, offset)
    events = [np.sort(np.append(np.array(trg_new), np.array(std_new)))]
    event_names = ["allevents"]
    a = run_deconvolve(offset, events, event_names, signal_hp, interval)
    response = np.array(a.betas_per_event_type[0]).ravel()
    # baseline the kernels:
    response = response - response[0].mean()
    trial_x = np.linspace(0, interval[1], len(response))
    irf_params = get_IRF_params()
    kernel = single_pupil_IRF(irf_params.valuesdict(), trial_x)
    return response, kernel


def offset_std_all_ls(params, trg_events, std_events, signal_hp, interval):
    offset = params['offset'].value
    trg_new = add_offset(trg_events, offset)
    std_new = add_offset(std_events, offset)
    events = [np.sort(np.append(np.array(trg_new), np.array(std_new)))]
    event_names = ["allevents"]
    a = run_deconvolve(offset, events, event_names, signal_hp, interval)
    response = np.array(a.betas_per_event_type[0]).ravel()
    # baseline the kernels:
    response = response - response[0].mean()
    trial_x = np.linspace(0, interval[1], len(response))
    irf_params = get_IRF_params()
    kernel = single_pupil_IRF(irf_params.valuesdict(), trial_x)
    return np.sum((response - kernel)**2)


def offset_fit_trg(offset, trg_events, signal_hp, interval):
    trg_new = ((trg_events - pd.Timestamp(0)).astype('timedelta64[ms]') / 1000.) + offset
    events = [trg_new]
    event_names = ["targets"]
    a = run_deconvolve(offset, events, event_names, signal_hp, interval)
    response = np.array(a.betas_per_event_type[0]).ravel()
    # baseline the kernels:
    response = response - response[0].mean()
    trial_x = np.linspace(0, interval[1], len(response))
    irf_params = get_IRF_params()
    result = minimize(single_pupil_IRF_ls, irf_params, method='powell', args=(trial_x, response))
    kernel = single_pupil_IRF(result.params, trial_x)
    return result, response, kernel


def offset_fit_trg_ls(params, trg_events, signal_hp, interval):
    offset = params['offset'].value
    trg_new = ((trg_events - pd.Timestamp(0)).astype('timedelta64[ms]') / 1000.) + offset
    events = [trg_new]
    event_names = ["targets"]
    a = run_deconvolve(offset, events, event_names, signal_hp, interval)
    response = np.array(a.betas_per_event_type[0]).ravel()
    # baseline the kernels:
    response = response - response[0].mean()
    trial_x = np.linspace(0, interval[1], len(response))
    irf_params = get_IRF_params()
    result = minimize(single_pupil_IRF_ls, irf_params, method='powell', args=(trial_x, response))
    kernel = single_pupil_IRF(result.params, trial_x)
    return np.sum((response - kernel)**2)


def offset_std_trg(offset, trg_events, signal_hp, interval):
    trg_new = ((trg_events - pd.Timestamp(0)).astype('timedelta64[ms]') / 1000.) + offset
    events = [trg_new]
    event_names = ["targets"]
    a = run_deconvolve(offset, events, event_names, signal_hp, interval)
    response = np.array(a.betas_per_event_type[0]).ravel()
    # baseline the kernels:
    response = response - response[0].mean()
    trial_x = np.linspace(0, interval[1], len(response))
    irf_params = get_IRF_params()
    kernel = single_pupil_IRF(irf_params.valuesdict(), trial_x)
    return response, kernel


def offset_std_trg_ls(params, trg_events, signal_hp, interval):
    offset = params['offset'].value
    trg_new = ((trg_events - pd.Timestamp(0)).astype('timedelta64[ms]') / 1000.) + offset
    events = [trg_new]
    event_names = ["targets"]
    a = run_deconvolve(offset, events, event_names, signal_hp, interval)
    response = np.array(a.betas_per_event_type[0]).ravel()
    # baseline the kernels:
    response = response - response[0].mean()
    trial_x = np.linspace(0, interval[1], len(response))
    irf_params = get_IRF_params()
    kernel = single_pupil_IRF(irf_params.valuesdict(), trial_x)
    return np.sum((response - kernel)**2)


def plot_response_kernel(response, kernel, interval, title, plotdir):
    # Turn interactive plotting off
    plt.ioff()
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
    
def create_regressors(series_time_idx, trg_events, std_events, offset, kernel):
    offsettime  = pd.Timedelta(offset, unit="s")
    trg_idx = [series_time_idx.get_loc(i-offsettime, method="nearest") for i in trg_events]
    std_idx = [series_time_idx.get_loc(i-offsettime, method="nearest") for i in std_events]
    trg_reg = np.zeros(len(series_time_idx))
    trg_reg[trg_idx] = 1
    trg_reg_conv = sci.signal.fftconvolve(trg_reg, kernel, 'full')[:-(len(kernel)-1)]
    std_reg = np.zeros(len(series_time_idx))
    std_reg[std_idx] = 1
    std_reg_conv = sci.signal.fftconvolve(std_reg, kernel, 'full')[:-(len(kernel)-1)]
    regs = [trg_reg_conv, std_reg_conv]
    return regs


def run_glm(signal_hp, regs):
    intercept = np.ones(len(signal_hp))
    design_matrix = np.matrix(np.vstack((intercept,[reg for reg in regs]))).T
    betas = np.array(((design_matrix.T * design_matrix).I * design_matrix.T) * np.matrix(signal_hp).T).ravel()
    pred_events = [betas[i+1]*regs[i] for i in range(len(betas)-1)]
    return betas, pred_events


def calc_sess_stats(noblink_series, blinktimes, ao_eprime, fit=True):
    """
    DIFF: Target max - Standard max
    CNR1: Target max / Standard SD
    CNR2: (Target max - Standard max) / Standard SD
    CNR3: Target SD / Standard SD
    CNR4: (Target max - Standard max) / Standard max
    """
    subid, sess = noblink_series.name
    noblink_series.index = pd.to_datetime(noblink_series.index)
    blinktime_series = blinktimes[noblink_series.name]
    eprime_sess = ao_eprime[(ao_eprime['Subject_ID']==subid) & 
                             (ao_eprime['Session']==sess)] 
    noblink_series = noblink_series.dropna()
    blinktime_series = blinktime_series.dropna()
    eprime_sess.index = pd.to_datetime(list(eprime_sess.Tone_Onset), unit='ms')
    trg_events, std_events = get_events(noblink_series, eprime_sess)
    signal_hp = hp_filter(noblink_series, sample_rate)
    # Estimate offset between e-prime and pupillometer timestamps  
    offset_params = get_offset_params()
    if fit==True:
        offset_result = minimize(offset_fit_all_ls, offset_params, method='powell', args=(trg_events, std_events, signal_hp, interval))
#        offset_result = minimize(offset_fit_trg_ls, offset_params, method='powell', args=(trg_events, signal_hp, interval))
    else:
        offset_result = minimize(offset_std_all_ls, offset_params, method='powell', args=(trg_events, std_events, signal_hp, interval))
#        offset_result = minimize(offset_std_trg_ls, offset_params, method='powell', args=(trg_events, signal_hp, interval))
    offset = float(offset_result.params['offset'])
    # Deconvolve and fit pupil IRF kernel
    if fit==True:
        result, response, kernel = offset_fit_all(offset, trg_events, std_events, signal_hp, interval)
#        result, response, kernel = offset_fit_trg(offset, trg_events, signal_hp, interval)
    else:
        response, kernel = offset_std_all(offset, trg_events, std_events, signal_hp, interval)
#        response, kernel = offset_std_trg(offset, trg_events, signal_hp, interval)
    # Plot response and kernel
    title = "%s_%s" %(subid, sess)
    plot_response_kernel(response, kernel, interval, title, plotdir)
    # create regressors:
    regs = create_regressors(noblink_series.index, trg_events, std_events, offset, kernel)
    # GLM:
    betas, pred_events = run_glm(signal_hp, regs)    
    plot_glm(signal_hp, pred_events, title, plotdir)
    resultdict = {'intercept_beta' : betas[0],
                  'targets_beta': betas[1],
                  'standards_beta': betas[2],
                  'offset': offset}
    return pd.Series(resultdict)


def calc_subj_stats(noblinkdata, blinktimes, ao_eprime, fit):
    subj_sess_snr = noblinkdata.apply(calc_sess_stats, args=(blinktimes, ao_eprime, fit))
    subj_sess_snr = subj_sess_snr.T
    subj_sess_snr.index = pd.MultiIndex.from_tuples(subj_sess_snr.index, names=("Subject","Session"))    
    subj_sess_snr['contrast_beta'] = subj_sess_snr['targets_beta'] - subj_sess_snr['standards_beta']
    return subj_sess_snr



#----------------------------------------------------------------------------#

# Set filenames
noblinkdata_fname = "/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/LCIP_Oddball_NoBlinkData20180212.csv"
blinktimes_fname = "/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/LCIP_Oddball_BlinkTimes20180212.csv"
behav_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/raw/oddball/OddballP300_LCI_Pilot_AllSubjects_12062017.csv'
plotdir = "/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/STDplots/all"
outdir = "/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/stats"
tstamp = datetime.now().strftime("%Y%m%d")
stats_fname = 'LCIP_Oddball_STDStats_AllEventsFit' + tstamp + '.csv'

# Set parameters
interval = [-1, 4]
sample_rate = 15.
fit = False

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

subj_sess_stats = calc_subj_stats(noblinkdata, blinktimes, ao_eprime, fit)
subj_sess_stats = subj_sess_stats.reset_index()
outfile = os.path.join(outdir, stats_fname)
subj_sess_stats.to_csv(outfile, index=False, header=True)    

