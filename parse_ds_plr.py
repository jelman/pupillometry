# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 10:24:32 2016

@author: jelman

This is the first script to run on raw data file from pupillometer.
It will parse the text file and output a spreadsheet with one row
per trial that can be used for input to proc scripts that will perform
cleaning and averaging.

This script can parse digit span and pupil light reflex trials.
"""

import itertools
import os
import re
import sys
import numpy as np
import pandas as pd
from datetime import datetime


def read_file(filename):
# Open file and read lines
    with open(filename, 'r') as f:
        lines = f.readlines()
    return lines


def clean_text(lines):
    # Strip newlines characters
    lines = [line.replace('\r\n', '') for line in lines]
    # Delete redundant 75% recovery time value reported on same line as latency
    lines = [re.sub(', 75%.*', '', line) for line in lines]
    return lines


def multi_delete(list_, indexes):
    indexes = sorted(list(indexes), reverse=True)
    for index in indexes:
        del list_[index]
    return list_

def join_multilines(lines):
    break_idx = [ i for i, item in enumerate(lines) if item.endswith("= ") ]
    for i in break_idx:
        lines[i] = ''.join(lines[i:i+2])
    del_idx = [x+1 for x in break_idx]
    lines = multi_delete(lines, del_idx)
    return lines
    
def split_trial_lists(lines):
    # Break list into sublists of trials
    sublists = [list(x[1]) for x in itertools.groupby(lines, lambda x: x=='') if not x[0]]
    return sublists


def find_lcip_subs(trial_lists):
    # Find and return only trials of LCIP subjects
    return [tlist for tlist in trial_lists for s in tlist if "Subject ID = LCIP" in s]


def get_task_lists(trial_lists):
    # Initialize list for digit span data and pupil light reflex
    sublistsDS = []
    sublistsPLR = []
    # Append data to digit span and pupil light reflex lists
    for sublist in trial_lists:
        if  'Measurement Duration = 15.000sec' in sublist:
            sublistsDS.append(sublist)
        elif 'Measurement Duration = 5.000sec' in sublist:
            sublistsPLR.append(sublist)
        else:
            continue
    return sublistsPLR, sublistsDS


def create_plr_df(sublistsPLR):
    dictPLRlist = []
    for sublistPLR in sublistsPLR:
        dictPLR = {}
        for item in sublistPLR:
            key, val = item.split(' = ')
            dictPLR[key.strip()] = val.strip()
        dictPLRlist.append(dictPLR)
    plr_df = pd.DataFrame(dictPLRlist)
    return plr_df


def create_ds_df(sublistsDS):
    # Create pandas dataframe from digit span sublists
    dictDSlist = []
    for sublistDS in sublistsDS:
        dictDS = {}
        for item in sublistDS:
            key, val = item.split(' = ')
            dictDS[key.strip()] = val.strip()
        dictDSlist.append(dictDS)
    ds_df = pd.DataFrame(dictDSlist)
    return ds_df


def save_to_csv(df, outfile, cols=None):
    if cols:
        df = df[cols]
    try:
        df.to_csv(outfile, index=False)
        print 'Saved to %s' %(outfile)
    except IOError:
        print 'Error: Could not save to %s' %(outfile)
    return outfile


def get_trial_nums(df):
   trialnums = df.groupby('ID').apply(lambda x: [(i+1) for i in range(x.shape[0])])
   trialnums = list(itertools.chain.from_iterable(trialnums))
   trialnums = [int(i) for i in trialnums]
   return trialnums


def parse_pupil_data(filename, outdir):
    lines = read_file(filename)
    clean_lines = clean_text(lines)
    clean_lines = join_multilines(clean_lines)
    trial_lists = split_trial_lists(clean_lines)
    lcip_lists = find_lcip_subs(trial_lists)
    sublistsPLR, sublistsDS = get_task_lists(lcip_lists)
    timestamp = datetime.now().strftime("%Y%m%d")
    # Pupil Light Reflex
    plr_df = create_plr_df(sublistsPLR)
    plr_df = plr_df.sort_values(['Subject ID','Time']).reset_index(drop=True)
    plrCols = ['Subject ID', 'Time', 'Device ID', 'Eye Measured', 'Record ID',
          'Profile Normal', 'Diameter', 'Measurement Duration', 'Mean/Max C. Vel',
          'dilation velocity', 'Lat', '75% recovery time', 'Pupil Profile']
    plr_fname = 'LCIP_PLR_Raw_' + timestamp + '.csv'
    plr_outfile = os.path.join(outdir,plr_fname)
    save_to_csv(plr_df, plr_outfile, plrCols)
    # Digit Span
    ds_df = create_ds_df(sublistsDS)
    ds_df = ds_df.sort_values(['Subject ID','Time']).reset_index(drop=True)
    dsCols = ['Subject ID', 'Time', 'Device ID', 'Eye Measured', 'Record ID',
              'Profile Normal', 'Diameter', 'Measurement Duration',
              'Pupil Profile']
    ds_fname = 'LCIP_DS_Raw_' + timestamp + '.csv'
    ds_outfile = os.path.join(outdir, ds_fname)
    save_to_csv(ds_df, ds_outfile, dsCols)



if __name__ == '__main__':


    if len(sys.argv) == 1:
        print 'USAGE: %s <file> <output directory>' % os.path.basename(sys.argv[0])
        print 'Takes pupillometer data text file outputs csv files for use'
        print 'in further analysis.'
    else:
        filename = sys.argv[1]
        outdir = sys.argv[2]
        parse_pupil_data(filename, outdir)
