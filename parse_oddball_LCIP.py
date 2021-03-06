# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:40:01 2016

@author: jelman

This is the first script to run on raw data file from pupillometer.
It will parse the text file and output a spreadsheet with one row
per trial that can be used for input to proc scripts that will perform
cleaning and averaging.

This script can parse dauditory oddball trials, and is specific to the
LCIP project.
"""

import itertools
import os
import re
import sys
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
    # Initialize list
    sublistsAO = []
    # Append data list
    for sublist in trial_lists:
        if 'Finite Measurement Duration = 0' in sublist:
            sublistsAO.append(sublist)
        else:
            continue
    return sublistsAO


def create_ao_df(sublistsAO):
    # Create pandas dataframe from auditory oddball sublists
    dictAOlist = []
    for sublistAO in sublistsAO:
        dictAO = {}
        for item in sublistAO:
            key, val = item.split(' = ')
            dictAO[key.strip()] = val.strip()
        dictAOlist.append(dictAO)
    ao_df = pd.DataFrame(dictAOlist)
    return ao_df


def get_session_nums(df):
   sessionnums = df.groupby('Subject ID').apply(lambda x: [(i+1) for i in range(x.shape[0])])
   sessionnums = list(itertools.chain.from_iterable(sessionnums))
   return sessionnums


def parse_pupil_data(filename, outdir=None):
    lines = read_file(filename)
    clean_lines = clean_text(lines)
    clean_lines = join_multilines(clean_lines)
    trial_lists = split_trial_lists(clean_lines)
    lcip_lists = find_lcip_subs(trial_lists)
    sublistsAO = get_task_lists(lcip_lists)
    # Auditory Oddball
    parsed_df = create_ao_df(sublistsAO)
    parsed_df = parsed_df.sort_values(['Subject ID', 'Time']).reset_index(drop=True)
    parsed_df['Session'] = get_session_nums(parsed_df)
    if outdir:
        timestamp = datetime.now().strftime("%Y%m%d")
        outname = 'LCIP_Oddball_Parsed_' + timestamp + '.csv'
        outfile = os.path.join(outdir, outname)
        try:
            parsed_df.to_csv(outfile, index=False)
            print 'Saved to %s' %(outfile)
            return parsed_df
        except IOError:
            print 'Error: Could not save to %s' %(outfile)
    else:
        return parsed_df


if __name__ == '__main__':

    if len(sys.argv) == 1:
        print 'USAGE: %s <file> <output directory>' % os.path.basename(sys.argv[0])
        print 'Takes pupillometer data text file as input. Extracts oddball'
        print 'trials and parses for further processing.'
    else:
        filename = sys.argv[1]
        outdir = sys.argv[2]
        parse_pupil_data(filename, outdir)
