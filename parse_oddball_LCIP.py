# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:40:01 2016

@author: jelman
"""

import itertools
import os
import re
import sys
import pandas as pd



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


def join_multilines(lines):
    # Add blank space to end so that loop will include last element in list
    lines.append('')
    # Create iterator
    iterlines = iter(lines)
    # Initialize results container
    results = []
    # Initialize previous line variable
    prev = next(iterlines)
    # Begin looping. Join multi-line elements and place breakpoints between trials
    for line in iterlines:
        if re.search("= $", prev):
            results.append(prev+line)
            prev = next(iterlines)
        elif (prev=='' and line==''):
            results.append('BREAK')
            prev = next(iterlines)
            continue
        else:
            results.append(prev)
            prev = line
    return results


def split_trial_lists(lines):
    # Break list into sublists of trials
    sublists = [list(x[1]) for x in itertools.groupby(lines, lambda x: x=='BREAK') if not x[0]]
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
        outfile = os.path.join(outdir, 'LCIP_Oddball_Parsed.csv')
        try:
            parsed_df.to_csv(outfile, index=False)
            print 'Saved to %s' %(outfile)
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
