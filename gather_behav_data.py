# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 09:30:17 2016

@author: jelman

Script to concatenate all subjects' auditory oddball behavioral files.
Needed for proc_oddball_LCIP.py to match dilations with trial types.
"""

import os
import pandas as pd
from glob import glob

datadir = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/raw/oddball'
# Find all behavioral files
globstr = os.path.join(datadir, 'OddballP300_LCI_Pilot_*.csv')
filelist = glob(globstr)
# Remove previously concatenated files from list
filelist = [ f for f in filelist if "AllSubjects" not in f ]

# Append subject files to master dataframe
allsubs = pd.DataFrame()
for fname in filelist:
    subdat = pd.read_csv(fname)
    allsubs = allsubs.append(subdat)
# Zero pad subject IDs
allsubs['Subject_ID'] = allsubs['Subject_ID'].astype(str).apply(lambda x: x.zfill(3))
allsubs['Subject_ID'] = 'LCIP' + allsubs['Subject_ID']
# Sort data
allsubs = allsubs.sort_values(by=['Subject_ID', 'Session'])
# Save out
outfile = os.path.join(datadir, 'OddballP300_LCI_Pilot_AllSubjects_12062017.csv')
allsubs.to_csv(outfile, index=False)
