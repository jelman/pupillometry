# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 09:30:17 2016

@author: jelman
"""

import os
import pandas as pd
from glob import glob

datadir = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/raw/oddball'
globstr = os.path.join(datadir, 'OddballP300_LCI_Pilot_*.csv')
filelist = glob(globstr)

allsubs = pd.DataFrame()

for fname in filelist:   
    subdat = pd.read_csv(fname)
    allsubs = allsubs.append(subdat)

allsubs['Subject_ID'] = allsubs['Subject_ID'].astype(str).apply(lambda x: x.zfill(3))
allsubs['Subject_ID'] = 'LCIP' + allsubs['Subject_ID']

allsubs = allsubs.sort(columns=['Subject_ID', 'Session'])

outfile = os.path.join(datadir, 'OddballP300_LCI_Pilot_AllSubjects.csv')
allsubs.to_csv(outfile, index=False)