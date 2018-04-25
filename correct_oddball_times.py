#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 11:33:10 2017

@author: jelman

Correct timing of oddball trials in which pupil tracking did not start 
until the first trial. Takes as input the paresed oddball data and adds 
3000ms to all times on affected sessions.
"""

import pandas as pd

infile = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/LCIP_Oddball_Parsed_20171129_LCIP001-028.csv'
outfile = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/LCIP_Oddball_Parsed_20171129_LCIP001-028_corrected.csv'
parsed_df = pd.read_csv(infile)



for subject, session in [('LCIP024',1), ('LCIP025',1), ('LCIP026',2)]:
    sess_df = parsed_df[(parsed_df['Subject ID']==subject) & (parsed_df['Session']==session)]
    timeprofile = sess_df['Time Profile'].str.split('\t').values[0]
    corrected_time = [str(float(t) + 3) for t in timeprofile]
    corrected_time = '\t'.join(corrected_time)
    parsed_df.loc[(parsed_df['Subject ID']==subject) & (parsed_df['Session']==session), 'Time Profile'] = corrected_time
    
parsed_df.to_csv(outfile, index=False)