#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Created on Tue Dec  5 13:24:25 2017

This is the fourth script to run on data from the pupillometer. It takes as 
input the stats output by Granholm lab matlab script. This output contains 
the baseline pupil diameter, raw diameter at each second for loads 3/6/9, and 
the baseline corrected diameter at each second. 
The output is the baseline diameter and data from the timepoint of interest for
each condition in both wide and long formats.
@author: jelman
"""

import pandas as pd
import re

########################   Set filenames    ###############################
pupilfile = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/stats/DS_30Hz_VETSA3_12052017_stats.csv'
behavfile = "/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/LCIP_ExptInfo.csv"
wideout = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/stats/LCIP_DS_Stats_wide_20171129.csv'
longout = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/stats/LCIP_DS_Stats_long_20171129.csv'
###########################################################################

# Load data 
pupildf = pd.read_csv(pupilfile)
behavdf = pd.read_csv(behavfile)

# Select columns of interest from behavioral session data
behavdf = behavdf[["Subject","Date","Device","Eye"]]

# Define and select columns of interest from pupil data
cols = ['ID', 'ntr3', 'AvgAmpBaseline3', 'AvgAmpTW4_3', 'AvgNormAmpTW4_3',
        'ntr6', 'AvgAmpBaseline6', 'AvgAmpTW7_6', 'AvgNormAmpTW7_6',
        'ntr9', 'AvgAmpBaseline9', 'AvgAmpTW10_9', 'AvgNormAmpTW10_9']
pupildf = pupildf[cols]
# Rename columns for easier reshape into long format
newcols = [re.sub('TW([47]|10)_','',word) for word in cols]
pupildf.columns = newcols
pupildf.rename(columns={"ID":"Subject"}, inplace=True)

# Merge pupil and behavioral session data
widedf = behavdf.merge(pupildf, on="Subject")

# Save out wide format data
widedf.to_csv(wideout, sep=",", index=False)


# Reshape into long format. Each row corresponds to one level of load per subject
longdf = pd.wide_to_long(widedf, ['ntr','AvgAmpBaseline','AvgAmp','AvgNormAmp'], i='Subject',j='Load')

# Sort by subject then load
longdf = longdf.sort_index()

# Save out long format data
longdf.to_csv(longout, index=True)