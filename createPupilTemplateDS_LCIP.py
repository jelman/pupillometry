# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 14:08:17 2016

@author: jelman
"""

import pandas as pd
import numpy as np


def rotate(strg,n):
    """ Create function to rotate characters in a string from front to back"""
    return strg[n:] + strg[:n]

def create_trialnames(df):
    letters = ['a','b','c','d']
    nrow = df.shape[0]
    letterseries = pd.Series(letters[:nrow], index=df.index)
    newtrials = df.Load.astype('str') + letterseries
    df['Bx Trial'] = newtrials
    return df

########################   Set filenames    ###############################
pupil_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/LCIP_DS_Raw_20171205_LCIP001-028.csv'
behav_fname = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/behavioral/PupilTestingData.xlsx'
outfile = '/home/jelman/netshare/VETSA_NAS/PROJ/LCIP/data/pupillometry/task_data/DS_30Hz_VETSA3_12052017.csv'
###########################################################################

### Load Data ###

# Load parsed pupil data
pupildf = pd.read_csv(pupil_fname, sep=",")


# Load behavioral data with trial times and scan quality
behavdflong = pd.read_excel(behav_fname, sheetname="DS")


### Behavioral Data ###

# Rename subject ID to match pupil data
behavdflong['Subject'] = 'LCIP' + behavdflong['Subject'].astype('str').str.zfill(3)

# Sort by subject id
behavdflong = behavdflong.sort_values(by=["Subject","Trial"])

# Rename subject ID, time and trials fields
behavdflong = behavdflong.rename(columns={"Subject":"ID",
                                          "Time":"Pupil Trial Time",
                                          "Trial":"Trial #"})
# Rename time field
behavdflong = behavdflong.rename(columns={"Time":"Pupil Trial Time"})

# Add Pair ID field
behavdflong['Pair ID'] = behavdflong['ID'].str.replace("[AB]", "")

# Convert to string
behavdflong['Pupil Trial Time'] = behavdflong['Pupil Trial Time'].astype('str')

# Create Task field and Bx Trial field
behavdflong['Task (DS)'] = 'DS' + behavdflong['Trial #'].astype('str')
grpdf = behavdflong.groupby(['ID','Load'])
behavdflong = grpdf.apply(create_trialnames)

# Specify bad data trials and filter out unacquired trials
behavdflong = behavdflong.dropna(axis=0, subset=["Pupil Trial Time"])
behavdflong.loc[behavdflong['Use']=='Bad','Pupil Trial Time'] = "bad data"

### Pupil Data ###

# Split pupil measurements into columnes
pupilprofile = pupildf.pop('Pupil Profile')
pupilprofile = pd.DataFrame(pupilprofile.str.split('\t').values.tolist())
pupilprofile.columns = pupilprofile.columns.astype('str')
pupildf = pd.concat([pupildf, pupilprofile], axis=1)

# Remove date from pupil trial time columns
pupildf['Time'] = pupildf['Time'].str.split(' ').str.get(1)

# Calculate Min and Max fields
# Insert min
pupildf.insert(pupildf.columns.get_loc("0"), 'Min', pupildf.loc[:,"0":].min(axis=1))
# Insert max
pupildf.insert(pupildf.columns.get_loc("0"), 'Max', pupildf.loc[:,"0":].max(axis=1))

# Rename columns
pupildf = pupildf.rename(columns={"Subject ID":"ID", "Eye Measured":"LR", "Time":"Pupil Trial Time"})
ndatacols = pupildf.columns[pupildf.columns.get_loc("0"):].shape[0]
datacolnames = [str(x +1) + 'st data pt' for x in range(ndatacols)]
pupildf.columns = list(pupildf.columns[:pupildf.columns.get_loc("0")]) + datacolnames


### Merge data ###
mergeddf = behavdflong.merge(pupildf, how="inner", on=["ID","Pupil Trial Time"])

# Get number of twins (paired or singleton)
mergeddf['Matched?'] = ""

templateDS = mergeddf[['ID','LR','Pair ID','Matched?','Task (DS)',
                       'Pupil Trial Time','Trial #', 'Bx Trial','Min',
                       'Max']+datacolnames]
templateDS.to_csv(outfile, index=False)
