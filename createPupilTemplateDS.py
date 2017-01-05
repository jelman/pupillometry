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

########################   Set filenames    ###############################
pupil_fname = '/home/jelman/netshare/K/data/Pupillometry/VETSA3/Pupillometry Data_DS_Data.csv'
behav_fname = '/home/jelman/netshare/M/PSYCH/KREMEN/PupillometryV3/data/raw/pupillometry_behavioral_data_01052017.csv'
outfile = '/home/jelman/netshare/M/PSYCH/KREMEN/PupillometryV3/data/DS_30Hz_VETSA3_01052017.xlsx'
###########################################################################

### Load Data ###

# Load parsed pupil data
pupildf = pd.read_csv(pupil_fname, sep=",")

if '0.000000000' in pupildf.columns:
    pupildf.rename(columns={'0.000000000':'0'}, inplace=True)

# Load behavioral data with trial times and scan quality
behavdf = pd.read_csv(behav_fname, sep=",")
behavdf.columns = behavdf.columns.str.replace("_v2","")

if 'DSTIM' in behavdf.columns:
    behavdf = behavdf.drop("DSTIM", 1)

if 'SUBJECTID' in behavdf.columns:
    behavdf.rename(columns={"SUBJECTID":"vetsaid"}, inplace=True)

### Behavioral Data ###

# Select columns of interest
dstcols = [col for col in behavdf.columns if ("DST" in col) & (("SCN" in col) or ("TIM" in col))]
behavdf = behavdf[["vetsaid"]+dstcols]

# Rotate variable names to make it easier to go from wide format to long
behavdf.columns = [rotate(col, -3) if col in dstcols else col for col in behavdf.columns]
# Convert from wide to long format
behavdflong = pd.wide_to_long(behavdf, stubnames=["TIM","SCN"], i="vetsaid", j="Trial #").reset_index()
behavdflong['TIM'].replace("999999",np.nan, inplace=True)

# Strip string from trial number and convert to int so it will sort in numerical order
behavdflong["Trial #"] = behavdflong["Trial #"].str.replace("DST","").astype(int)
# Sort by subject id
behavdflong = behavdflong.sort_values(by=["vetsaid","Trial #"])
# Rename subject ID field
behavdflong = behavdflong.rename(columns={"vetsaid":"ID"})

# Add Pair ID field
behavdflong['Pair ID'] = behavdflong['ID'].str.replace("[AB]", "")

# Get rid of decimal and convert to string
behavdflong['TIM'] = behavdflong['TIM'].astype(str).str.replace("\.0","").str.zfill(6)

# Convert to timestamp
behavdflong.ix[behavdflong['TIM'].str.contains("nan"),'TIM'] = np.nan

behavdflong['Pupil Trial Time'] = behavdflong['TIM'].str[:2] + ":" + behavdflong['TIM'].str[2:4] + ":" + behavdflong['TIM'].str[4:]

# Create Task field
behavdflong['Task (DS)'] = 'DS' + behavdflong['Trial #'].astype(str)
bxtrialnames = [''.join([i,j,k]) for i in ['S'] for j in ['3','6','9'] for k in ['a','b','c','d']]
# Make sure all subjects have 12 rows
assert (behavdflong.groupby('ID')['ID'].count()==12).all()
behavdflong['Bx Trial'] = bxtrialnames * behavdflong.ID.nunique()

# Specify bad data trials and filter out unacquired trials
behavdflong = behavdflong.dropna(axis=0, subset=["Pupil Trial Time"])
behavdflong.loc[behavdflong['SCN']<>1,'Pupil Trial Time'] = "bad data"

### Pupil Data ###

# Calculate Min and Max fields
# Insert min
pupildf.insert(pupildf.columns.get_loc("0"), 'Min', pupildf.ix[:,'0':].min(axis=1))
# Insert max
pupildf.insert(pupildf.columns.get_loc("0"), 'Max', pupildf.ix[:,'0':].max(axis=1))

# Rename columns
pupildf = pupildf.rename(columns={"Subject ID":"ID", "Eye Measured":"LR", "Time":"Pupil Trial Time"})
ndatacols = pupildf.columns[pupildf.columns.get_loc("0"):].shape[0]
datacolnames = [str(x +1) + 'st data pt' for x in range(ndatacols)]
pupildf.columns = list(pupildf.columns[:pupildf.columns.get_loc("0")]) + datacolnames


### Merge data ###
mergeddf = behavdflong.merge(pupildf, how="inner", on=["ID","Pupil Trial Time"])

# Get number of twins (paired or singleton)
ntwins = mergeddf.groupby(['Pair ID']).apply(lambda x : x.ID.nunique())
mergeddf['Matched?'] = ""

templateDS = mergeddf[['ID','LR','Pair ID','Matched?','Task (DS)',
                       'Pupil Trial Time','Trial #', 'Bx Trial','Min',
                       'Max']+datacolnames]
templateDS.to_excel(outfile, index=False)
