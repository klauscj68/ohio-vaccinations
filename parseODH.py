# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 23:07:24 2021

Parse the ODH master data sheet into daily incidence by age cohort

@author: Colin
"""
#%% Import needed packages
import pandas as pd
import numpy as np
from datetime import datetime
import os

#%% Initialize
mydir = 'ODH_Data'
DAY = '0616'

#%% Download and save ODH data
DF_ODH = pd.read_csv('https://coronavirus.ohio.gov/static/dashboards/COVIDDeathData_CountyOfResidence.csv')
DF_ODH.to_csv('ODHraw_' + DAY + '.csv')
os.rename('ODHraw_'+DAY+'.csv',mydir + '/ODHraw_'+DAY+'.csv')
DF_VAX = pd.read_csv('https://coronavirus.ohio.gov/static/dashboards/vaccine_data.csv')
DF_VAX.to_csv('VAXraw_' + DAY + '.csv')
os.rename('VAXraw_'+DAY+'.csv',mydir + '/VAXraw_'+DAY+'.csv')

#%% Date formatter
def date_fmt(day):
    """
    Convert from ODH date format to np.datetime64 format and give
    integer index for 2020-01-01 is day 0

    Parameters
    ----------
    day : String in format M/D/YYYY

    Returns
    -------
    day : String in format YYYY-MM-DD

    """
    # # ODH switched their case reports to match the YYYY-MM-DD format
    # ram = day.split('/')
    # if len(ram[0]) == 1:
    #     ram[0] = '0'+ram[0]
    # if len(ram[1]) == 1:
    #     ram[1] = '0'+ram[1]
    # day = ram[2]+'-'+ram[0]+'-'+ram[1]
    day = np.datetime64(day)

    tdelta = (day - np.datetime64('2020-01-01')).astype(int)
    return day, tdelta

#%% Parse for dates and ages for incidence
AGES = []
DATES = []

for i in range(DF_ODH.shape[0]):
    # Buld age keys
    flag_new = True
    ram = DF_ODH['Age Range'][i]
    if ram == 'Total':
        flag_new = False
    else:
        for j in range(len(AGES)):
            if ram == AGES[j]:
                flag_new = False
                break
    if flag_new:
        AGES.append(ram)

    # Build date keys
    flag_new = True
    ram = DF_ODH['Onset Date'][i]
    if (ram == 'Total')|(str(ram) == 'nan')|(len(str(ram)) == 0):
        flag_new = False
    else:
        for j in range(len(DATES)):
            if ram == DATES[j]:
                flag_new = False
                break
    if flag_new:
        DATES.append(ram)

# Convert dates to np datetime format and numeric type
DATES64 = [date_fmt(x)[0] for x in DATES]
taxis = np.array([date_fmt(x)[1] for x in DATES])


#%% Create master array storing incidence data
MST = np.zeros((taxis.max(),9))
AGEDICT = {'0-19':[0,1],
           '20-29':[2],
           '30-39':[3],
           '40-49':[4],
           '50-59':[5],
           '60-69':[6],
           '70-79':[7],
           '80+':[8]}

# To count if ODH changed its file format and you missed anything
mycheck = 0
for i in range(DF_ODH.shape[0]):
    # Find the age columns
    agenow = AGEDICT[DF_ODH['Age Range'][i]]

    ram = DF_ODH['Onset Date'][i]
    if (ram == 'Total')|(str(ram) == 'nan')|(len(str(ram)) == 0):
        mycheck += int(DF_ODH['Case Count'][i])
        continue

    # Find the date rel Jan 1, 2020
    tnow = date_fmt(ram)[1]

    new = float(DF_ODH['Case Count'][i])
    if len(agenow) == 2:
        MST[tnow-1,0:2] += np.array([0.49171424240893247*new,
                                     0.5082857575910675*new])
    else:
        MST[tnow-1,agenow[0]] += new

if mycheck != 0:
    print('Warning: ' + str(mycheck)
          + ' many cases did not have dates ...')

# Add two extra columns for Daily All and Cum All
daily_all = np.sum(MST,axis=1)
cum_all = np.cumsum(daily_all)

MST = np.hstack((MST,
                 daily_all.reshape((-1,1)),
                 cum_all.reshape((-1,1))))

#%% Turn into data frame
COLS = ['0-9','10-19','20-29','30-39','40-49',
        '50-59','60-69','70-79','80+',
        'daily_confirm', 'cum_confirm']
DAY0 = np.datetime64('2020-01-01')
t0 = taxis.min()
tfin = taxis.max()
taxis = [(DAY0 + np.timedelta64(i,'D')).astype(str)
         for i in range(t0,tfin+1)]

DF = pd.DataFrame(MST,index=taxis,columns=COLS)
DF.index.name = 'time'
DF.to_csv('ODH_'+DAY+'.csv')
os.rename('ODH_'+DAY+'.csv',mydir+'/ODH_'+DAY+'.csv')

#%% Parse for dates and ages for deaths
AGES = []
DATES = []

for i in range(DF_ODH.shape[0]):
    # Buld age keys
    flag_new = True
    ram = DF_ODH['Age Range'][i]
    if ram == 'Total':
        flag_new = False
    else:
        for j in range(len(AGES)):
            if ram == AGES[j]:
                flag_new = False
                break
    if flag_new:
        AGES.append(ram)

    # Build date keys
    flag_new = True
    ram = DF_ODH['Date Of Death'][i]
    if (ram == 'Total')|(str(ram) == 'nan')|(len(str(ram)) == 0):
        flag_new = False
    else:
        for j in range(len(DATES)):
            if ram == DATES[j]:
                flag_new = False
                break
    if flag_new:
        DATES.append(ram)

# Convert dates to np datetime format and numeric type
DATES64 = [date_fmt(x)[0] for x in DATES]
taxis = np.array([date_fmt(x)[1] for x in DATES])


#%% Create master array storing incidence data
MST = np.zeros((taxis.max(),9))
AGEDICT = {'0-19':[0,1],
           '20-29':[2],
           '30-39':[3],
           '40-49':[4],
           '50-59':[5],
           '60-69':[6],
           '70-79':[7],
           '80+':[8]}

# To count if ODH changed its file format and you missed anything
mycheck = 0
for i in range(DF_ODH.shape[0]):
    # Find the age columns
    agenow = AGEDICT[DF_ODH['Age Range'][i]]

    ram = DF_ODH['Date Of Death'][i]
    if (ram == 'Total')|(str(ram) == 'nan')|(len(str(ram)) == 0):
        mycheck += int(DF_ODH['Death Due To Illness Count - County Of Residence'][i])
        continue

    # Find the date rel Jan 1, 2020
    tnow = date_fmt(ram)[1]

    new = float(DF_ODH['Death Due To Illness Count - County Of Residence'][i])
    if len(agenow) == 2:
        MST[tnow-1,0:2] += np.array([0.49171424240893247*new,
                                     0.5082857575910675*new])
    else:
        MST[tnow-1,agenow[0]] += new

if mycheck != 0:
    print('Warning: ' + str(mycheck)
          + ' many cases did not have dates ...')

# Add two extra columns for Daily All and Cum All
daily_all = np.sum(MST,axis=1)
cum_all = np.cumsum(daily_all)

MST = np.hstack((MST,
                 daily_all.reshape((-1,1)),
                 cum_all.reshape((-1,1))))

#%% Turn into data frame
COLS = ['0-9','10-19','20-29','30-39','40-49',
        '50-59','60-69','70-79','80+',
        'daily_confirm', 'cum_confirm']
DAY0 = np.datetime64('2020-01-01')
t0 = taxis.min()
tfin = taxis.max()
taxis = [(DAY0 + np.timedelta64(i,'D')).astype(str)
         for i in range(t0,tfin+1)]

# Truncate MST to run from tmin to tfin
MST = MST[(t0-1):,:]

DF = pd.DataFrame(MST,index=taxis,columns=COLS)
DF.index.name = 'time'
DF.to_csv('ODH_'+DAY+'_death.csv')
os.rename('ODH_'+DAY+'_death.csv',mydir+'/ODH_'+DAY+'_death.csv')

#%% Parse for dates for vaccinations
DATES = []

for i in range(DF_VAX.shape[0]):
    # Build date keys
    flag_new = True
    ram = DF_VAX['date'][i]
    if (ram == 'Total')|(str(ram) == 'nan')|(len(str(ram)) == 0):
        flag_new = False
    else:
        for j in range(len(DATES)):
            if ram == DATES[j]:
                flag_new = False
                break
    if flag_new:
        DATES.append(ram)

# Convert dates to np datetime format and numeric type
DATES64 = [date_fmt(x)[0] for x in DATES]
taxis = np.array([date_fmt(x)[1] for x in DATES])


#%% Create master array storing vaccination data
MST = np.zeros((taxis.max(),2))

# To count if ODH changed its file format and you missed anything
# mycheck = 0
for i in range(DF_VAX.shape[0]):

    ram = DF_VAX['date'][i]
    # if (ram == 'Total')|(str(ram) == 'nan')|(len(str(ram)) == 0):
    #    mycheck += int(DF_ODH['Case Count'][i])
    #    continue

    # Find the date rel Jan 1, 2020
    tnow = date_fmt(ram)[1]

    new = float(DF_VAX['vaccines_started'][i])
    MST[tnow-1,0] += new

    new = float(DF_VAX['vaccines_completed'][i])
    MST[tnow-1,1] += new

#if mycheck != 0:
#    print('Warning: ' + str(mycheck)
#          + ' many cases did not have dates ...')

# Add two extra columns for cumulative
cum_start = np.cumsum(MST[:,0])
cum_finish = np.cumsum(MST[:,1])

MST = np.hstack((MST,
                 cum_start.reshape((-1,1)),
                 cum_finish.reshape((-1,1))))

#%% Turn into data frame
COLS = ['vax_start','vax_finish',
        'cum_start', 'cum_finish']
DAY0 = np.datetime64('2020-01-01')
t0 = 1
tfin = taxis.max()
taxis = [(DAY0 + np.timedelta64(i,'D')).astype(str)
         for i in range(t0,tfin+1)]

DF = pd.DataFrame(MST,index=taxis,columns=COLS)
DF.index.name = 'time'
DF.to_csv('ODH_'+DAY+'_vax.csv')
os.rename('ODH_'+DAY+'_vax.csv',mydir+'/ODH_'+DAY+'_vax.csv')
