# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 15:14:02 2015

@author: Ferriss

Figure roughly following Padron-Navarta et al. 2014 Figure 8
Comparing olivine and cpx rim formation over time

Generate data and save to xenolithRims.txt in json format
"""
import pynams.diffusion as diff
import numpy as np
import json

################### User input variables #####################################
savefolder = '../../../CpxPaper/figures/'

lengths_microns = [1000.] * 3
logD_list = [-12.5, -12., -11.5, -11., -10.5, -10.]
#logD_list = [-10.]

time_minutes = [0.25, 0.5, 0.75] + list(np.arange(1, 11, 0.5)) + range(15, 65, 1) + range(65, 5*65, 5)

water_fraction = 0.90 # How much does concentration drop to qualify as a rim

direction = 0 # 0=[100]*, 1=[010], 2=[001]
points_in_calc = 30

################## Set up and run calculation ################################
v_sat = np.sum(np.ones([points_in_calc, points_in_calc, points_in_calc]))
center_microns = lengths_microns[0] / 2.
savedata = [time_minutes]

for D_m2s in logD_list:
    D3 = [D_m2s]*3

    percent_water_remaining = np.zeros_like(time_minutes)
    rim_location_microns = np.zeros_like(time_minutes)

    idx = 0
    for minutes in time_minutes:
        time_seconds =  minutes * 60.
    
        v, x, y = diff.diffusion3Dnpi(lengths_microns, D3, time_seconds,
                                      plot3=False, points=points_in_calc)
        percent_water_remaining[idx] = 100. * (v_sat - np.sum(v)) / v_sat
        idx_rim = (np.abs(y[direction][0:points_in_calc/2] - water_fraction)).argmin()
        rim_location_microns[idx] = x[direction][idx_rim]

        print ''.join(('logD=', '{:.1f}'.format(D_m2s), ', ',
                       '{:.1f}'.format(minutes), ' minutes done, ',
                       '{:.1f}'.format(rim_location_microns[idx]), ' microns',))

        if (center_microns - rim_location_microns[idx]) < 20.:
            print 'hit center!'
            rim_location_microns[idx+1:] = center_microns
#            print rim_location_microns
            break
        
        idx = idx + 1

    savedata.append(list(rim_location_microns))


savedata.append(logD_list)

# Save data to file
workfile = ''.join((savefolder, '-xenolithRims_90.txt'))
with open(workfile, 'w') as diff_file:
    diff_file.write(json.dumps(savedata))

