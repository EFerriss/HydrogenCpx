# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 15:14:02 2015

@author: Ferriss

Figure roughly following Padron-Navarta et al. 2014 Figure 8
Comparing olivine and cpx rim formation over time

for cube with 1 mm sides

Generate data and save to xenolithRims.txt in json format
"""
import pynams.diffusion as diff
import numpy as np
import json

################### User input variables #####################################

lengths_microns = [1000.] * 3
logD_list = [-12.5, -12., -11.5, -11., -10.5, -10.]

# 0 to 24 hours
time_minutes = [0.25, 0.5, 0.75] + list(np.arange(1, 11, 0.5)) + range(15, 65, 1) + range(65, 5*65, 5)

water_fraction = 0.80 # How much does concentration drop to qualify as a rim

direction = 0 # 0=[100]*, 1=[010], 2=[001]
points_in_calc = 300

################## Set up and run calculation ################################
fig, ax, v, x, y = diff.diffusion3Dnpi(lengths_microns, [-20.]*3, 
                              time_seconds=0.000000000001,
                              plot3=True, points=points_in_calc,
                              centered=False)
v_sat = np.sum(v) 

#%%

savedata = [time_minutes]

for D_m2s in logD_list:
    D3 = [D_m2s]*3

    percent_water_remaining = np.zeros_like(time_minutes)
    rim_location_microns = np.zeros_like(time_minutes)

    for idx, minutes in enumerate(time_minutes):
        time_seconds =  minutes * 60.
    
        # determine internal concentration (v), positions through the 
        # center (x), and concentration profiles through the center (y)
        # for each time step        
        v, x, y = diff.diffusion3Dnpi(lengths_microns, D3, time_seconds,
                                      plot3=False, points=points_in_calc,
                                      centered=False)

        # check percent water remaining
        v_sum = np.sum(v)
        
        # should ever actually get to this point
        if v_sum <= 0:
            print '   water gone!'
            break
            
        percent_water_remaining[idx] = 100. * (v_sat - v_sum) / v_sat

        # find the rim
        idx_rim = (np.abs(y[direction][0:points_in_calc/2] - water_fraction)).argmin()
        rim_location_microns[idx] = x[direction][idx_rim]

        print ''.join(('logD=', '{:.1f}'.format(D_m2s), ', ',
                       '{:.1f}'.format(minutes), ' minutes, rim at ', 
                       '{:.1f}'.format(rim_location_microns[idx]) ))

        if idx_rim == (points_in_calc/2)-1:
            print '   diffusion hit center!'
            rim_location_microns[idx+1:-1] = rim_location_microns[idx]
            break

    savedata.append(list(rim_location_microns))


savedata.append(logD_list)

#%% Save data to file
workfile = 'xenolithRims-1mm-80.txt'
with open(workfile, 'w') as diff_file:
    diff_file.write(json.dumps(savedata))

