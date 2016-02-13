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
import matplotlib.pyplot as plt
import json

################### User input variables #####################################
savefolder = 'C://Users//Ferriss//Documents//CpxPaper//figures//'

lengths_microns = [2000.] * 3
logD_list = [-15., -14., -13., -12., -11., -10.]
time_minutes = np.linspace(0.00001, 744., 75) * 60.
time_hours = time_minutes / 60.

direction = 0 # 0=[100]*, 1=[010], 2=[001]
points_in_calc = 25

################## Set up and run calculation ################################

# initial concentration summed over all cells
#v_sat = np.sum(np.ones([points_in_calc, points_in_calc, points_in_calc]))
v, x, y = diff.diffusion3Dnpi(lengths_microns, [-20.]*3, 
                              time_seconds=0.0000000001,
                              plot3=False, points=points_in_calc)
v_sat = np.sum(v) 

# looping through each diffusivity and time                                 
data = [list(time_minutes)]
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

        print ''.join(('logD=', '{:.1f}'.format(D_m2s), ', ',
                       '{:.1f}'.format(minutes), ' minutes done',
                       ' % remaining ', '{:.1f}'.format(percent_water_remaining[idx])))
        idx = idx + 1

    data.append(list(percent_water_remaining))
print data


#%%### Save data to file
data.append(logD_list)

workfile = ''.join((savefolder, '-xenolithEquilibration_2mm_logD14.txt'))
with open(workfile, 'w') as diff_file:
    diff_file.write(json.dumps(list((data))))
