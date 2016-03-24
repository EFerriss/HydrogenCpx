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

lengths_microns = [3000.] * 3
logD_list = [-12.5, -12., -11.5, -11., -10.5, -10.]

#time_minutes = [0.25, 0.5, 0.75] + list(np.arange(1, 11, 0.5)) + range(15, 65, 1) + range(65, 5*65, 5)
time_minutes = list(np.linspace(0.00001, 5., 100) * 60.)

direction = 0 # 0=[100]*, 1=[010], 2=[001]
points_in_calc = 100

#%%################# Set up and run calculation ################################
fig, ax, v, x, y = diff.diffusion3Dnpi(lengths_microns, [-20.]*3, 
                              time_seconds=0.000000000001,
                              plot3=True, points=points_in_calc,
                              centered=False)
v_sat = np.sum(v) 

data = [time_minutes]
#data = []

for D_m2s in logD_list:
    print D_m2s
    D3 = [D_m2s]*3

    percent_equilibrated = np.zeros_like(time_minutes)
    rim_location_microns = np.zeros_like(time_minutes)

    for idx, minutes in enumerate(time_minutes):
        time_seconds =  minutes * 60.
    
        v, x, y = diff.diffusion3Dnpi(lengths_microns, D3, time_seconds,
                                      plot3=False, points=points_in_calc)
        v_sum = np.sum(v)
        percent_equilibrated[idx] = 100. * (v_sat-v_sum) / v_sat

        print ''.join(('logD=', '{:.1f}'.format(D_m2s), ', ',
                       '{:.1f}'.format(minutes), ' minutes done'))
                       
        if v_sum <= 0:
            print '   water gone!'
            percent_equilibrated[idx+1:-1] = 100.
            break

    data.append(list(percent_equilibrated))

print data
data.append(logD_list)

#%% Save data to file
workfile = 'xenolithEquilibration-3mm.txt'
with open(workfile, 'w') as diff_file:
    diff_file.write(json.dumps(data))

#%% Plotting
fig = plt.figure()
fig.set_size_inches(6, 5)
ax = fig.add_subplot(111)
plt.style.use('paper')

colors = ['green', 'b', 'purple', 'red', 'black']
time_hours = np.array(time_minutes) / 60.

for idx, percentRemaining in enumerate(data[1:-1]):
    ax.plot(time_hours, percentRemaining, '-', mew=1, linewidth=2, )
#            color=colors[idx_D],
#            label=''.join(('logD=', '{:.1f}'.format(logD_list[idx]))

volume_mm3 = (lengths_microns[0] * lengths_microns[1] * lengths_microns[2]) / 1E9
tit = ''.join(('{:.1f}'.format(volume_mm3), ' mm$^3$cube'))

#ylab = 'Preservation of hydrogen (%)'
ylab = '% equilibration'
ax.set_ylabel(ylab)
ax.set_xlabel('Time (hours)')
ax.set_title(tit)
ax.set_ylim(0., 100.)
ax.set_xlim(0, max(time_hours))

