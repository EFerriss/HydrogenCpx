# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 15:14:02 2015

@author: Ferriss

Figure roughly following Padron-Navarta et al. 2014 Figure 8
Comparing olivine and cpx rim formation over time
"""
import pynams.diffusion as diff
import numpy as np
import matplotlib.pyplot as plt
import json

################### User input variables #####################################
savefolder = '../../../CpxPaper/figures/'

celsius = 1100.

lengths_microns = [1000.] * 3

cpx_log10D_m2s = [-12.] * 3
mid_log10D_m2s = [-11.] * 3
#oli_log10D_m2s = np.log10(np.array([3.1E-10, 1.8E-11, 1.0E-11])) # KM98 1100C fast mechanism
oli_log10D_m2s = [-10.] * 3

#time_hours = [0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05]
time_hours = [0.005, 0.01, 0.03, 0.05, 0.1, 0.175, 0.25, 0.375, 0.5, 0.75, 1.]
#time_hours = [0.005, 0.05, 0.5, 1., 4., 8.]
#time_hours = [0.05, 0.5, 1., 2., 4., 8., 10., 15., 20.]
#time_hours = [0.5]

water_fraction = 0.9 # How much does concentration drop to qualify as a rim
resolution_microns = 100.

direction = 0 # 0=[100]*, 1=[010], 2=[001]
points_in_calc = 50

################## Set up and run calculation ################################
v_sat = np.sum(np.ones([points_in_calc, points_in_calc, points_in_calc]))
volume_mm3 = lengths_microns[0] * lengths_microns[1] * lengths_microns[2] / 1E9
center_microns = lengths_microns[direction] / 2.

cpx_percent_water_remaining = np.zeros_like(time_hours)
cpx_rim_location_microns = np.zeros_like(time_hours)
mid_percent_water_remaining = np.zeros_like(time_hours)
mid_rim_location_microns = np.zeros_like(time_hours)
oli_percent_water_remaining = np.zeros_like(time_hours)
oli_rim_location_microns = np.zeros_like(time_hours)

idx = 0
for hours in time_hours:
    time_seconds = hours * 3600.

    cpx_v, cpx_x, cpx_y = diff.diffusion3Dnpi(lengths_microns, cpx_log10D_m2s, time_seconds,
                                              plot3=False, points=points_in_calc)
    cpx_percent_water_remaining[idx] = 100. * (v_sat - np.sum(cpx_v)) / v_sat
    cpx_idx_rim = (np.abs(cpx_y[direction][0:points_in_calc/2] - water_fraction)).argmin()
    cpx_rim_location_microns[idx] = cpx_x[direction][cpx_idx_rim]

    mid_v, mid_x, mid_y = diff.diffusion3Dnpi(lengths_microns, mid_log10D_m2s, time_seconds,
                                              plot3=False, points=points_in_calc)
    mid_percent_water_remaining[idx] = 100. * (v_sat - np.sum(mid_v)) / v_sat
    mid_idx_rim = (np.abs(mid_y[direction][0:points_in_calc/2] - water_fraction)).argmin()
    mid_rim_location_microns[idx] = mid_x[direction][mid_idx_rim]

    oli_v, oli_x, oli_y = diff.diffusion3Dnpi(lengths_microns, oli_log10D_m2s, time_seconds,
                                              plot3=False, points=points_in_calc)
    oli_percent_water_remaining[idx] = 100. * (v_sat - np.sum(oli_v)) / v_sat
    oli_idx_rim = (np.abs(oli_y[direction][0:points_in_calc/2] - water_fraction)).argmin()
    oli_rim_location_microns[idx] = oli_x[direction][oli_idx_rim]

    # print out cpx info
    string = ' '.join(('hours:',
                        '{:.2f}'.format(hours), 
                        'rim observed at',
                        '{:.1f}'.format(cpx_rim_location_microns[idx]), 
                        'microns;',
                        '{:.1f}'.format(cpx_percent_water_remaining[idx]),
                        '% H lost'))
    print string
    idx = idx + 1

#plot final outcome - olivine
#f, ax, v, x, y = diff.diffusion3Dnpi(lengths_microns, oli_log10D_m2s, time_seconds,
#                    plot3=True, points=points_in_calc)

#### Save data to file
a = [time_hours, list(cpx_rim_location_microns), list(mid_rim_location_microns),
     list(oli_rim_location_microns)]

workfile = ''.join((savefolder, '-xenolithRims.txt'))
with open(workfile, 'w') as diff_file:
    diff_file.write(json.dumps(a))


#%%############## Plot #########################
workfile = ''.join((savefolder, '-xenolithRims.txt'))

with open(workfile, 'r') as rimfile:
    a = rimfile.read()
data = json.loads(a)

time_hours = data[0]
cpx_rim_location_microns = data[1]
mid_rim_location_microns = data[2]
oli_rim_location_microns = data[3]

fig = plt.figure()
fig.set_size_inches(6, 5)
ax = fig.add_subplot(111)
plt.style.use('paper')


time = np.array(time_hours) * 60.
ax.set_xlabel('Time (minutes)')


ax.plot(time, cpx_rim_location_microns, '+', markersize=10, mew=3,
        label=''.join(('logD=', '{:.1f}'.format(cpx_log10D_m2s[direction]),
                       ' m$^2$/s\n~min for Fe-rich cpx')))

ax.plot(time, mid_rim_location_microns, 'x', markersize=8, mew=2,
        label=''.join(('logD=', '{:.1f}'.format(mid_log10D_m2s[direction]),
                       ' m$^2$/s\n~Fe-rich cpx')))

ax.plot(time, oli_rim_location_microns, 'o', markersize=10, mew=1,
        label=''.join(('logD=', '{:.1f}'.format(oli_log10D_m2s[direction]),
                       ' m$^2$/s\n~olivine proton-polaron mech.')))

ylab = ''.join(('Distance from rim in center of cube face where\nhydrogen concentration <', 
                '{:.0f}'.format(water_fraction*100), '% of initial ($\mu$m)'))

tit_volume = ''.join(('{:.1f}'.format(volume_mm3), ' mm$^3$ cube'))
tit_temp = ''.join((' ', '{:.0f}'.format(celsius), ' $\degree$C'))
tit = ', '.join((tit_volume, tit_temp))

ax.set_ylabel(ylab)
ax.set_title(tit)
ax.set_ylim(center_microns, 0)

plt.axhspan(0, resolution_microns, facecolor='red', alpha=0.2)
ax.text(ax.get_xlim()[1]/2., resolution_microns/2., 
        'Diffusion profiles unlikely to be observed', ha='center', va='center')
        
plt.axhspan(center_microns, center_microns - 0.05*center_microns, 
            facecolor='yellow', alpha=0.5)
ax.text(ax.get_xlim()[0], center_microns - 0.07*center_microns, 
        '   loss at the center', ha='left', va='center')

ax.legend(loc=4)

fig.savefig(''.join((savefolder, 'xenolithRims.png')), dpi=200)