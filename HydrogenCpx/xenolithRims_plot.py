# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 15:14:02 2015

@author: Ferriss

Figure roughly following Padron-Navarta et al. 2014 Figure 8
Comparing olivine and cpx rim formation over time

Data generated and saved to xenolithRims.txt in json format 
in xenolithRims_data.py
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import json

water_fraction = 0.80 # How much does concentration drop to qualify as a rim
resolution_microns = 100.

lengths_microns = [1000.] * 3

savefolder = '../../../../CpxPaper/figures/'
workfile_80 = ''.join((savefolder, '-xenolithRims.txt'))
workfile_90 = ''.join((savefolder, '-xenolithRims_90.txt'))

with open(workfile_80, 'r') as rimfile:
    a = rimfile.read()
data = json.loads(a)

with open(workfile_90, 'r') as rimfile:
    a = rimfile.read()
data_90 = json.loads(a)

time_minutes = data[0]
time_minutes_90 = data_90[0]
logD_list = data[-1]
logD_list_90 = data_90[-1]

center_microns = lengths_microns[0] / 2.
volume_mm3 = (lengths_microns[0] * lengths_microns[1] * lengths_microns[2]) / 1E9

# delete logD = -10.5 in data
idx_to_delete = -2
data = np.delete(data, idx_to_delete, 0)
logD_list = np.delete(logD_list, idx_to_delete)

#%% Plotting
time_hours = np.array(time_minutes) / 60.
fig = plt.figure()
fig.set_size_inches(6, 4)
ax = fig.add_subplot(111)
plt.style.use('paper')

colors = ['blue', 'green', 'red', 'teal', 'black']
fillcolors = ['blue', 'green', 'red', 'teal', 'yellow']
#colors = ['blue']*5 + ['black']
#fillcolors = ['lightblue']*4 + ['yellow']
hatches = [None, None, None, '\\', '*']

idx_D = 0
for rimloss in data[1:-1]:
    ax.plot(time_hours, rimloss, '-', mew=1, linewidth=2, color=colors[idx_D],
            label=''.join(('logD=', '{:.1f}'.format(logD_list[idx_D]))))
    idx_D = idx_D + 1

idx_D = 0
for rimloss in data_90[1:-1]:
    ax.plot(time_hours, rimloss, '-', mew=1, linewidth=2, color=colors[idx_D],
            label=''.join(('logD=', '{:.1f}'.format(logD_list_90[idx_D]))))
    idx_D = idx_D + 1

for idx in range(0, 5):
    plt.fill_between(time_hours, data[idx+1], data_90[idx+1],
                     facecolor=fillcolors[idx], alpha=0.4, interpolate=True,
                     hatch=hatches[idx])

#ylab = ''.join(('Distance from rim in center of cube face where\nhydrogen concentration <', 
#                '{:.0f}'.format(water_fraction*100), '% of initial ($\mu$m)'))
ylab = ''.join(('Distance from center of cube face where\nH concentration <', 
                '80-90% of initial ($\mu$m)'))

ax.set_xlabel('Time (hours)')
tit_volume = ''.join(('{:.1f}'.format(volume_mm3), ' mm$^3$cube'))
tit = ''.join(('A. ', tit_volume))

ax.set_ylabel(ylab)
ax.set_title(tit)
ax.set_ylim(center_microns, 0)
ax.set_xlim(0, 5.)

plt.axhspan(0, resolution_microns, facecolor='green', alpha=0.2)
ax.text(ax.get_xlim()[1]/2., resolution_microns/2., 
        'Diffusion profiles unlikely to be observed', ha='center', va='center')
        
bottom_box_top = center_microns - 0.06*center_microns
plt.axhspan(center_microns, bottom_box_top, 
            facecolor='orange', alpha=0.5)
#yt = center_microns - 0.07*center_microns
yt = 485.
ax.text(290 / 60., yt, 'hydrogen loss at the center', ha='right', va='center')

#ax.legend(loc=4)
labels = []
for idx in range(5):
    labels.append(''.join(('logD=', '{:.1f}'.format(logD_list[idx]), ' m$^2$','/s')))

ax.text(150 / 60., 120, labels[0], rotation=-10)
ax.text(140 / 60., 200, labels[1], rotation=-20)
ax.text(85 / 60., 300, labels[2], rotation=-60)
ax.text(30 / 60., 270, labels[3], rotation=-80)
ax.text(9 / 60., 300, labels[4], rotation=-87)

for y in [resolution_microns, bottom_box_top]:
    ax.plot(ax.get_xlim(), [y, y], '-k', linewidth=3)

### time ranges for comparison
#depth_km = np.array([40., 50.])
ascent_rate_m_per_s = [0.2, 0.5]
depth_km = np.array([2., 3.])
#ascent_rate_m_per_s = [15., 20., 30., 50., 100.]
#ytloc = [320, 175, 300]
ytloc = [450] * len(ascent_rate_m_per_s)

depth_m = depth_km * 1E3
idx = 0
for rate in ascent_rate_m_per_s:
    ascent_time_s = depth_m / rate
    ascent_time_m = ascent_time_s / 60.

#    ax.axvspan(ascent_time_m[0], ascent_time_m[1], facecolor='b', alpha=0.05)
    rect = patches.Rectangle((ascent_time_m[0] / 60., resolution_microns), 
                             ascent_time_m[1]  / 60.- ascent_time_m[0] / 60., 
                             bottom_box_top - resolution_microns, 
                             color='grey', alpha=0.2)
    ax.add_patch(rect)
    ax.text(np.mean(ascent_time_m) / 60., ytloc[idx], 
        ''.join(('{:.1f}'.format(rate), ' m/s\n(',
                 '{:.0f}'.format(depth_km[0]), '-', 
                 '{:.0f}'.format(depth_km[1]), ' km)')),
        rotation=-90., ha='center', va='bottom')
    idx = idx + 1
#####

#for minutes in [60.]:
#    ax.plot([minutes, minutes], ax.get_ylim(), '-r')


#plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig15A.eps', 
#            format='eps', dpi=1000)
#plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig15A.tif', dpi=600)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig15A.tif')

plt.show(fig)
print 'Finished'
