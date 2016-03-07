# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 15:14:02 2015

@author: Ferriss

Figure roughly following Padron-Navarta et al. 2014 Figure 8
Comparing olivine and cpx rim formation over time

Data generated and saved to xenolithRims.txt in json format 
in xenolithRims_data.py

Updated for 3mm sided cube instead of 1mm 
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import json

water_fraction = 0.80 # How much does concentration drop to qualify as a rim
resolution_microns = 100.

lengths_microns = [3000.] * 3

savefolder = '\\Users\\Ferriss\\Documents\\Code\\Python\\'
workfile_80 = ''.join((savefolder, '-xenolithRims-3mm-80.txt'))
workfile_90 = ''.join((savefolder, '-xenolithRims-3mm-90.txt'))

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
#idx_to_delete = -2
#data = np.delete(data, idx_to_delete, 0)
#logD_list = np.delete(logD_list, idx_to_delete)

#%% Plotting

fig = plt.figure()
fig.set_size_inches(6, 4)
ax = fig.add_subplot(111)
plt.style.use('paper')

colors = ['red', 'teal', 'black', 'blue', 'green']
fillcolors = ['blue', 'green', 'red', 'teal', 'yellow']
hatches = [None, None, None, '\\', '*']

idx_D = 0
for rimloss in data[3:-1]:
    ax.plot(time_minutes, rimloss, '-', mew=1, linewidth=2, color=colors[idx_D],
            label=''.join(('logD=', '{:.1f}'.format(logD_list[idx_D]))))
    idx_D = idx_D + 1

idx_D = 0
for rimloss in data_90[3:-1]:
    ax.plot(time_minutes, rimloss, '-', mew=1, linewidth=2, color=colors[idx_D],
            label=''.join(('logD=', '{:.1f}'.format(logD_list_90[idx_D]))))
    idx_D = idx_D + 1

for idx in range(2, 5):
    plt.fill_between(time_minutes, data[idx+1], data_90[idx+1],
                     facecolor=fillcolors[idx], alpha=0.4, interpolate=True,
                     hatch=hatches[idx])

#ylab = ''.join(('Distance from rim in center of cube face where\nhydrogen concentration <', 
#                '{:.0f}'.format(water_fraction*100), '% of initial ($\mu$m)'))
ylab = ''.join(('Distance from center of cube face where\nH concentration <', 
                '80-90% of initial ($\mu$m)'))

ax.set_xlabel('Time (minutes)')
tit_volume = ''.join(('{:.1f}'.format(volume_mm3), ' mm$^3$cube'))
tit = ''.join(('B. ', tit_volume, ', short timescales'))
#tit = tit_volume

ax.set_ylabel(ylab)
ax.set_title(tit)
ax.set_ylim(center_microns, 0)
#ax.set_xlim(0, 30.)

plt.axhspan(0, resolution_microns, facecolor='green', alpha=0.2)
ax.text(ax.get_xlim()[1]/2., resolution_microns/2., 
        'Diffusion profiles unlikely to be observed', ha='center', va='center')
        
bottom_box_top = center_microns - 0.06*center_microns
plt.axhspan(center_microns, bottom_box_top, 
            facecolor='orange', alpha=0.5)
#yt = center_microns - 0.07*center_microns
yt = 485.
ax.text(29, yt, 'hydrogen loss at the center', ha='right', va='center')

#ax.legend(loc=4)
labels = []
for idx in range(5):
    labels.append(''.join(('logD=', '{:.1f}'.format(logD_list[idx]), ' m$^2$','/s')))

#labels[3] = ''.join(('~clinopyroxene, ', labels[3]))
#ax.text(150, 120, labels[0], rotation=-10)
#ax.text(140, 200, labels[1], rotation=-20)
ax.text(20, 135, labels[2], rotation=-8)
ax.text(20, 240, labels[3], rotation=-17)
ax.text(10.8, 300, labels[4], rotation=-64)
ax.text(6, 190, '~olivine', rotation=-35)
ax.text(20, 175, '~clinopyroxene', rotation=-9)

for y in [resolution_microns, bottom_box_top]:
    ax.plot(ax.get_xlim(), [y, y], '-k', linewidth=3)

### time ranges for comparison
#ascent_rate_m_per_s = [10., 2.]
depth_km = np.array([2., 3.])
#depth_km = np.array([30., 35.])
ascent_rate_m_per_s = [10., 2.]
#ax.set_xlim(0, 10)
ytloc = [450] * len(ascent_rate_m_per_s)

depth_m = depth_km * 1E3
idx = 0
for rate in ascent_rate_m_per_s:
    ascent_time_s = depth_m / rate
    ascent_time_m = ascent_time_s / 60.

#    ax.axvspan(ascent_time_m[0], ascent_time_m[1], facecolor='b', alpha=0.05)
    rect = patches.Rectangle((ascent_time_m[0], resolution_microns), 
                             ascent_time_m[1] - ascent_time_m[0], 
                             bottom_box_top - resolution_microns, 
                             color='grey', alpha=0.2)
    ax.add_patch(rect)
#    ax.text(np.mean(ascent_time_m), ytloc[idx], 
#        ''.join(('{:.1f}'.format(rate), ' m/s (',
#                 '{:.0f}'.format(depth_km[0]), '-', 
#                 '{:.0f}'.format(depth_km[1]), ' km)')),
#        rotation=-90., ha='center', va='bottom')
    idx = idx + 1
#####
ax.text(21., 450., '2 m/s (2-3 km)', rotation=-90., ha='center', va='bottom')
ax.text(2., 450., '$\uparrow$ 10 m/s\n(2-3 km)', rotation=-90., ha='center', va='bottom')

#for minutes in [60.]:
#    ax.plot([minutes, minutes], ax.get_ylim(), '-r')

#plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig15B.eps', 
#            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig15B-3mm.tif')

plt.show(fig)
print 'Finished'

