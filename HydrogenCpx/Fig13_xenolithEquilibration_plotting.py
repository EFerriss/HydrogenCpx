# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 15:14:02 2015

@author: Ferriss

Figure roughly following Padron-Navarta et al. 2014 Figure 8
Comparing olivine and cpx rim formation over time

Data generated and saved to xenolithEquilibrium.txt in json format 
in xenolithEquilibrium_data*.py

See also xenolithRims_*.py

Two panels, for 1mm (cpx-y) and 3mm (olivine-y) sided cubes
"""
import matplotlib.pyplot as plt
import json
import numpy as np

length_microns = np.array([1000., 3000.])
center_microns = length_microns / 2.
volume_mm3 = (length_microns**3) / 1E9

savefolder = 'C://Users//Ferriss//Documents//Code//Python//HydrogenCpx//HydrogenCpx//'
workfile = ''.join((savefolder, 'xenolithEquilibration-1mm.txt'))
workfile3mm = ''.join((savefolder, 'xenolithEquilibration-3mm.txt'))

with open(workfile, 'r') as rimfile:
    a = rimfile.read()
data = json.loads(a)
time_minutes = data[0]
time_hours = np.array(time_minutes) / 60.

with open(workfile3mm, 'r') as rimfile:
    a3mm = rimfile.read()
data3mm = json.loads(a3mm)
time_minutes3mm = data3mm[0]
time_hours3mm = np.array(time_minutes3mm) / 60.

logD_list = [-12.5, -12., -11.5, -11., -10.5, -10.]


#%% Plotting
fig = plt.figure()
fig.set_size_inches(6, 5)
plt.style.use('paper')

ax_big = fig.add_subplot(111)
ax_big.spines['top'].set_color('none')
ax_big.spines['bottom'].set_color('none')
ax_big.spines['left'].set_color('none')
ax_big.spines['right'].set_color('none')
ax_big.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

ax1mm = fig.add_subplot(121)
ax3mm = fig.add_subplot(122) 

colors = ['black', 'blue', 'green', 'darkgoldenrod', 'red', 'teal', ]

for idx, percentRemaining in enumerate(data[1:-1]):
    ax1mm.plot(time_hours, percentRemaining, '-k', mew=1, linewidth=2, 
#               color=colors[idx], 
               label=''.join(('logD=', '{:.1f}'.format(logD_list[idx]))))

for idx, percentRemaining in enumerate(data3mm[1:-1]):
    ax3mm.plot(time_hours3mm, percentRemaining, '-k', mew=1, linewidth=2, 
#               color=colors[idx], 
               label=''.join(('logD=', '{:.1f}'.format(logD_list[idx]))))

# make cpx estimated range green
ax1mm.fill_between(time_hours, data[2], data[3],facecolor='green', 
                   alpha=0.5, interpolate=True, )

ax3mm.fill_between(time_hours3mm, data3mm[2], data3mm[3],facecolor='green', 
                   alpha=0.5, interpolate=True, )


ylab = 'Extent of "water"\nre-equilibration with surroundings (%)'
ax1mm.set_ylabel(ylab)
ax_big.set_xlabel('Time (hours)')
ax3mm.set_yticklabels('')
fig.autofmt_xdate()

for idx, ax in enumerate([ax1mm, ax3mm]):
    tit = ''.join(('{:.0f}'.format(volume_mm3[idx]), ' mm$^3$cube'))
    ax.set_title(tit)
    ax.set_ylim(0., 100.)
    ax.set_xlim(0, max(time_hours))

### time ranges for comparison
ascent_rate_m_per_s = [10., 5.]
depth_km = np.array([40., 50.])
ytloc = [2] * len(ascent_rate_m_per_s)

depth_m = depth_km * 1E3

## ascent times
for ax in [ax1mm, ax3mm]:
    for idx, rate in enumerate(ascent_rate_m_per_s):
        ascent_time_s = depth_m / rate
        ascent_time_m = ascent_time_s / 60.
        ascent_time_h = ascent_time_m / 60.
    
        ax.axvspan(ascent_time_h[0], ascent_time_h[1], facecolor='grey', alpha=0.2)
    
        ax.text(np.mean(ascent_time_h), ytloc[idx], 
            ''.join(('{:.1f}'.format(rate), ' m/s',
                     ' (',
                     '{:.0f}'.format(depth_km[0]), '-', 
                     '{:.0f}'.format(depth_km[1]), ' km)'
                     )),
            rotation=-90., ha='center', va='bottom')
    ax.set_xlim(0, 4)

# diffusivity labels
labels = []
for idx in range(5, -1, -1):
    labels.append('{:.1f}'.format(logD_list[idx]))
ax1mm.text(0.2, 97., labels[0], rotation=30)
ax1mm.text(0.85, 96, labels[1], rotation=10)
ax1mm.text(3.3, 95, labels[2], rotation=0)
ax1mm.text(3.3, 79, labels[3], rotation=0)
ax1mm.text(3.3, 55, labels[4], rotation=0)
ax1mm.text(3.3, 34, labels[5], rotation=0)

labels = []
for idx in range(5, -1, -1):
    labels.append('{:.1f}'.format(logD_list[idx]))
ax3mm.text(3.2, 95., labels[0], rotation=0)
ax3mm.text(3.2, 80., labels[1], rotation=0)
ax3mm.text(3.2, 54., labels[2], rotation=0)
ax3mm.text(3.2, 32., labels[3], rotation=0)
ax3mm.text(3.2, 16., labels[4], rotation=0)
ax3mm.text(3.2, 7, labels[5], rotation=0)

ax1mm.text(1.5, 68, '~clinopyroxene', rotation=30)
ax3mm.text(3.0, 25, '~cpx', rotation=0)
ax1mm.text(0.1, 90, '~olivine', rotation=0)
ax3mm.text(1.2, 90, '~olivine', rotation=30)

plt.savefig('Fig13_xenolithEquilib.eps', format='eps', dpi=1000)
fig.savefig('Fig13_xenolithEquilib.tif', format='tif', dpi=300)

plt.show(fig)
print 'Finished'