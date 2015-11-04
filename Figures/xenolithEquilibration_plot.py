# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 15:14:02 2015

@author: Ferriss

Figure roughly following Padron-Navarta et al. 2014 Figure 8
Comparing olivine and cpx rim formation over time

Data generated and saved to xenolithEquilibrium.txt in json format 
in xenolithEquilibrium_data.py

See also xenolithRims_*.py
"""
import matplotlib.pyplot as plt
import json
import numpy as np

lengths_microns = [1000.] * 3

savefolder = '../../../CpxPaper/figures/'
workfile = ''.join((savefolder, '-xenolithEquilibration.txt'))

with open(workfile, 'r') as rimfile:
    a = rimfile.read()
data = json.loads(a)

time_minutes = data[0]
time_hours = np.array(time_minutes) / 60.
logD_list = data[-1]

#%% Plotting
fig = plt.figure()
fig.set_size_inches(6, 5)
ax = fig.add_subplot(111)
plt.style.use('paper')

#colors = ['green', 'b', 'purple', 'red', 'black']
#colors = ['blue']*5 + ['black']

#ax.plot(time_hours, data[1])


#ax.plot(time_hours, data[k])
idx_D = 0
for percentRemaining in data[1:-1]:
    ax.plot(time_hours, percentRemaining, '-', mew=1, linewidth=2, 
#            color=colors[idx_D],
#            color=colors[idx_D],
            label=''.join(('logD=', '{:.1f}'.format(logD_list[idx_D]))))
    idx_D = idx_D + 1


volume_mm3 = (lengths_microns[0] * lengths_microns[1] * lengths_microns[2]) / 1E9
tit = ''.join(('{:.1f}'.format(volume_mm3), ' mm$^3$cube'))

ylab = 'Extent of "water"\nre-equilibration with surroundings (%)'
ax.set_ylabel(ylab)
ax.set_xlabel('Time (hours)')
ax.set_title(tit)
ax.set_ylim(0., 100.)
ax.set_xlim(0, max(time_hours))

### time ranges for comparison
ascent_rate_m_per_s = [20., 10., 6., 4.]
depth_km = np.array([40., 50.])
ytloc = [2] * len(ascent_rate_m_per_s)

depth_m = depth_km * 1E3

idx = 0
for rate in ascent_rate_m_per_s:
    ascent_time_s = depth_m / rate
    ascent_time_m = ascent_time_s / 60.
    ascent_time_h = ascent_time_m / 60.

    plt.axvspan(ascent_time_h[0], ascent_time_h[1], facecolor='grey', alpha=0.2)

    ax.text(np.mean(ascent_time_h), ytloc[idx], 
        ''.join(('{:.1f}'.format(rate), ' m/s (',
                 '{:.0f}'.format(depth_km[0]), '-', 
                 '{:.0f}'.format(depth_km[1]), ' km)')),
        rotation=-90., ha='center', va='bottom')
    idx = idx + 1

ax.set_xlim(0, 4)

labels = []
for idx in range(5, -1, -1):
#    labels.append(''.join(('logD=', '{:.1f}'.format(logD_list[idx]), ' m$^2$','/s')))
    labels.append(''.join(('10$^{','{:.1f}'.format(logD_list[idx]), '}$ m$^2$','/s')))
#labels[0] = ''.join(('~ol.\n', labels[0]))
ax.text(0.1, 97., labels[0], rotation=30)
ax.text(0.8, 96, labels[1], rotation=10)
ax.text(3.2, 96, labels[2], rotation=3)
ax.text(3.2, 82, labels[3], rotation=10)
ax.text(3.2, 57, labels[4], rotation=10)
ax.text(3.2, 36.5, labels[5], rotation=10)

ax.text(0.05, 97, '~olivine', rotation=30)
ax.text(1., 88, '~clinopyroxene', rotation=25)
#ax.legend(loc=4)

fig.savefig(''.join((savefolder, 'xenolithEquilibration.png')), dpi=200)
plt.show(fig)
print 'Finished'