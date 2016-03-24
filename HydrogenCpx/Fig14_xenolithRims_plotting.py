# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 15:14:02 2015

@author: Ferriss

Figure roughly following Padron-Navarta et al. 2014 Figure 8
Comparing olivine and cpx rim formation over time for 1mm and 3mm side cubes
First hour only 

Data generated and saved to xenolithRims*.txt in json format 
in xenolithRims_data*.py. File names are a little screwy, sorry.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import json

### set up and get data
### data order: 1mm cube 90% loss; 1mm cube 80% loss; 
###             3mm cube 90% loss; 3mm cube 80% loss

resolution_microns = 100. # resolution limit for most FTIR studies

length_microns = np.array([1000., 3000.])
center_microns = length_microns / 2.
volume_mm3 = (length_microns**3) / 1E9

# How much does concentration drop to qualify as a rim
water_fractions = [0.90, 0.80, 0.90, 0.8]

savefolder = 'C://Users//Ferriss//Documents//Code//Python//HydrogenCpx//HydrogenCpx//'

workfiles = []
time_minutes = []
time_hours = []
logD_list = []
distances = []

workfiles.append(''.join((savefolder, 'xenolithRims-1mm-80.txt')))
workfiles.append(''.join((savefolder, 'xenolithRims-1mm-90.txt')))
workfiles.append(''.join((savefolder, 'xenolithRims-3mm-80-1hr.txt')))
workfiles.append(''.join((savefolder, 'xenolithRims-3mm-90-1hr.txt')))

# -10.5 curves look weird for 1 mm - main_idx [0, 1], diffusivity idx = 4
workfiles.append(''.join((savefolder, 'xenolithRims-1mm-80-1hr-log10pt5.txt')))
workfiles.append(''.join((savefolder, 'xenolithRims-1mm-90-1hr-log10pt5.txt')))

# get data out of the workfiles and dump it all into data array
for idx, workfile in enumerate(workfiles):
    with open(workfile, 'r') as rimfile:
        a = rimfile.read()
        data = json.loads(a)

        time_minutes.append([])
        time_hours.append([])
        for each_distance in data[1:-1]:
            time_minutes[idx].append(data[0])
            time_hours[idx].append(np.array(data[0])/60.)

        distances.append(data[1:-1])
        logD_list.append(data[-1])

# swap out new data [4, 5][0] for old [0, 1][4] log10 -10.5 m2/s in 1 mm sample
time_minutes[0][4] = time_minutes[4][0]
time_minutes[1][4] = time_minutes[5][0]
time_hours[0][4] = time_hours[4][0]
time_hours[1][4] = time_hours[5][0]
distances[0][4] = distances[4][0]
distances[1][4] = distances[5][0]


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

ylab = ''.join(('Distance from center of cube face where\n',
                'H concentration <80-90% of initial ($\mu$m)'))
ax1mm.set_ylabel(ylab)
ax_big.set_xlabel('Time (minutes)')
ax3mm.set_yticklabels('')
fig.autofmt_xdate()

for idx, ax in enumerate([ax1mm, ax3mm]):
    tit = ''.join(('{:.0f}'.format(volume_mm3[idx]), ' mm$^3$cube'))
    ax.set_title(tit)
    ax.set_ylim(500., 0)
    ax.set_xlim(0, 60.)

#colors = ['black', 'blue', 'green', 'darkgoldenrod', 'red', 'teal', ]
#fillcolors = ['grey', 'blue', 'green', 'darkgoldenrod', 'red', 'teal', 'yellow']
colors = ['k']*6
fillcolors = ['w', 'w', 'w', 'w', 'w', 'cornflowerblue']
hatches = ['/', '/', '/', '/',  '/', None,  None]
green_alpha = 0.45

# 1mm cube
for main_idx in [0, 1]: 
    for idx, rimloss in enumerate(distances[main_idx]):
        ydata = rimloss
            
        if ydata[-1] == 0.0:
            ydata[-1] = ydata[-2]

        ax1mm.plot(time_minutes[main_idx][idx], ydata, mew=1, linewidth=1, 
                 color=colors[idx],
                 label=''.join(('logD=', 
                                '{:.1f}'.format(logD_list[main_idx][idx]))))

    for idx in xrange(0, 6):
        ax1mm.fill_between(time_minutes[main_idx][idx], distances[0][idx], 
                           distances[1][idx],facecolor=fillcolors[idx], 
                           alpha=1., interpolate=True,
                           hatch=hatches[idx])
    
ax1mm.fill_between(time_minutes[0][0], distances[0][1], distances[1][2],
                   facecolor='green', 
                   alpha=green_alpha, interpolate=True, )

# 3mm cube
for main_idx in [2, 3]: 
    for idx, rimloss in enumerate(distances[main_idx]):
        ydata = rimloss
        if ydata[-1] == 0.0:
            ydata[-1] = ydata[-2]

        ax3mm.plot(time_minutes[main_idx][idx], ydata, mew=1, linewidth=1, 
                 color=colors[idx],
                 label=''.join(('logD=', 
                                '{:.1f}'.format(logD_list[main_idx][idx]))))

for idx in xrange(0, 6):
    ax3mm.fill_between(time_minutes[main_idx][idx], distances[2][idx], 
                       distances[3][idx],facecolor=fillcolors[idx], 
                       alpha=1., interpolate=True,
                       hatch=hatches[idx])

ax3mm.fill_between(time_minutes[2][0], distances[2][1], distances[3][2],
                   facecolor='green', 
                   alpha=green_alpha, interpolate=True, )

for idx, ax in enumerate([ax1mm, ax3mm]):
    ax.axhspan(0, resolution_microns, facecolor='red', alpha=0.2)
    bottom_box_top = center_microns[idx] - 0.06*center_microns[0]

xmin = 0.
ymin = 478.
xmax = 60.
ymax = 500.
rect = patches.Rectangle((xmin,ymin), xmax-xmin, ymax-ymin, facecolor='w', 
                         alpha=1., zorder=10) #transform=fig.transFigure)
ax1mm.add_patch(rect)
ax1mm.text(30., 500, 'diffusion reaches the center', ha='center', va='bottom', zorder=11)

### time ranges for comparison
ascent_rate_m_per_s = [2., 1.]
depth_km = np.array([2., 3.])
ytloc = [[350., 270.],[320., 420.]]

depth_m = depth_km * 1E3

for main_idx, ax in enumerate([ax1mm, ax3mm]):
    ax.text(59., 30, 'Diffusive loss not\nobserved', ha='right', va='center')
    for idx, rate in enumerate(ascent_rate_m_per_s):
        ascent_time_s = depth_m / rate
        ascent_time_m = ascent_time_s / 60.
    
#        ax.axvspan(ascent_time_m[0], ascent_time_m[1], facecolor='grey', alpha=0.4)
        rect = patches.Rectangle((ascent_time_m[0], resolution_microns), 
                                 ascent_time_m[1] - ascent_time_m[0], 
                                 bottom_box_top - resolution_microns, 
                                 color='grey', alpha=0.4)
        ax.add_patch(rect)
        ax.text(np.mean(ascent_time_m), ytloc[main_idx][idx],
            ''.join(('{:.1f}'.format(rate), ' m/s\n(',
                     '{:.0f}'.format(depth_km[0]), '-', 
                     '{:.0f}'.format(depth_km[1]), ' km)')),
            rotation=-90., ha='center', va='center')


# diffusivity labels
labels = []
for idx in range(5, -1, -1):
    labels.append('{:.1f}'.format(logD_list[0][idx]))
ax1mm.text(4, 420., labels[0], rotation=-90)
ax1mm.text(15, 460, labels[1], rotation=0.)
ax1mm.text(47, 460, labels[2], rotation=0)
ax1mm.text(50, 270, labels[3], rotation=0)
ax1mm.text(50., 160, labels[4], rotation=0)
ax1mm.text(50, 95, labels[5], rotation=0)

ax3mm.text(12, 490., labels[0], rotation=0.)
ax3mm.text(38, 490, labels[1], rotation=0.)
ax3mm.text(50, 455, labels[2], rotation=0)
ax3mm.text(50, 270, labels[3], rotation=0)
ax3mm.text(50., 160, labels[4], rotation=0)
ax3mm.text(50, 96, labels[5], rotation=0)

ax1mm.text(7, 70, '~clinopyroxene', rotation=-27)
ax3mm.text(7, 70, '~clinopyroxene', rotation=-27)
ax1mm.text(3, 300, '~olivine', rotation=-80)
ax3mm.text(4, 300, '~olivine', rotation=-80)

plt.tight_layout()

plt.savefig('Fig14_xenolithRims.eps', format='eps', dpi=1000)
fig.savefig('Fig14_xenolithRims.tif', format='tif', dpi=300)

print 'Finished'
