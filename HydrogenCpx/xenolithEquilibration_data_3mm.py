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

lengths_microns = [3000.] * 3
logD_list = [-12.5, -12., -11.5, -11., -10.5, -10.]
time_minutes = np.linspace(0.00001, 25., 75) * 60.
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

workfile = ''.join((savefolder, '-xenolithEquilibration_3mm.txt'))
with open(workfile, 'w') as diff_file:
    diff_file.write(json.dumps(list((data))))

#%% Plotting
fig = plt.figure()
fig.set_size_inches(6, 5)
ax = fig.add_subplot(111)
plt.style.use('paper')

colors = ['green', 'b', 'purple', 'red', 'black']

#ax.plot(time_hours, data[1])

#ax.plot(time_hours, data[k])
idx_D = 0
for percentRemaining in data[1:-1]:
    ax.plot(time_hours, percentRemaining, '-', mew=1, linewidth=2, )
#            color=colors[idx_D],
#            label=''.join(('logD=', '{:.1f}'.format(logD_list[idx_D]))))
    idx_D = idx_D + 1


volume_mm3 = (lengths_microns[0] * lengths_microns[1] * lengths_microns[2]) / 1E9
tit = ''.join(('{:.1f}'.format(volume_mm3), ' mm$^3$cube'))

ylab = 'Extent of "water" re-equilibration (%)'
ax.set_ylabel(ylab)
ax.set_xlabel('Time (hours)')
ax.set_title(tit)
ax.set_ylim(0., 100.)
ax.set_xlim(0, max(time_hours))

### time ranges for comparison
#ascent_rate_m_per_s = [0.1, 0.2, 0.5, 1., 2., 5., 10.]
#depth_km = np.array([2., 3.])
#ytloc = [300] * len(ascent_rate_m_per_s)
#
#depth_m = depth_km * 1E3
#idx = 0
#for rate in ascent_rate_m_per_s:
#    ascent_time_s = depth_m / rate
#    ascent_time_m = ascent_time_s / 60.


fig.savefig(''.join((savefolder, 'xenolithEquilibration_3mm.png')), dpi=200)
plt.show(fig)
print 'Finished'