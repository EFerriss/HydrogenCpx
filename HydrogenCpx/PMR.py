# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 14:37:34 2015

@author: Ferriss
"""

import pynams.pynams as nams
import pynams.styles as styles
import numpy as np
import matplotlib.pyplot as plt

#%%
diff = nams.diffusion

PMR_dehydrated = nams.Sample()
PMR_dehydrated.sample_thick_microns = np.mean([885., 867., 873., 878., 879.])

class PMR53(nams.TimeSeries):
    def __init__(self, fnames, hours):
        self.sample = PMR_dehydrated
        self.thick_microns = np.mean(self.sample.sample_thick_microns)
        self.fname_list = fnames
        self.times_hours = hours
        
PMR_unpol = PMR53(fnames=['P_0_unpol', 
                          'P_1_unpol', 
                          'P_2_unpol',
                          'P_3_unpol', 
                          'P_4_unpol', 
                          'P_5_unpol', 
                          'P_6_unpol',],
                  hours=[0., 0.25, 0.5, 0.75, 1., 2., 3.])

PMR_unpol.make_spectra_list()
PMR_unpol.initial_profile = PMR_unpol

PMR_unpol.get_baselines()
PMR_unpol.make_area_list()
PMR_unpol.get_peakfit()

#%%
plt.style.use('paper')


D = -11.02
Drange = 0.2

Ds = [-11.5, -11., -10.5]
sty = {'color' : 'k', 'linestyle' : 'None', 'marker' : 'o', 'alpha' : 0.5}

xmax = 4.
fig, ax = PMR_unpol.plot_timeseries(D_list=Ds, style=sty, max_hours=xmax)

ax.text(0.2, 1.25, 'A. $augite$\n$PMR-53$ \n $800$ $\degree C$', 
        ha='left', fontsize=11, va='center')

xloc =  3.9
ylocs = [0.5, 0.2, 0.02]

for k in range(len(Ds)):
    ax.text(xloc, ylocs[k], Ds[k], ha='right')
#ax.text(1.2, 0.35, Ds[2], rotation=-30)
#ax.text(1.4, 0.45, Ds[1], rotation=-30)
#ax.text(1.6, 0.55, Ds[0], rotation=-20)

plt.tight_layout()

# bulk H
PMR_unpol.make_area_list()
C = np.array(PMR_unpol.areas_list)
C0 = PMR_unpol.areas_list[0]
y = C / C0
ax.plot(PMR_unpol.times_hours, y, 'o', color='k', alpha=0.5,
        clip_on=False, markersize=8, label='bulk H')

# peak at 3620
C = PMR_unpol.peak_areas[0, :]
C_0 = C[0]
y = C / C_0

ax.plot(PMR_unpol.times_hours, y, 
        '+', color='orange', clip_on=False, markersize=10, 
        label='3620 cm$^{-1}$', alpha=1, mew=2)

# peak at 3550
C = PMR_unpol.peak_areas[1, :]
C_0 = C[0]
y = C / C_0
ax.plot(PMR_unpol.times_hours, y, 
        's', color='teal', clip_on=False, markersize=8, 
        label='3540 cm$^{-1}$', alpha=0.5, mew=1)

# peak at 3460 cm-1
C = PMR_unpol.peak_areas[2, :]
C_0 = C[0]
y = C / C_0
ax.plot(PMR_unpol.times_hours, y, color='g', clip_on=False, markersize=10, 
        label='3460 cm$^{-1}$', alpha=1, mew=2, marker='3', linestyle='none')

# peak at 3355 cm-1
C = PMR_unpol.peak_areas[3, :]
C_0 = C[0]
y = C / C_0
ax.plot(PMR_unpol.times_hours, y, 'd', color='violet', 
        clip_on=False, markersize=8, 
        label='3355 cm$^{-1}$', alpha=0.7)

ax.legend(loc=1, ncol=1, fontsize=9)
ax.set_xlim(0, xmax)
ax.set_ylim(0, 1.5)


fig.savefig('../../../../CpxPaper/figures/PMR_unpol.png', dpi=1200)

fig.show()
print 'Finished'

