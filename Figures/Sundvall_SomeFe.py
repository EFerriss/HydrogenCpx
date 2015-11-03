# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 08:54:58 2015

@author: Ferriss

Peak-specific hydrogen diffusivity in diopside with some Fe
Data from Sundvall et al. 2009 via PULI
# manual_peakpos = [3443 3355 3305 3250];
"""

import pynams.pynams as nams
import matplotlib.pyplot as plt
import pynams.diffusion as diffusion
import numpy as np

plt.style.use('paper')

#%%
Sundvall = nams.TimeSeries()
Sundvall.sample = 'diopside with some Fe'
Sundvall.thick_microns = 1000.
Sundvall.times_hours = [0., 4., 24., 72., 192., 528., 1428]
Sundvall.fname_list = ['Sundvall09_SomeFe_219_0', 'Sundvall09_SomeFe_219_4',
                       'Sundvall09_SomeFe_219_24',  'Sundvall09_SomeFe_219_72',
                       'Sundvall09_SomeFe_219_192', 'Sundvall09_SomeFe_219_528',
                       'Sundvall09_SomeFe_219_1428']
Sundvall.make_spectra_list()
Sundvall.initial_profile = Sundvall
Sundvall.change_baseline(highwn=3500., lowwn=3000.)

Sundvall.get_baselines()
Sundvall.get_peakfit()

#Sundvall.plot_spectra(True, False, False, True)
#Sundvall.plot_peakfits(top=20)
# Make baselines
#Sundvall.make_baselines(1)
#Sundvall.plot_spectra(True, False, False, True)
#Sundvall.save_baselines()
# manual_peakpos = [3443 3355 3305 3250];


#%% Change in area over time and diffusivities
fig = plt.figure()
fig.set_size_inches(3, 3)
ax = fig.add_subplot(111)
#ax.set_title('Sundvall et al. 2009 Fe-poor diopside\n#219 dehydrated at 800 C')
ax.set_xlabel('Time (hours)', fontsize=12)
ax.set_ylabel('Concentration/\nMaximum Concentration', fontsize=12)
fig.autofmt_xdate()

ax.text(100, 0.9, 'B.', 
        ha='left', fontsize=11, va='center')

ax.text(1550, 0.5, '$Synthetic$ $diopside$\n $800$ $\degree C$', 
        ha='right', fontsize=11, va='center')
ax.set_xlim(0, 1600)

## diffusivity curves
Ds = [-16.1, -15.9, -15.6, -14.7]
for D in Ds:
     t, cc = diffusion.diffusionThinSlab(log10D_m2s=D, 
            thickness_microns=55., max_time_hours=ax.get_xlim()[1])
     ax.plot(t, cc, '-k', linewidth=1)

ax.text(900, 0.35, Ds[0], rotation=-20)
ax.text(800, 0.25, Ds[1], rotation=-20)
ax.text(600, 0.15, Ds[2], rotation=-30)
ax.text(100, 0.1, Ds[3], rotation=-45)

# bulk H
Sundvall.make_area_list()
C = np.array(Sundvall.areas_list)
C0 = Sundvall.areas_list[0]
y = C / C0
ax.plot(Sundvall.times_hours, y, 'o', color='k', alpha=0.5,
        clip_on=False, markersize=8, label='bulk H')

# peak at 3450 cm-1
C = Sundvall.peak_areas[4, :]
C_0 = C[1]
y = C / C_0
ax.plot(Sundvall.times_hours, y, '3', color='b', clip_on=False, markersize=10, 
        label='3443 cm$^{-1}$', alpha=1, mew=2)

# peak at 3350 cm-1
C = Sundvall.peak_areas[5, :]
C_0 = C[0]
y = C / C_0
ax.plot(Sundvall.times_hours, y, 'd', color='violet', 
        clip_on=False, markersize=8, 
        label='3355 cm$^{-1}$', alpha=0.7)

plt.tight_layout()
ax.legend(loc=1, fontsize=9)

fig.savefig('../../../CpxPaper/figures/SundvallD.png', dpi=200)
#
#
#%% Spectra comparison figure
##fig = plt.figure()
##fig.set_size_inches(6.5, 4)
##
##ncol = 1
##nrow = 1
##gs = gridspec.GridSpec(nrow, ncol)
##
##axes = []
##for row in range(nrow):
##    for col in range(ncol):
##        axes.append(plt.subplot(gs[row, col]))
##
##xtickgrid=250
##ytickgrid=0.5
##xmajorLocator = MultipleLocator(xtickgrid)
##ymajorLocator = MultipleLocator(ytickgrid)
##
##axes[0].set_title('Synthetic Fe-poor diopside heated at 800 C\n Sundvall et al. 2009 #219')
##axes[0].set_xlabel('Wavenumber (cm$^{-1}$)', fontsize=14)
##axes[0].set_ylabel('Absorbance (cm$^{-1}$)', fontsize=14)
##
###axes[0].xaxis.set_major_locator(xmajorLocator)
###axes[0].yaxis.set_major_locator(ymajorLocator)
##axes[0].set_xlim(3500, 3100)
##axes[0].set_ylim(0, 12.)
##
##x  = Sundvall.spectra_list[0].base_wn
##y0 = Sundvall.spectra_list[0].abs_nobase_cm
##y1 = Sundvall.spectra_list[1].abs_nobase_cm
##y2 = Sundvall.spectra_list[2].abs_nobase_cm
##y3 = Sundvall.spectra_list[3].abs_nobase_cm
##y4 = Sundvall.spectra_list[4].abs_nobase_cm
##y5 = Sundvall.spectra_list[5].abs_nobase_cm
##y6 = Sundvall.spectra_list[6].abs_nobase_cm
##
##axes[0].plot(x, y0, label='Initial', color='b', linewidth=4)
##axes[0].plot(x, y1, label=' '.join((str(Sundvall.times_hours[1]), 'hr')), 
##                    color='m', linewidth=3.5, linestyle='--')
##axes[0].plot(x, y2, label=' '.join((str(Sundvall.times_hours[2]), 'hr')), 
##                    color='g', linewidth=3,)
##axes[0].plot(x, y3, label=' '.join((str(Sundvall.times_hours[3]), 'hr')), 
##                    color='k', linewidth=2.5, linestyle='--')
##axes[0].plot(x, y4, label=' '.join((str(Sundvall.times_hours[4]), 'hr')), 
##                    color='r', linewidth=2)
##axes[0].plot(x, y5, label=' '.join((str(Sundvall.times_hours[5]), 'hr')), 
##                    color='cadetblue', linewidth=1.5, linestyle='--')
##axes[0].plot(x, y6, label=' '.join((str(Sundvall.times_hours[6]), 'hr')), 
##                    color='saddlebrown', linewidth=1, linestyle='-')
##
##axes[0].fill_between(x, y0, y1, where=y0>=y1, facecolor='b', 
##        interpolate=True, alpha=0.2)
##axes[0].fill_between(x, y4, y5, where=y4>=y5, facecolor='r', 
##        interpolate=True, alpha=0.2)
##axes[0].fill_between(x, 0, y6, facecolor='y', alpha=0.2) 
##
##
##axes[0].legend(loc=0, fontsize=12)
##plt.tight_layout()
##fig.savefig('../../../CpxPaper/figures/RimsSundvall.png', dpi=100)
