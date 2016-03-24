# -*- coding: utf-8 -*-
"""
@author: Ferriss

"""
import cpx_spectra
import numpy as np
import matplotlib.pyplot as plt
import pynams.diffusivity_library as dlib
import pynams.styles as styles

#%% Import the data 
reload(cpx_spectra)

# diffusivities
Kunlun_91hr = cpx_spectra.D_K91
Kunlun_154hr = cpx_spectra.D_K154
Jaipur = cpx_spectra.D_J

# errors
Kunlun_91hr_errors = cpx_spectra.e_K91
Kunlun_154hr_errors = cpx_spectra.e_K154
Jaipur_errors = cpx_spectra.e_J

# move bulk H from the end to the beginning for all
for x in [Kunlun_154hr, Kunlun_154hr_errors, Kunlun_91hr, Kunlun_91hr_errors,
          Jaipur, Jaipur_errors]:
    x.insert(0, x.pop())

#%% The plot
fig = plt.figure(figsize=(6., 3.3))
ax = fig.add_subplot(111)

ax.set_xlim(-0.5, 11.5)
ax.set_ylim(-15, -8)

# temperature of interest
temp_celcius = 904.

# Set where they plot on the x-axis
x = np.arange(0.5, 11.5, 2)

# Kunlun data
ax.errorbar(x, Kunlun_91hr, yerr= Kunlun_91hr_errors, 
            ecolor='black', fillstyle='none', marker='.', color='k', mew=1,
            linestyle='none', label='Kunlun diopside\n(this study)')

# Jaipur data
ax.errorbar(x, Jaipur, yerr=Jaipur_errors,
            ecolor='r', fillstyle='none', marker='s', color='r', mew=1,
            linestyle='none', markersize=7, label='Jaipur diopside\n(this study)')

Woods = dlib.H_diopside_Woods00.whatIsD(904, orient='x')
ax.plot(0.5, Woods, markersize=12, alpha=0.4, marker='s',
        linestyle='none',
        label='Jaipur diopside calculated\nfrom Arrhenius relation\nof Woods et al. 2000')

ax.legend(loc=8, fancybox=True, ncol=3, bbox_to_anchor=(0.45, -0.185),
          columnspacing=0.2)

# More Kunlun data after legend formed
ax.errorbar(x, Kunlun_154hr, yerr= Kunlun_154hr_errors, 
            ecolor='black', fillstyle='none', marker='.', color='k', mew=1,
            linestyle='none', label='Kunlun diopside\n(this study)')

### Labels and commentary ###
dividers = np.arange(1.5, 51.5, 2)
for k in range(5):
    ax.plot([dividers[k], dividers[k]], ax.get_ylim(), 'k')

# no x labels at the bottom
plt.xticks([100])

ax.set_ylabel("log$_{10}$D$_c$ (m$^{2}$/s) at 904 $\degree$C")

# peak labels at the top
peaklabels = ['bulk H', '3645\ncm$^{-1}$', '3617\ncm$^{-1}$',
              '3540\ncm$^{-1}$',
#              '~3450\ncm$^{-1}$', 
              '3443\ncm$^{-1}$', 
              '3350\ncm$^{-1}$']
x_peaklabels = np.arange(0.5, 100.5, 2)
for k in range(len(peaklabels)):
    ax.text(x_peaklabels[k], -8.7, peaklabels[k], ha='center',
            va='center',
            backgroundcolor='none', fontsize=12)

Jstart = -10.4
Kstart = -14.2
Jdiff = 1.
Kdiff = 2.
plt.axhspan(Jstart, Jstart + Jdiff, facecolor='palegreen')
plt.axhspan(Kstart, Kstart + Kdiff, facecolor='khaki')

tgap = 0.3
ax.text(6.5, Kstart + tgap, '${Kunlun}$', fontsize=16, ha='center', va='center')
ax.text(6.5, Jstart + tgap, '${Jaipur}$', fontsize=16, ha='center', va='center')
      
plt.savefig('Fig8_PeakComparison.eps', format='eps', dpi=1000)
fig.savefig('Fig8_PeakComparison.tif', format='tif', dpi=300)
