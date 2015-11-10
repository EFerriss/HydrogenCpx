# -*- coding: utf-8 -*-
"""
@author: Ferriss

"""
import my_spectra
import numpy as np
import matplotlib.pyplot as plt
import pynams.diffusivity_library as dlib
import pynams.styles as styles
#import diffusivity_library 

#%% Import the data and style information
Jaipur_bulkH = my_spectra.D_Jaipur_bulkH
Jaipur_peak0 = my_spectra.D_Jaipur_peak0
Jaipur_peak1 = my_spectra.D_Jaipur_peak1
Jaipur_peak2 = my_spectra.D_Jaipur_peak2
Jaipur_peak3 = my_spectra.D_Jaipur_peak3
Jaipur_peak4 = my_spectra.D_Jaipur_peak4
Jaipur_peak5 = my_spectra.D_Jaipur_peak5

Kunlun_bulkH = my_spectra.K_bulk
Kunlun_peak0 = my_spectra.K_3645
Kunlun_peak1 = my_spectra.K_3617
Kunlun_peak2 = my_spectra.K_3540
Kunlun_peak4 = my_spectra.K_3443
Kunlun_peak5 = my_spectra.K_3350

#%%

# temperature of interest
temp_celcius = 904.

dlist = [
         Kunlun_bulkH, 
         Jaipur_bulkH, 
         Kunlun_peak0, 
         Jaipur_peak0, # VERY SMALL AREA
         Kunlun_peak1, 
         Jaipur_peak1, 
         Kunlun_peak2,
         Jaipur_peak2,
         Kunlun_peak4,
#         Jaipur_peak3,
         Jaipur_peak4,
         Kunlun_peak5, # SMALL AREA ON EDGE; MIXED UP WITH BACKGROUND AND PEAK 4
         Jaipur_peak5,
         ]

# Set where they plot on the x-axis
#x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11]
x = range(12)

labels = {Kunlun_bulkH : 'Kunlun',
          Kunlun_peak0 : 'Kunlun',
          Kunlun_peak1 : 'Kunlun',
          Kunlun_peak2 : 'Kunlun',
#          Kunlun_peak3 : 'Kunlun 3443 cm$^{-1}$',
          Kunlun_peak4 : 'Kunlun',
          Kunlun_peak5 : 'Kunlun',        
          Jaipur_bulkH : 'Jaipur',
          Jaipur_peak0 : 'Jaipur',
          Jaipur_peak1 : 'Jaipur',
          Jaipur_peak2 : 'Jaipur',
#          Jaipur_peak3 : 'Jaipur 3460 cm$^{-1}$',
#          Jaipur_peak3 : 'Jaipur',
          Jaipur_peak4 : 'Jaipur',
          Jaipur_peak5 : 'Jaipur',}
          
dividers = np.arange(1.5, 51.5, 2)

fig = plt.figure(figsize=(6., 3.3))
ax = fig.add_subplot(111)
#ax.yaxis.grid()

ax.set_xlim(-0.5, 11.5)
ax.set_ylim(-15, -8)

k = 0
xlabels = []
styleDx = []
x_explicit = []
D_explicit = []
error_high_explicit = []
error_low_explicit = []

for group in dlist:
    a = group.celsius_all
    idx_list = [item for item in range(len(a)) if a[item] == temp_celcius]

    Dx = []    
    Dxe = []
    for idx in idx_list:
        Dx.append(group.logDx[idx])
        Dxe.append(group.logDx_error[idx])

    style = styles.style_points
    style['markersize'] = 10
    
#    ax.errorbar(np.ones_like(Dx)*x[k], Dx, yerr=Dxe, ecolor='black', **style)

    x_explicit = x_explicit + list(np.ones_like(Dx)*x[k])
    D_explicit = D_explicit + Dx
    error_high_explicit = error_high_explicit + Dxe
    error_low_explicit = error_low_explicit + Dxe
    
    xlabels.append(labels[group])

    k = k + 1

###### Explicitly set values, including errors ############
### Default x_explicit and D_explicit from loops and files, but can be set here.
#error_high_explicit[3] = 2. # K 3645 91hr
# peak_D = [-13.6, -12.6, -12.5, 1., -13.1, -13.] #K4 91
# peak_D = [-13.6, -12.6, -12.5, 1., -13.1, -13.] # K4 154 
# peak_D = [-10.3, -10., -9.5, 0., -9.9, -9.75] # J

### See CpxErrorbars.xlsx
x = [0.5, 2.5, 4.5, 6.5, 8.5, 10.5]

D_K154 = my_spectra.D_K154[:]
D_K91 = my_spectra.D_K91[:]
D_J = my_spectra.D_J[:]
e_K154 = my_spectra.e_K154[:] 
e_K91 = my_spectra.e_K91[:] 
e_J = my_spectra.e_J[:] 

# Move bulk from end to beginning of list
for D in [D_K154, e_K154, D_K91, e_K91, D_J, e_J]:
    D.insert(0, D.pop(-1))

### Bulk H from linear mixing of peaks: MixingPeaks.py ###
x_mixJ = [0.5]
x_mixK = [0.5, 0.5]
mixJ = [-9.83]
mixJ_e = [0.12]
mixK = [-13.02, -13.10]
mixK_e = [0.17, 0.17]

x_K = x + x #+ x_mixK
D_K = D_K154 + D_K91 #+ mixK
e_K = e_K154 + e_K91 #+ mixK_e

ax.errorbar(x_K, D_K, yerr= e_K, 
            ecolor='black', fillstyle='none', marker='.', color='k', mew=1,
            linestyle='none', label='Kunlun diopside\n(this study)')

ax.errorbar(x, D_J, yerr=e_J,
            ecolor='r', fillstyle='none', marker='s', color='r', mew=1,
            linestyle='none', markersize=7, label='Jaipur diopside\n(this study)')

Woods = dlib.H_diopside_Woods00.whatIsD(904, orient='x')
ax.plot(0.5, Woods, markersize=12, alpha=0.4, marker='s',
        linestyle='none',
        label='Jaipur diopside calculated\nfrom Arrhenius relation\nof Woods et al. 2000')

ax.legend(loc=8, fancybox=True, ncol=3, bbox_to_anchor=(0.45, -0.185),
          columnspacing=0.2)
###############################################################

for k in range(5):
    ax.plot([dividers[k], dividers[k]], ax.get_ylim(), 'k')

# label Kunlun or Jaipur
xlabels.append('Kunlun')
plt.xticks(x, xlabels, rotation=30, ha='right')

# or no x labels at the bottom
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
            
#%%            
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig8.eps', 
            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig8.png', dpi=300)
