# -*- coding: utf-8 -*-
"""
@author: Ferriss

diffusivities as a function of Fe content
Following Woods et al 2000 Figure 6

For Alex Lloyd

"""
Celsius = 1030.

bottomD = -14.
topD = -9.

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
import numpy as np
import my_spectra
import pynams.diffusivity_library as dlib
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost

#%% Bulk H data
Ingrin = dlib.H_CrDiopside_Ingrin95
Woods = dlib.H_diopside_Woods00
Xia = dlib.H_cpxBasanite_Xia00
NoFe = dlib.H_diopside_noFe
Sundvall = dlib.H_diopside_Sundvall
Kunlun = my_spectra.K_bulk
PMR = my_spectra.D_PMR

direction = 'z'

#%%################# DATA ################################
Names = [#'Kunlun diopside',
         'Jaipur diopside || a, c',
#         'augite PMR-53',
         'Russian diopside',
         unicode('N\374shan cpx || a', 'latin-1'),
         '$\leftarrow$ Fe-free diopside',
         'synthetic diopside',
         'Fuego\nphenocryst'
         ]

Fe = np.array([#0.012 + 0.013,
               0.036 + 0.039,
#               0.049 + 0.167,
               0.036, 
               0.212,
               0.,
               0.02,
               0.295])

FeOwt = np.array([#0.83,
                  2.48,
#                  7.12,
                  np.nan,
                  np.nan,
                  0.,
                  0.7,
                  8.99])

Mg = np.array([#0.015 + 0.944,
               0.012 + 0.918,
#               0.267 + 0.720,
               0.030 + 0.940,
               0.190 + 0.613,
               1., #CHECK, but doesn't matter for Mg number
               0.069 + 0.981,
               1. # didn't do, just placeholder
               ])

MgNumber = 100. * Mg / (Mg + Fe)

bulk = np.array([#Kunlun.whatIsD(Celsius, orient=direction),
        Woods.whatIsD(Celsius, orient=direction), 
#        Woods.whatIsD(Celsius, orient=direction) - 0.1, # PMR: -11.02 at 800C
        Ingrin.whatIsD(Celsius, orient=direction),
        Xia.whatIsD(Celsius, orient=direction),
        NoFe.whatIsD(Celsius, orient='x'),
        Sundvall.whatIsD(Celsius, orient='y'),
        Woods.whatIsD(Celsius, orient=direction), 
        ])

JaipurSlow = Woods.whatIsD(Celsius, orient='y')

#%%############## SETUP PLOT ###########################
x = np.log10(Fe)
fig = plt.figure(figsize=(6, 5))
gs = gridspec.GridSpec(1,1)
ax = SubplotHost(fig, 1,1,1)
fig.add_subplot(ax)
ax.set_ylim(bottomD, topD)

ax_Fe = ax.twin()
Fe_labels = list(np.arange(0.0, 0.1, 0.02)) + list(np.arange(0.1, 0.9, 0.1))
parasite_tick_locations = np.log10(Fe_labels)
ax_Fe.set_xticks(parasite_tick_locations)
ax_Fe.set_xticklabels(Fe_labels)
ax_Fe.axis["top"].set_label("Fe (a.p.f.u.)")
ax_Fe.axis["top"].label.set_visible(True)
ax_Fe.axis["right"].major_ticklabels.set_visible(False)

ax.set_ylabel(''.join(('log$_{10}$ diffusivity$_{H}$ $(m^{-2}/s)$ at ',
                       '{:.0f}'.format(Celsius), '$\degree$C')))
ax.set_xlabel('log$_{10}$ Fe (a.p.f.u.)')

label = []

### plot Jaipur slow direction
ax.plot(x[0], JaipurSlow, 'rx') #**style[1])


x[3] = -1.8 # To plot Fe free

for idx in range(5):
    # Main plotting
    ax.plot(x[idx], bulk[idx], 'rx', clip_on=False) #**style[idx])

xtext = np.array([-1.1, -1.3, -0.6, -1.8, -1.4, -0.48, -1.2])
ytext = np.array([-10.2, -11.7, -11.4, -12.35, -13.8, -10., -11.3])
hatext = ['center', 'right', 'center', 'left', 'center', 'left', 'left', 'left']

for idx in range(6):
    # create individual label names
    if np.isnan(FeOwt[idx]):
        label.append(''.join((Names[idx], '\n[Mg# ', 
                         '{:.1f}'.format(MgNumber[idx]), ']')))
    elif 'Fuego' in Names[idx]:
        label.append(''.join((Names[idx], '\n[Mg# ', '{:.1f}'.format(MgNumber[idx]),
                     ',\n', '{:.2f}'.format(FeOwt[idx]), '% FeO]')))
    else:
        label.append(''.join((Names[idx], '\n[Mg# ', '{:.1f}'.format(MgNumber[idx]),
                     ', ', '{:.2f}'.format(FeOwt[idx]), '% FeO]')))

label.append('Jaipur || b')

for idx, lab in enumerate(label):
    ax.text(xtext[idx], ytext[idx], lab, ha=hatext[idx])   

### Fuego
FuegoLoc = (x[-1], -9.5)
FuegoRange = 1.
ax.add_artist(Ellipse(FuegoLoc, 0.06, FuegoRange, facecolor='none', 
                      edgecolor='r', hatch='x'))

        
ax.set_xlim(-1.8, -0.2) 
ax.legend(loc=4, ncol=1, fontsize=10, fancybox=True)


##fig.autofmt_xdate()
plt.tight_layout()
#plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig10.eps', 
#            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\LloydCpxCoAuthor\\Fe', dpi=300)

plt.show(fig)
print 'Finished'
