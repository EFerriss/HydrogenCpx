# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 23:29:20 2015

@author: Ferriss

diffusivities as a function of Fe content
Following Woods et al 2000 Figure 6
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
import numpy as np
import cpx_spectra
import pynams.diffusivity_library as dlib
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost

reload(cpx_spectra)

# Bulk H data
Ingrin = dlib.H_CrDiopside_Ingrin95
Woods = dlib.H_diopside_Woods00
Xia = dlib.H_cpxBasanite_Xia00
NoFe = dlib.H_diopside_noFe
Sundvall = dlib.H_diopside_Sundvall
Kunlun = cpx_spectra.K_bulk
PMR = cpx_spectra.D_PMR

# Peak-specific differences: 3645, 3617, 3540, 3443, 3350, BULK H
K3 = cpx_spectra.D_K3
K4_154 = cpx_spectra.D_K154
K4_91 = cpx_spectra.D_K91
K5 = cpx_spectra.D_K5
J = cpx_spectra.D_J

Kdifferences = np.ones((4, 5))
idx = 0
for K in [K3, K4_154, K4_91, K5]:
    Kdifferences[idx, :] = K[-1] - np.array(K[0:5])
    idx = idx + 1
pwK = np.mean(Kdifferences, axis=0)
pwJ = J[-1] - np.array(J[0:5]) 

###print '\nKunlun average difference in log10D from bulk, high to low wn'
for p in pwJ:
    print -1*p
#print pwK
#print '\nJaipur average difference in log10D from bulk, high to low wn'
#print pwJ

#%%################# DATA ################################
Celsius = 800.
direction = 'z'

Names = ['Kunlun diopside',
         'Jaipur diopside',
         'augite PMR-53',
         'Russian diopside',
         unicode('N\374shan cpx', 'latin-1'),
         'Fe-free diopside',
         'synthetic diopside',
         'Fuego\nphenocryst'
         ]

Fe = np.array([0.012 + 0.013,
               0.036 + 0.039,
               0.049 + 0.167,
               0.036, 
               0.212,
               0.,
               0.02,
               0.295])

FeOwt = np.array([0.83,
                  2.48,
                  7.12,
                  np.nan,
                  np.nan,
                  0.,
                  0.7,
                  8.99])

Mg = np.array([0.015 + 0.944,
               0.012 + 0.918,
               0.267 + 0.720,
               0.030 + 0.940,
               0.190 + 0.613,
               1., #CHECK, but doesn't matter for Mg number
               0.069 + 0.981,
               1. # didn't do, just placeholder
               ])

MgNumber = 100. * Mg / (Mg + Fe)

style = [Kunlun.basestyle, 
         Woods.basestyle,
         PMR.basestyle,
         Ingrin.basestyle, 
         Xia.basestyle,
         NoFe.basestyle,
         Sundvall.basestyle,
         {}
#         Woods.basestyle.copy()
         ]

for st in style:
    st['alpha'] = 0.3
    st['markersize'] = 18
    st['mew'] = 1
    st['fillstyle'] = 'full'

style[0]['marker'] = 'o'
style[1]['marker'] = 's'
style[3]['marker'] = '^'
style[4]['marker'] = '>'
style[5]['marker'] = '<'
style[6]['marker'] = 'v'

style[2]['color'] = 'orange'
style[2]['markersize'] = 16
#style[7]['marker'] = '8'
#style[7]['markerfacecolor'] = 'k'
#style[7]['alpha'] = 0.7
#style[7]['markersize'] = 12

bulk = np.array([Kunlun.whatIsD(Celsius, orient=direction),
        Woods.whatIsD(Celsius, orient=direction), 
        -11.02,
        Ingrin.whatIsD(Celsius, orient=direction),
        Xia.whatIsD(Celsius, orient=direction),
        NoFe.whatIsD(Celsius, orient='x'),
        Sundvall.whatIsD(Celsius, orient='y'),
        Woods.whatIsD(Celsius, orient=direction), 
        ])

Peak3645 =  np.array([Kunlun.whatIsD(Celsius, orient=direction) - pwK[0],
            Woods.whatIsD(Celsius, orient=direction) - pwJ[0],
            -11.04,
            Ingrin.whatIsD(Celsius, orient=direction) - pwK[0],
            Xia.whatIsD(Celsius, orient=direction) - pwJ[0],
            np.nan,
            np.nan,
            np.nan,
            ])

Peak3540 =  np.array([Kunlun.whatIsD(Celsius, orient=direction) - pwK[2],
            Woods.whatIsD(Celsius, orient=direction) - pwJ[2],
            -11.06,
            Ingrin.whatIsD(Celsius, orient=direction) + 0.3,
            Xia.whatIsD(Celsius, orient=direction)  - pwJ[2],
            np.nan,
            np.nan,
            np.nan,
            ])

Peak3450 =  np.array([Kunlun.whatIsD(Celsius, orient=direction)  - pwK[3],
            Woods.whatIsD(Celsius, orient=direction)  - pwJ[3],
            -11.01,
            np.nan,
            np.nan,
            np.nan,
            -16.1,
            np.nan,
            ])

Peak3330 =  np.array([Kunlun.whatIsD(Celsius, orient=direction)  - pwK[4],
            Woods.whatIsD(Celsius, orient=direction)  - pwJ[4],
            np.nan,
            np.nan,
            np.nan,
            NoFe.whatIsD(Celsius, orient='x'),    
            -15.6,
            np.nan,
            ])

JaipurSlow = Woods.whatIsD(Celsius, orient='y')
#%%############## SETUP PLOT ###########################
#plt.style.use('paper')
#
#    fig = plt.figure(figsize=figsize_inches)
#    ax = SubplotHost(fig, 1,1,1)
#    ax_celsius = ax.twin()
#    parasite_tick_locations = 1e4/(celsius_labels + 273.15)
#    ax_celsius.set_xticks(parasite_tick_locations)
#    ax_celsius.set_xticklabels(celsius_labels)
#    fig.add_subplot(ax)
#    ax.axis["bottom"].set_label("10$^4$/Temperature (K$^{-1}$)")
#    ax.axis["left"].set_label("log$_{10}$diffusivity (m$^{2}$/s)")
#    ax_celsius.axis["top"].set_label("Temperature ($\degree$C)")
#    ax_celsius.axis["top"].label.set_visible(True)
#    ax_celsius.axis["right"].major_ticklabels.set_visible(False)
#    ax.set_xlim(low, high)
#    ax.set_ylim(bottom, top)
x = np.log10(Fe)
fig = plt.figure(figsize=(6, 5))
gs = gridspec.GridSpec(1,1)
#ax = plt.subplot(gs[0, 0])
ax = SubplotHost(fig, 1,1,1)
fig.add_subplot(ax)
ax.set_ylim(-17, -10)

ax_Fe = ax.twin()
Fe_labels = list(np.arange(0.0, 0.1, 0.02)) + list(np.arange(0.1, 0.9, 0.1))
parasite_tick_locations = np.log10(Fe_labels)
ax_Fe.set_xticks(parasite_tick_locations)
ax_Fe.set_xticklabels(Fe_labels)
ax_Fe.axis["top"].set_label("Fe (a.p.f.u.)")
ax_Fe.axis["top"].label.set_visible(True)
ax_Fe.axis["right"].major_ticklabels.set_visible(False)


ax.set_ylabel('log$_{10}$ diffusivity$_{H}$ $(m^{-2}/s)$ at 800 $\degree$C')
ax.set_xlabel('log$_{10}$ Fe (a.p.f.u.)')


label = ['name'] * 10

#for idx in range(len(Names)):
for idx in [0, 1, 2, 3, 4, 6, 7]:
    ax.plot(x[idx], bulk[idx], #label=Names[idx], 
            clip_on=False, **style[idx])

# individual labels
    xyloc = (x[idx], bulk[idx])
    if np.isnan(FeOwt[idx]):
        label[idx] = ''.join((Names[idx], '\n[Mg# ', 
                         '{:.1f}'.format(MgNumber[idx]), ']'))
    elif 'Fuego' in Names[idx]:
        label[idx] = ''.join((Names[idx], '\n[Mg# ', '{:.1f}'.format(MgNumber[idx]),
                     ',\n', '{:.2f}'.format(FeOwt[idx]), '% FeO]'))
    else:
        label[idx] = ''.join((Names[idx], '\n[Mg# ', '{:.1f}'.format(MgNumber[idx]),
                     ', ', '{:.2f}'.format(FeOwt[idx]), '% FeO]'))

#        label = ''.join((Names[idx], '\nMg# ', '{:.1f}'.format(MgNumber[idx]),
#                         '\nFeO wt%: ', 
#                        '{:.2f}'.format(FeOwt[idx])))
#    else:
#        xyloc = (x[idx], -11.5)
#        label = ''.join((Names[idx], '\nMg# 74-78',
#                         '\nFeO wt %: ', 
#                        '{:.2f}'.format(FeOwt[idx])))

#    ax.annotate(label[idx], xy=xyloc, xytext=xytextloc[idx],
#                arrowprops=dict(facecolor='black', arrowstyle='->'))

ax.text(-1.55, -14.4, label[0], ha='left')
ax.text(-1.175, -11.2, label[1], ha='right')
ax.text(-0.7, -10.65, label[2], ha='center')
ax.text(-1.38, -13.2, label[3], ha='left')
ax.text(-0.7, -12.9, label[4], ha='center')
ax.text(-1.65, -16.1, label[6], ha='left')
ax.text(-0.49, -11.95, label[7], ha='left')

### Jaipur slow direction
ax.plot(x[1], JaipurSlow, **style[1])
ax.text(-1.175, JaipurSlow, 'Jaipur || b', ha='right', va='center')

### Fuego
FuegoLoc = FuegoLoc = (x[-1], (bulk[1] + JaipurSlow)/2.)
ax.add_artist(Ellipse(FuegoLoc, 0.06, 1., facecolor='none'))

#### Xenoliths
xenoMin = np.log10(0.08)
xenoMax = np.log10(0.1)
xenoDrange = 1.
xenoDmin = -12.0638 # from Al2.py
fbox = mpatches.FancyBboxPatch([xenoMin, xenoDmin], xenoMax-xenoMin, xenoDrange,
                               boxstyle=mpatches.BoxStyle("Round", pad=0.02),
                               fill=True, alpha=0.25, color='r')
#                               label='peridodite mantle cpx?')
ax.add_patch(fbox)
ax.text(-1.07, xenoDmin + 0.5*xenoDrange, 
        'peridotite\nmantle\ncpx?', ha='left', va='center')
        
#
#### Diopside
#
#
########### plot all peak-specific together #################
ax.plot(x, Peak3645, label='~3645 cm$^{-1}$', markersize=12, mew=2,
        color='r', marker='x', linestyle='none', clip_on=False)
ax.plot(x, Peak3540, label='~3540 cm$^{-1}$', markersize=7, mew=1,
        color='teal', marker='s', linestyle='none', clip_on=False, alpha=0.6)
ax.plot(x, Peak3450, '3', label='~3450 cm$^{-1}$', markersize=12, mew=2,
        color='blue', clip_on=False)
ax.plot(x, Peak3330, 'd', label='~3330 cm$^{-1}$', markersize=8, mew=1,
        color='violet', clip_on=False, alpha=0.5)

ax.set_xlim(-1.8, -0.2) 
ax.legend(loc=4, ncol=1, fontsize=10, fancybox=True)


##fig.autofmt_xdate()
plt.tight_layout()
plt.savefig('Fig10_Fe.eps', format='eps', dpi=1000)
fig.savefig('Fig10_Fe.tif', format='tif', dpi=300)

plt.show(fig)
print 'Finished'
