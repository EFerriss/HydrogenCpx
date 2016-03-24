# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 23:29:20 2015

@author: Ferriss

diffusivities as a function of TETRAHEDRAL Al content
Focus is only on samples with Fe content > 0.05 apfu
See supplementary material for microprobe data and normalization.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
import numpy as np
from HydrogenCpx import cpx_spectra
import pynams.diffusivity_library as dlib
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost

# Bulk H data for big 4: Jaipur, Nushan, PMR, Fuego
Jaipur = dlib.H_diopside_Woods00
Nushan = dlib.H_cpxBasanite_Xia00
PMR = cpx_spectra.D_PMR

#%%
Celsius = 800.
direction = 'z'

Names = ['Jaipur diopside\n(|| a, c*)',
         unicode('N\374shan cpx', 'latin-1'),
         'augite PMR-53',
         'Fuego\nphenocryst',
         'Jaipur diopside\n(|| b)',
         ]

Al = np.array([
               0.016,
               0.204,
               0.026,
               0.133,
               0.016,])

style = [Jaipur.basestyle,
         Nushan.basestyle,
         PMR.basestyle,
         {},
         Jaipur.basestyle,
         ]

for st in style:
    st['alpha'] = 0.3
    st['markersize'] = 18
    st['mew'] = 1
    st['fillstyle'] = 'full'

style[0]['color'] = 'orange'
style[0]['markersize'] = 16
style[1]['marker'] = '>'

bulk = np.array([Jaipur.whatIsD(Celsius, orient=direction), 
                 Nushan.whatIsD(Celsius, orient=direction),
                 -11.02,
                 -11.1, 
                 Jaipur.whatIsD(Celsius, orient='y')
                 ])

x = np.log10(Al)

#%% Plotting

# setup plot, twin axes, labels
fig = plt.figure(figsize=(3, 4))
gs = gridspec.GridSpec(1,1)
ax = SubplotHost(fig, 1,1,1)
fig.add_subplot(ax)
ax.set_ylim(-12.5, -10.5)
ax_Al = ax.twin()

#Al_labels = list(np.arange(0.0, 0.1, 0.02)) + list(np.arange(0.1, 0.9, 0.1))
Al_labels = [0.01, 0.02, 0.05, 0.1, 0.2]
parasite_tick_locations = np.log10(Al_labels)
ax_Al.set_xticks(parasite_tick_locations)
ax_Al.set_xticklabels(Al_labels)

ax_Al.axis["top"].set_label("IV-Al (a.p.f.u.)")
ax_Al.axis["top"].label.set_visible(True)
ax_Al.axis["right"].major_ticklabels.set_visible(False)
ax.set_ylabel('log$_{10}$ diffusivity$_{H}$ $(m^{-2}/s)$ at 800 $\degree$C')
ax.set_xlabel('log$_{10}$ IV-Al (a.p.f.u.)')
ax.set_xlim(-2.0, -0.4) 

Names = ['Jaipur diopside\n(|| a, c*)',
         unicode('N\374shan cpx', 'latin-1'),
         'augite PMR-53',
         'Fuego\nphenocryst',
         'Jaipur diopside\n(|| b)',
         ]

Al = np.array([
               0.016,
               0.204,
               0.026,
               0.133,
               0.016,])

bulk = np.array([Jaipur.whatIsD(Celsius, orient=direction), 
                 Nushan.whatIsD(Celsius, orient=direction),
                 -11.02,
                 -11.1, 
                 Jaipur.whatIsD(Celsius, orient='y')
                 ])

# Jaipur diopside fast direction
ax.plot(np.log10(0.016), Jaipur.whatIsD(Celsius, orient='x'), **Jaipur.basestyle)
ax.text(-1.55, -10.8, 'Jaipur diopside\nfast direction', ha='center')

# augite PMR-53
ax.plot(np.log10(0.026), -11.02, **PMR.basestyle)
ax.text(-1.75, -11.3, 'augite\nPMR-53',)

# Nushan cpx
ax.plot(np.log10(0.204), Nushan.whatIsD(Celsius, orient='z'), **Nushan.basestyle)
ax.text(-0.7, -12.35, unicode('N\374shan\ncpx', 'latin-1'), ha='center')

# Fuego cpx
FuegoLoc = (np.log10(0.133), (Jaipur.whatIsD(Celsius, orient='x')+Jaipur.whatIsD(Celsius, orient='y'))/2.)
ax.add_artist(Ellipse(FuegoLoc, 0.06, 1., facecolor='none', edgecolor='k'))
ax.text(-0.95, -11.5, 'Fuego\nphenocryst', ha='right')

# Xenoliths
xenoMin = np.log10(0.11 - 0.035)
xenoMax = np.log10(0.11 + 0.035)
xenoDrange = 1.
xenoDmin = Nushan.whatIsD(Celsius, orient='z')-0.05
fbox = mpatches.FancyBboxPatch([xenoMin, xenoDmin], xenoMax-xenoMin, xenoDrange,
                               boxstyle=mpatches.BoxStyle("Round", pad=0.02),
                               fill=True, alpha=0.25, color='r')
ax.add_patch(fbox)
ax.text(xenoMin, xenoDmin, '~mantle\ncpx', ha='left', va='bottom',
        rotation=90.)


#ax.grid(True)
#plt.tight_layout()
plt.subplots_adjust(right=0.95, left=0.26, bottom=0.0, top=0.9)
fig.autofmt_xdate()
            
plt.savefig('Fig12_Al.eps', format='eps', dpi=1000)
fig.savefig('Fig12_Al.tif', format='tif', dpi=300)

plt.show(fig)
print 'Finished'
