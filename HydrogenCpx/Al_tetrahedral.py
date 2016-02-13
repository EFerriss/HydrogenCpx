# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 23:29:20 2015

@author: Ferriss

diffusivities as a function of TETRAHEDRAL Al content
Focus only on Fe content > 0.05 apfu
See Al.py and Al2.py for more Al contents
See Fe.py for peak-specific
See spreadsheet cpx_normalization for data
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
import numpy as np
from HydrogenCpx import my_spectra
import pynams.diffusivity_library as dlib
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost

# Bulk H data for big 4: Jaipur, Nushan, PMR, Fuego
Jaipur = dlib.H_diopside_Woods00
Nushan = dlib.H_cpxBasanite_Xia00
PMR = my_spectra.D_PMR

#%%
Celsius = 800.
direction = 'z'

Names = ['Jaipur diopside\n(|| a, c*)',
         unicode('N\374shan cpx', 'latin-1'),
         'augite PMR-53\n(not oriented)',
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

#%%
bulk = np.array([Jaipur.whatIsD(Celsius, orient=direction), 
                 Nushan.whatIsD(Celsius, orient=direction),
                 -11.02,
                 -11.1, 
                 Jaipur.whatIsD(Celsius, orient='y')
                 ])

x = np.log10(Al)
#%% Plot
fig = plt.figure(figsize=(6, 5))
gs = gridspec.GridSpec(1,1)
ax = SubplotHost(fig, 1,1,1)
fig.add_subplot(ax)
ax.set_ylim(-12.5, -10.5)

ax_Al = ax.twin()
Al_labels = list(np.arange(0.0, 0.1, 0.02)) + list(np.arange(0.1, 0.9, 0.1))
parasite_tick_locations = np.log10(Al_labels)
ax_Al.set_xticks(parasite_tick_locations)
ax_Al.set_xticklabels(Al_labels)
ax_Al.axis["top"].set_label("IV-Al (a.p.f.u.)")
ax_Al.axis["top"].label.set_visible(True)
ax_Al.axis["right"].major_ticklabels.set_visible(False)

ax.set_ylabel('log$_{10}$ diffusivity$_{H}$ $(m^{-2}/s)$ at 800 $\degree$C')
ax.set_xlabel('log$_{10}$ IV-Al (a.p.f.u.)')

xytextloc = [(-1., -15.)] * 8
xytextloc[0] = (-2., -11.2) # Jaipur // a, c*
xytextloc[4] = (-1.9, -11.6) # Jaipur // b
xytextloc[1] = (-0.6, -12.35) # Nushan
xytextloc[2] = (-1.5, -10.95) # PMR
xytextloc[3] = (-1.2, -11.4) # Fuego


for idx in range(5):
    ax.plot(x[idx], bulk[idx], #label=Names[idx], 
            clip_on=False, **style[idx])
    xyloc = (x[idx], bulk[idx])
    label = Names[idx]
    ax.annotate(label, xy=xyloc, xytext=xytextloc[idx],
                arrowprops=dict(facecolor='black', arrowstyle='->'))

### Xenoliths
#xenoMin = np.log10(0.2 - 0.07)
#xenoMax = np.log10(0.2 + 0.07)
#xenoDrange = 0.5
#xenoDmin = min(bulk)
#fbox = mpatches.FancyBboxPatch([xenoMin, xenoDmin], xenoMax-xenoMin, xenoDrange,
#                               boxstyle=mpatches.BoxStyle("Round", pad=0.02),
#                               fill=True, alpha=0.25, color='r')
#ax.add_patch(fbox)
#ax.text(xenoMax, xenoDmin, 
#        'peridotite\nmantle\ncpx?', ha='right', va='bottom')


### Fuego
FuegoLoc = (np.log10(Al[3]), (bulk[0]+bulk[-1])/2.)
ax.add_artist(Ellipse(FuegoLoc, 0.06, 1., facecolor='none'))

ax.set_xlim(-2.0, -0.2) 
ax.grid(True)
ax.legend(loc=4, ncol=1, fontsize=10, fancybox=True)
plt.tight_layout()

#plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig12.eps', 
#            format='eps', dpi=1000)
#plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig12.tif', dpi=600)

# Add in Demouchy cpx in xenolith
xDem = np.log10(0.096)
yDem = -12.
plt.plot(xDem, yDem, 'xk', markeredgewidth=3)
ax.text(-1.3, -12.25, 'Demouchy et al. 2006\nxenolith cpx')

# D'Orazio Pali-Aiki phenocryst cpx
xPali = np.log10(np.array([0.14, 0.119, 0.094, 0.114]))
yPali = [-11.]*len(xPali)
plt.plot(xPali, yPali, '+r', markeredgewidth=3, label='Pali-Aiki pheno. core')

xPali = np.log10(np.array([0.267, 0.272, 0.144, 0.144]))
yPali = [-11.]*len(xPali)
plt.plot(xPali, yPali, 'og', markeredgewidth=1, label='Pali-Aiki pheno. rim',
         alpha=0.5)

xPali = np.log10(np.array([0.302, 0.244]))
yPali = [-11.]*len(xPali)
plt.plot(xPali, yPali, '^b', markeredgewidth=1, label='Pali-Aiki groundmass',
         alpha=0.5)

ax.legend(loc=3)
            
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Al_tetrahedral', dpi=300)

plt.show(fig)
print 'Finished'
