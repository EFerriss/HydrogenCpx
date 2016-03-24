# -*- coding: utf-8 -*-
"""
CreateJd on Thu Jul 02 08:37:08 2015

@author: Ferriss

Arrhenius plot of Jaipur, PMR, and Fe-bearing olivine
"""

import pynams.diffusion as diff
import numpy as np
import pynams.diffusivity_library as dlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import cpx_spectra
plt.style.use('paper')


# create list of diffusivity groups to plot
Kunlun = cpx_spectra.K_bulk
Jaipur = dlib.H_diopside_Woods00
my_Jaipur = cpx_spectra.J_bulk
KM98_fast = dlib.KM98_fast
PMR = cpx_spectra.D_PMR
Xia = dlib.H_cpxBasanite_Xia00
KM98_slow = dlib.KM98_slow
Fo = dlib.DM03
pnav_Mg = dlib.pnav_Mg
pnav_Si = dlib.pnav_Si
pnav_Ti = dlib.pnav_Ti
Sundvall = dlib.H_diopside_Sundvall
Ingrin = dlib.H_CrDiopside_Ingrin95

Dlist = [
         KM98_fast, 
         Fo, 
         pnav_Mg, 
         pnav_Ti,
         pnav_Si,
         Ingrin, 
         Sundvall,
         Kunlun,
         Jaipur, 
         my_Jaipur, 
         PMR, 
         Xia,          
         ]

#%%## Plotting

# setup figure
top = -8
sunk = -1.5
fig, ax, hleg = diff.Arrhenius_outline(low=7., high=10.5, 
                                       bottom=-17, top=top,
                                       figsize_inches = (6.5, 5), 
                                       shrinker_for_legend = 0.,
                                       generic_legend=True, sunk=sunk, ncol=3)
ax.grid('off')

# Plotting styles
Kunlun.basestyle = {'marker' : 'o', 'markersize' : 14, 'color': 'chocolate', 'alpha':1.}

Jaipur.description = 'Jaipur\ncpx'
Jaipur.basestyle = {'marker' : 'o', 'markersize' : 8, 'color': 'orange', 'alpha':0.5}

my_Jaipur.description = 'Jaipur di. (this work)'
my_Jaipur.basestyle = {'marker' : 'o', 'markersize' : 14, 'color': 'orange', 'alpha':1}

PMR.description = 'PMR\ncpx'
PMR.basestyle = {'marker' : 'o', 'markersize' : 20, 'color': 'darkred', 'alpha':1, 'mew':3}

Xia.description = unicode('N\374shan cpx', 'latin-1')
Xia.basestyle = {'marker' : 'o', 'markersize' : 8, 'color': 'yellow', 'alpha':1}

Ingrin.description = 'Russian cpx'
Ingrin.basestyle = {'marker' : 'o', 'markersize' : 8, 'color': 'salmon', 'alpha':1}

Sundvall.description = 'synthetic cpx'
Sundvall.basestyle = {'marker' : 'o', 'markersize' : 8, 'color': 'firebrick', 'alpha':1}

KM98_fast.description = "ol.\n[p.p.]"
KM98_fast.basestyle = {'marker' : 's', 'markersize' : 8, 'color': 'blue', 'alpha':0.5}

Fo.description = 'forsterite'
Fo.basestyle = {'marker' : 's', 'markersize' : 7, 'color': 'sage'}

pnav_Mg.description = 'ol. [Mg]'
pnav_Mg.basestyle = {'marker' : 's', 'color': 'green', 'markersize' : 6, 'mew' : 1.5}

pnav_Ti.description = 'ol. [Ti]'
pnav_Ti.basestyle = {'marker' : 's', 'color': 'darkcyan', 'markersize' : 6, 'mew' : 3.}

pnav_Si.description = 'ol. [Si]'
pnav_Si.basestyle = {'marker' : 's', 'color': 'slategray', 'markersize' : 7, 'mew' : 1.}

# Fe-rich cpx fill
rotate_main = -14.5
ax.text(8.3, -9.3, 'cpx with Fe > 0.05 a.p.f.u.', rotation=rotate_main)
m = -0.875
b1 = -2.2
gap = 2.
p1 = (m, b1)
p2 = (m, b1-gap)
x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 50)
y1 = np.polyval(p1, x)
y2 = np.polyval(p2, x)
plt.fill_between(x, y1, y2, facecolor='grey', alpha=0.2, interpolate=True)

# apply labels
for D in Dlist:
    D.xyloc = [7., -17.]
    D.xytextloc = [8., -18.]
    D.ha = 'left'
    D.rotate = 0

PMR.xyloc = [9.35, -11.]
PMR.xytextloc = [9.45, -11.4]
Jaipur.xyloc = [9.8, -11.25]
Jaipur.xytextloc = [9.95, -11.7]
KM98_fast.xyloc = [7.8, -9.4]
KM98_fast.xytextloc = [7.9, -9.9]
Xia.xyloc = [8.175, -11.3]
Xia.xytextloc = [7.9, -10.8]
Xia.rotate = rotate_main
Ingrin.xyloc = [9.35, -13.]
Ingrin.xytextloc = [8.7, -12.8]
Ingrin.rotate = rotate_main + 5.
Kunlun.xyloc = [9.2, -13.65]
Kunlun.xytextloc = [9.3, -13.7]
Kunlun.rotate = 0.
pnav_Mg.xyloc = [9.2, -14.3]
pnav_Mg.xytextloc = [9.4, -14.6]
pnav_Ti.xyloc = [9.2, -15.3]
pnav_Ti.xytextloc = [9.4, -15.5]
Sundvall.xyloc = [9.2, -15.9]
Sundvall.xytextloc = [9.4, -16.3]
pnav_Si.xyloc = [7.9, -15.7]
pnav_Si.xytextloc = [7.9, -16.]
Fo.xyloc = [7.5, -11.5]
Fo.xytextloc = [7.1, -11.2]

Dlist_just_text = [PMR, KM98_fast, Kunlun, pnav_Mg, Xia,
                  pnav_Ti, Sundvall, pnav_Si, Jaipur, 
                  Fo, Ingrin]   

for D in Dlist_just_text:    
    label = D.description
    try:
        ax.text(D.xytextloc[0], D.xytextloc[1], label, ha=D.ha, rotation=D.rotate)
    except AttributeError:
        D.ha = 'center'
        D.rotate = 0.
        ax.text(D.xytextloc[0], D.xytextloc[1], label, ha=D.ha, rotation=D.rotate)
    
# add lines and data
def ArrheniusLine(Ea, D0, figax, celsiusMin=700., celsiusMax=1000., scale=1e4,
                  style=None):
    GAS_CONSTANT = 0.00831 # kJ/mol K
    Tmin = scale / (celsiusMin + 273.15)
    Tmax = scale / (celsiusMax + 273.15)
    Dmin = D0 * np.exp(Ea / (GAS_CONSTANT * Tmin))
    Dmax = D0 * np.exp(Ea / (GAS_CONSTANT * Tmax))
    figax.plot([Tmin, Tmax], [Dmin, Dmax], **style)

for D in Dlist:
    D.plotD(ax, show_error=True, ecolor='k')   
    style = D.basestyle.copy()

# add Fuego
FuegoLoc = (1E4/(1030.+273.15), -9.5)
Fuego = Ellipse(FuegoLoc, 0.1, 1.1, facecolor='none', edgecolor='r')
ax.add_artist(Fuego)
ax.text(7.6, -9.45, 'Fuego\ncpx', ha='right')

# Make second legend with orientation information
para, = ax.plot(0, 0, 'sk', alpha=0.5, fillstyle='left', label='|| a or a*', linestyle='--')
parb, = ax.plot(0, 0, 'sk', alpha=0.5, fillstyle='bottom', label='|| b', linestyle='-.')
parc, = ax.plot(0, 0, 'sk', alpha=0.5, fillstyle='right', label='|| c or c*', linestyle=':')
paru, = ax.plot(0, 0, 'sk', alpha=0.5, fillstyle='none', label='unoriented', linestyle='-')
plt.legend(handles=[para, parb, parc, paru], loc=1)

plt.savefig('Fig11_Arrhenius.eps', format='eps', dpi=1000)
fig.savefig('Fig11_Arrhenius.tif', format='tif', dpi=300)

plt.show(fig)
print 'Finished'
