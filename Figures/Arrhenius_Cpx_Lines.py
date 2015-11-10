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
import CpxCode.my_spectra as my_spectra
plt.style.use('paper')

Jaipur = dlib.H_diopside_Woods00

my_Jaipur = my_spectra.D_Jaipur_bulkH

Jaipur.celsius_x = Jaipur.celsius_x + my_Jaipur.celsius_all
Jaipur.celsius_y = Jaipur.celsius_y + my_Jaipur.celsius_all
Jaipur.celsius_z = Jaipur.celsius_z + my_Jaipur.celsius_all

Jaipur.logDx = np.append(Jaipur.logDx, my_Jaipur.logDx)
Jaipur.logDy = np.append(Jaipur.logDy, my_Jaipur.logDy)
Jaipur.logDz = np.append(Jaipur.logDz, my_Jaipur.logDz)

KM98_fast = dlib.KM98_fast
PMR = my_spectra.D_PMR
Xia = dlib.H_cpxBasanite_Xia00
KM98_slow = dlib.KM98_slow
Fo = dlib.DM03
pnav_Mg = dlib.pnav_Mg
pnav_Si = dlib.pnav_Si
pnav_Ti = dlib.pnav_Ti
Kunlun = my_spectra.K_bulk
Sundvall = dlib.H_diopside_Sundvall
Ingrin = dlib.H_CrDiopside_Ingrin95

#%%## Plotting
top = -8
sunk = -1.5
fig, ax, hleg = diff.Arrhenius_outline(low=7., high=10.5, 
                                       bottom=-17, top=top,
                      figsize_inches = (6.5, 5), shrinker_for_legend = 0.,
                      generic_legend=True, sunk=sunk, ncol=3)



#PMR.basestyle['alpha'] = 0.
Kunlun.basestyle['markerfacecolor'] = 'red'
Kunlun.basestyle['alpha'] = 1
Kunlun.description = 'Kunlun cpx'
PMR.description = 'PMR\ncpx'
PMR.basestyle['fillstyle'] = 'none'
PMR.basestyle['markeredgewidth'] = 2
PMR.basestyle['markeredgecolor'] = 'r'
Jaipur.description = 'Jaipur\ncpx'
my_Jaipur.description = 'Jaipur di. (this work)'
#Xia.description = 'High-Al cpx (Xia et al. 2000)'
Xia.description = unicode('N\374shan cpx', 'latin-1')
Xia.basestyle['markerfacecolor'] = 'purple'
Xia.basestyle['marker'] = 'h'
Xia.basestyle['markersize'] = 10
#Fo.description = '[Mg] mechanism in olivine\n(Demouchy & Mackwell, 2003'
Fo.description = 'forsterite'
Fo.basestyle['marker'] = 'p'
#pnav_Mg.description = '[Mg] mechanism in olivine\n(Padron-Navarta et al. 2014)'
pnav_Mg.description = 'ol. [Mg]'
pnav_Si.description = 'ol. [Si]'
pnav_Ti.description = 'ol. [Ti]'
pnav_Mg.basestyle['marker'] = 'p'
pnav_Si.basestyle['marker'] = '+'
pnav_Ti.basestyle['marker'] = 'x'
pnav_Si.basestyle['color'] = 'darkcyan'
#pnav_Si.description = '[Si] mechanism in olivine\n(Padron-Navarta et al. 2014)'
#pnav_Ti.description = '[Ti] mechanism in olivine\n(Padron-Navarta et al. 2014)'
#KM98_fast.description = 'proton-polaron in olivine\n(Kohlstedt & Mackwell, 1998'
KM98_fast.description = "ol.\n[p.p.]"
KM98_fast.basestyle['markerfacecolor'] = 'sienna'
KM98_fast.basestyle['alpha'] = 0.7
KM98_slow.description = '[Mg] mechanism in olivine\n(Kohlstedt & Mackwell, 1998'
#Sundvall.description = 'synthetic Fe-poor diopside\n(Sundvall et al., 2009)'
Sundvall.description = 'synthetic cpx'
Sundvall.basestyle['markersize'] = 10
#Ingrin.description = 'High-Cr diopside\n(Ingrin et al. 1995)'
Ingrin.description = 'Russian cpx'
Ingrin.basestyle['markersize'] = 8
Ingrin.basestyle['color'] = 'k'


#ax.text(10.4, top+0.04*top, 'H diffusivities in olivine and\nFe-bearing clinopyroxene', 
#        ha='right', va='top', backgroundcolor='white')
rotate_main = -14.5
ax.text(8.3, -9.3, 'cpx with Fe > 0.05 a.p.f.u.', rotation=rotate_main)

### Fe-rich cpx fill
m = -0.875
b1 = -2.2
gap = 2.
p1 = (m, b1)
p2 = (m, b1-gap)
x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 50)
y1 = np.polyval(p1, x)
y2 = np.polyval(p2, x)
plt.fill_between(x, y1, y2, facecolor='grey', alpha=0.2, interpolate=True)

Dlist = [
         Kunlun,
         Jaipur, 
         my_Jaipur, 
         PMR, 
         Xia, 
         Ingrin, 
         Sundvall,
         KM98_fast, 
         Fo, 
         pnav_Mg, 
         pnav_Ti,
         pnav_Si
         ]

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
#Ingrin.xytextloc = [9.4, -13.3]
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

Dlist_to_label = []
Dlist_just_text = [PMR, KM98_fast, Kunlun, pnav_Mg, Xia,
                  pnav_Ti, Sundvall, pnav_Si, Jaipur, Fo, Ingrin]   
idx = 0          
for D in Dlist_to_label:
    label = D.description
    ax.annotate(label, xy=D.xyloc, xytext=D.xytextloc, ha=D.ha,
            arrowprops=dict(facecolor='black', arrowstyle='->',),
            rotation=D.rotate)
    idx = idx + 1
    
for D in Dlist_just_text:    
    label = D.description
    ax.text(D.xytextloc[0], D.xytextloc[1], label, ha=D.ha, rotation=D.rotate)
    
def ArrheniusLine(Ea, D0, figax, celsiusMin=700., celsiusMax=1000., scale=1e4,
                  style=None):
    GAS_CONSTANT = 0.00831 # kJ/mol K
    Tmin = scale / (celsiusMin + 273.15)
    Tmax = scale / (celsiusMax + 273.15)
    Dmin = D0 * np.exp(Ea / (GAS_CONSTANT * Tmin))
    Dmax = D0 * np.exp(Ea / (GAS_CONSTANT * Tmax))
    figax.plot([Tmin, Tmax], [Dmin, Dmax], **style)

for D in Dlist:
#    D.add_to_legend(ax, hleg, ncol=2, sunk=sunk)
    D.plotD(ax, er=True, ecolor='k')
    
    style = D.basestyle.copy()

### Fuego
FuegoLoc = (1E4/(1030.+273.15), -9.5)
ax.add_artist(Ellipse(FuegoLoc, 0.1, 1.1, facecolor='none'))
ax.text(7.6, -9.45, 'Fuego\ncpx', ha='right')
#ax.annotate('Fuego\ncpx', xy=FuegoLoc, xytext=(7.55, -9.4), 
#            ha='right',
#            arrowprops=dict(facecolor='black', arrowstyle='->'))

### Make second legend with orientation information
para, = ax.plot(0, 0, 'sk', alpha=0.5, fillstyle='left', label='|| a or a*', linestyle='--')
parb, = ax.plot(0, 0, 'sk', alpha=0.5, fillstyle='bottom', label='|| b', linestyle='-.')
parc, = ax.plot(0, 0, 'sk', alpha=0.5, fillstyle='right', label='|| c or c*', linestyle=':')
paru, = ax.plot(0, 0, 'sk', alpha=0.5, fillstyle='none', label='unoriented', linestyle='-')
plt.legend(handles=[para, parb, parc, paru], loc=1)

ax.grid('off')

plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig11.eps', 
            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig11.png', dpi=300)


plt.show(fig)
print 'Finished'
