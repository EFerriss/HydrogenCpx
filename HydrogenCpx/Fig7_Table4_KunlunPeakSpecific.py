# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 08:37:08 2015

@author: Ferriss

Arrhenius plot of various diffusivities in clinopyroxene
"""

from pynams import diffusion as diff
from uncertainties import ufloat
import numpy as np
import matplotlib.pyplot as plt
import cpx_spectra
plt.style.use('paper')

reload(diff)
reload(cpx_spectra)

K_bulk = cpx_spectra.K_bulk
K_3645 = cpx_spectra.K_3645
K_3617 = cpx_spectra.K_3617
K_3540 = cpx_spectra.K_3540
K_3443 = cpx_spectra.K_3443
K_3350 = cpx_spectra.K_3350

### Plotting
sunk = 0.9
fig, ax, hleg = diff.Arrhenius_outline(low=7.5, high=9.5, bottom=-16, top=-11.,
                      figsize_inches = (6.5, 5), shrinker_for_legend = 0.,
                      generic_legend=True, sunk=sunk, ncol=3)


#offsets = [0, 10, 5, 0, -5, -10]
offsets = [0]*6
alphas = [0.5, 1., 1., 0.7, 1., 0.7]
celsius_annotate = [925., 975., 850, 825., 950., 875]
xatext = [8.2, 7.7, 8.85, 9.1, 7.95, 8.6]
yatext = [-13.9, -14.5, -11.9, -11.6, -14.2, -12.4, -12.7]
idx_offset = 0
Ea_list = []
D0_list = []

for D in [
          K_bulk, 
          K_3645, 
          K_3617, 
          K_3540, 
          K_3443, 
          K_3350
          ]:
    style = D.basestyle
    style['fillstyle'] = 'full'
    style['alpha'] = alphas[idx_offset]
    ecolor = style['color']

    #### Solve for each activation energy Ea in kJ/mol and D0 in m2/s ###
    Ea, D0 = D.solve_Ea_D0(orient='x')
    
    # Fit also to highest and lowest possible diffusivities
    logD_list_high = D.logD[0] + D.logD_error[0]
    Ea_high, D0_high = diff.solve_Ea_D0(logD_list_high, D.celsius_all)

    logD_list_low = D.logD[0] - D.logD_error[0]
    Ea_low, D0_low = diff.solve_Ea_D0(logD_list_low, D.celsius_all)

    Ea_mean = np.mean(np.array([Ea.n, Ea_high.n, Ea_low.n]))
    D0_mean = np.mean(np.array([D0.n, D0_high.n, D0_low.n]))

    Ea_std = np.std(np.array([Ea.n, Ea_high.n, Ea_low.n]))
    D0_std = np.std(np.array([D0.n, D0_high.n, D0_low.n]))

    # go with the higher number for the error
    if Ea_std > Ea.s:
        Ea = ufloat(Ea.n, Ea_std)
    
    if abs(D0_std) > abs(D0.s):
        D0 = ufloat(D0.n, D0_std)

    # print out activation energies for Table
#    print
#    print D.description
    print '{:.3f}'.format(Ea.n), '{:.3f}'.format(Ea.s), \
          '{:.3f}'.format(np.log10(D0.n)), '{:.3f}'.format(np.log10(D0.s))   
    
    # and plot!
    x_lines = np.array(ax.get_xlim())
    celsius_lines = (1E4 / x_lines) - 273.15
    y_lines = [0., 0.]
    y_lines[0] = diff.whatIsD(Ea.n, D0.n, celsius_lines[0], printout=False)
    y_lines[1] = diff.whatIsD(Ea.n, D0.n, celsius_lines[1], printout=False)
    
    if D == K_bulk:
        lw = 3
    else:
        lw = 1.5
    ax.plot(x_lines, y_lines, '-', color=ecolor, alpha=0.8, linewidth=lw)
    ### Annotate lines
    xa = 1E4 / (celsius_annotate[idx_offset] + 273.15)
    ya = diff.whatIsD(Ea.n, D0.n, celsius_annotate[idx_offset], printout=False)
    ax.annotate(D.description, xy=(xa, ya), 
                xytext=(xatext[idx_offset], yatext[idx_offset]),
                arrowprops=dict(facecolor='red', arrowstyle='->',))

    ### Plot data points ###

    if offsets[idx_offset] != 0:
        D.description = ''.join((D.description, ', offset ', 
                             str(offsets[idx_offset]) ,' $\degree$C'))
    D.add_to_legend(ax, hleg, ncol=3, sunk=sunk)

    D.plotD(ax, style=style, show_error=True, 
            offset_celsius=offsets[idx_offset], 
            ecolor=ecolor, plotline=False)

    idx_offset = idx_offset + 1

ax.text(7.6, -11.5, '$Kunlun$ $diopside$', fontsize=16, backgroundcolor='w')

plt.savefig('Fig7_Arrhenius_KunlunPeakSpecific.eps', format='eps', dpi=1000)
plt.savefig('Fig7_Arrhenius_KunlunPeakSpecific.tif', format='tif', dpi=300)

plt.show(fig)
print 'Finished'
