# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 08:37:08 2015

@author: Ferriss

Arrhenius plot of various diffusivities in clinopyroxene
"""

import pynams.diffusion as diff
import numpy as np
import matplotlib.pyplot as plt
from HydrogenCpx import my_spectra
plt.style.use('paper')

K_bulk = my_spectra.K_bulk
K_3645 = my_spectra.K_3645
K_3617 = my_spectra.K_3617
K_3540 = my_spectra.K_3540
K_3443 = my_spectra.K_3443
K_3350 = my_spectra.K_3350

### Plotting
sunk = 1.25
fig, ax, hleg = diff.Arrhenius_outline(low=7.5, high=9.5, bottom=-16, top=-11.,
                      figsize_inches = (6.5, 5), shrinker_for_legend = 0.,
                      generic_legend=True, sunk=sunk, ncol=3)


offsets = [0, 10, 5, 0, -5, -10]
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
    if style['color'] == 'yellow':
        style['color'] = 'teal'
    ecolor = style['color']

    #### Solve for each activation energy Ea in kJ/mol and D0 in m2/s ###
    Ea, D0 = D.solve_Ea_D0(orient='x')
    # print out for excel
    print '{:.3f}'.format(Ea.n), '{:.3f}'.format(Ea.s), \
          '{:.3f}'.format(np.log10(D0.n)), '{:.3f}'.format(np.log10(D0.s))
    # and plot!
    x_lines = np.array(ax.get_xlim())
    celsius_lines = (1E4 / x_lines) - 273.15
    y_lines = [0., 0.]
    y_lines[0] = diff.whatIsD(Ea.n, D0.n, celsius_lines[0])
    y_lines[1] = diff.whatIsD(Ea.n, D0.n, celsius_lines[1])
    
    if D == K_bulk:
        lw = 3
    else:
        lw = 1.5
    ax.plot(x_lines, y_lines, '-', color=ecolor, alpha=0.8, linewidth=lw)
    ### Annotate lines
    xa = 1E4 / (celsius_annotate[idx_offset] + 273.15)
    ya = diff.whatIsD(Ea.n, D0.n, celsius_annotate[idx_offset])
    ax.annotate(D.description, xy=(xa, ya), 
                xytext=(xatext[idx_offset], yatext[idx_offset]),
                arrowprops=dict(facecolor='red', arrowstyle='->',))
#                                 relpos=(0, 1), color='k'))
    ### Plot data points ###

    if offsets[idx_offset] != 0:
        D.description = ''.join((D.description, ', offset ', 
                             str(offsets[idx_offset]) ,' $\degree$C'))
    D.add_to_legend(ax, hleg, ncol=2, sunk=sunk)

    D.plotDx(ax, style_x=style, er=True, offset_celsius=offsets[idx_offset], 
             ecolor=ecolor, plotline=False)

    idx_offset = idx_offset + 1

#K_bulk.plotDu(ax, er=True, offset_celsius=0)
#K_3645.plotDu(ax, er=True, offset_celsius=5)

ax.text(7.6, -11.5, '$Kunlun$ $diopside$', fontsize=16, backgroundcolor='w')

#plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig7.eps', 
#            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig7.tif', dpi=600)

plt.show(fig)
print 'Finished'
