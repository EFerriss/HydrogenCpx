# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 14:37:34 2015

@author: Ferriss
"""
import pynams.pynams as nams
import pynams.styles as st
import cpx_spectra
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('paper')

diff = nams.diffusion

PMR_unpol = cpx_spectra.PMR_unpol

PMR_unpol.make_spectra_list()
PMR_unpol.initial_profile = PMR_unpol
PMR_unpol.get_baselines()
PMR_unpol.make_area_list()
PMR_unpol.get_peakfit()

#%% ##### Make figure #########
specs = PMR_unpol.initial_profile.spectra_list

fig = plt.figure(figsize=(6.5, 3.))
ax = fig.add_subplot(121)

specs[0].plot_spectrum_outline(figaxis=ax)

colors = ['k', 'magenta', 'purple', 'blue', 'green', 'orangered', 'r']
times = ['initial', '15min', '30min', '45min', '1hr', '2hr', '3hr']
xtext = [3400]*7

x_arrow = [3435, 3435, 3435, 3435, 3560, 3550, 3550]
x_text = [3380.]*4 + [3635]*3
y_text = [2.5, 2.0, 1.5, 1.2, 0.6, 0.4, 0.2]

for idx, spec in enumerate(specs):
    style = st.style_1.copy()
    style['color'] = colors[idx]
    spec.plot_subtractbaseline(style=style, figaxis=ax, label=times[idx])
    
    arrow_idx = (np.abs(spec.base_wn - x_arrow[idx])).argmin()
    y_arrow = spec.abs_nobase_cm[arrow_idx]
    ax.annotate(times[idx], xy=(x_arrow[idx], y_arrow), 
                xytext=(x_text[idx], y_text[idx]),
                arrowprops=dict(edgecolor=colors[idx], arrowstyle='->'))

plt.grid(False)    
ax.set_ylim(0, 3.)
ax.set_title('')
plt.tight_layout()

ax.text(3680, 2.7, 'A.', fontsize=11,)

######### Diffusion side of figure ###########
ax = fig.add_subplot(122)

D = -11.02
Drange = 0.2

Ds = [-11.5, -11., -10.5]
sty = {'color' : 'k', 'linestyle' : 'None', 'marker' : 'o', 'alpha' : 0.5}

xmax = 4.
PMR_unpol.plot_timeseries(D_list=Ds, style=sty, max_hours=xmax, figaxis=ax)

ax.text(0.3, 1.3, 'B.', fontsize=11)

#ax.text(0.2, 1.25, 'B. $augite$\n$PMR-53$ \n $800$ $\degree C$', 
#        ha='left', fontsize=11, va='center')

xloc =  3.9
ylocs = [0.5, 0.2, 0.02]

for k in range(len(Ds)):
    ax.text(xloc, ylocs[k], Ds[k], ha='right')

plt.tight_layout()

# bulk H
PMR_unpol.make_area_list()
C = np.array(PMR_unpol.areas_list)
C0 = PMR_unpol.areas_list[0]
y = C / C0
ax.plot(PMR_unpol.times_hours, y, 'o', color='k', alpha=0.5,
        clip_on=False, markersize=8, label='bulk H')

# peak at 3620
C = PMR_unpol.peak_areas[0, :]
C_0 = C[0]
y = C / C_0

ax.plot(PMR_unpol.times_hours, y, 
        '+', color='orange', clip_on=False, markersize=10, 
        label='3620 cm$^{-1}$', alpha=1, mew=2)

# peak at 3550
C = PMR_unpol.peak_areas[1, :]
C_0 = C[0]
y = C / C_0
ax.plot(PMR_unpol.times_hours, y, 
        's', color='teal', clip_on=False, markersize=8, 
        label='3540 cm$^{-1}$', alpha=0.5, mew=1)

# peak at 3460 cm-1
C = PMR_unpol.peak_areas[2, :]
C_0 = C[0]
y = C / C_0
ax.plot(PMR_unpol.times_hours, y, color='g', clip_on=False, markersize=10, 
        label='3460 cm$^{-1}$', alpha=1, mew=2, marker='3', linestyle='none')

# peak at 3355 cm-1
C = PMR_unpol.peak_areas[3, :]
C_0 = C[0]
y = C / C_0
ax.plot(PMR_unpol.times_hours, y, 'd', color='violet', 
        clip_on=False, markersize=8, 
        label='3355 cm$^{-1}$', alpha=0.7)

ax.legend(loc=1, ncol=1, fontsize=9)
ax.set_xlim(0, xmax)
ax.set_ylim(0, 1.5)

fig.autofmt_xdate()

plt.savefig('Fig9_augitePMR.eps', format='eps', dpi=1000)
fig.savefig('Fig9_augitePMR.tif', format='tif', dpi=300)
