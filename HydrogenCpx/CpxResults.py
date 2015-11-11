# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 17:55:02 2015

@author: Ferriss

Before and after heating comparison of FTIR spectra on the rims of cpx samples
"""

from HydrogenCpx import my_spectra
import string
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator

#%%################ DATA SETUP ################################
prof_idx = 2
specs = []
ispecs = []

wblist = [
          my_spectra.K4wb_quench,
          my_spectra.K3wb_init,
          my_spectra.K3wb_700C_19hr,
          my_spectra.K3wb_700C_35hr,
          my_spectra.K3wb_800C_15hr,
          my_spectra.K3wb_6days,
          my_spectra.K4wb_1hr, 
          my_spectra.K4wb_91hr,
          my_spectra.K4wb_154hr, 
          my_spectra.K5wb_75hr,
          my_spectra.J1wb,
          ]
proflist = []
iproflist = []
for wb in wblist:
    wb.setupWB()
    proflist.append(wb.profiles[prof_idx])
    iproflist.append(wb.initial_profiles[prof_idx])

for prof in proflist + iproflist:
    prof.make_spectra_list()

for prof in proflist:
    specs.append(prof.spectra_list[-1])
    
for iprof in iproflist:
    ispecs.append(iprof.spectra_list[-1])  
  
  
specs.append(my_spectra.PMR_6)
ispecs.append(my_spectra.PMR_0)

for sp in specs + ispecs:
    sp.get_baseline()
    sp.subtract_baseline()

for k in range(len(specs)):
    print ispecs[k].fname
    print specs[k].fname
    print '\n'

#%% ####################### Figure #####################################
### setup and labeling
fig = plt.figure()
fig.set_size_inches(6., 7.5)

ncol = 3
nrow = 4
gs = gridspec.GridSpec(nrow, ncol)

axes = []
for row in range(nrow):
    for col in range(ncol):
        axes.append(plt.subplot(gs[row, col]))

xtickgrid=250
ytickgrid=0.5
xmajorLocator = MultipleLocator(xtickgrid)
ymajorLocator = MultipleLocator(ytickgrid)

#fig.suptitle('Baseline-subtracted FTIR spectra', fontsize=14)
fig.text(0.04, 0.5, 'Absorbance (cm$^{-1}$)', 
         ha='center', va='center', rotation='vertical', fontsize=14)
axes[10].set_xlabel('Wavenumber (cm$^{-1}$)', fontsize=14)

top = 1.4
legfont = 9
for ax in axes:
#    ax.xaxis.set_major_locator(xmajorLocator)
#    ax.yaxis.set_major_locator(ymajorLocator)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.set_xlim(3700, 3200)
    ax.set_ylim(0, top)

for idx in [0, 3, 6, 9]:
    plt.setp(axes[idx].get_yticklabels(), visible=True)

for idx in [9, 10, 11]:
    plt.setp(axes[idx].get_xticklabels(), visible=True, rotation=45)
fig.autofmt_xdate()

fig.subplots_adjust(top=0.98, bottom=0.1, right=0.98)

labels = [
          'K4\n480$\degree$C\n0.6hr',
          'K3, 696$\degree$C, 2hr',
          'K3\n696$\degree$C\n   19hr',
          'K3, 696$\degree$C, 35hr',
          'K3, 796$\degree$C, 15.7hr',
          'K3, 812$\degree$C, 140hr',
          'K4, 904$\degree$C, 0.7hr',
          'K4, 904$\degree$C, 91hr',
          'K4, 904$\degree$C, 154hr',
          'K5\n1000$\degree$C\n75hr',
          'Jaipur diopside\n    904$\degree$C, 0.6hr',
          'PMR, 800$\degree$C, 3hr\n    Divided by 5',
          ]

for idx in range(len(labels)):
    lab = '. '.join((string.ascii_uppercase[idx], labels[idx]))
    axes[idx].text(3660, top - (0.06*top), lab, va='top')
    
#fig.savefig('../../../CpxPaper/figures/CpxResults.png', dpi=100)

def labelpeaks(axes, idx, x, y1, y2, peakpos=[3645, 3540, 3460, 3355]):
    """Label peak at positions peakpos. x is wavenumber range.
    Label hits just above highest value in y1 and y2"""
    for peakwn in peakpos:        
        peak_idx = (np.abs(x-peakwn)).argmin()
        if y1[peak_idx] > y2[peak_idx]:
            ypeak = y1[peak_idx]
        else:
            ypeak = y2[peak_idx]
        ypeak = ypeak + 0.01
        xpeak = peakwn
        if peakwn == 3617:
            xtext = 3600
        elif peakwn == 3443:
            xtext = 3423
        else:
            xtext = xpeak
        axes[idx].annotate(str(peakwn), xy=(xpeak, ypeak), 
                       xytext=(xtext, ypeak+0.2), fontsize=8,
                        rotation=90, ha='center', va='bottom',
                        arrowprops=dict(facecolor='black', arrowstyle='->',
                                        linewidth=0.7))
#        axes[idx].text(xpeak, ypeak, ''.join(('$\leftarrow$', 
#                       '{:.0f}'.format(peakwn))), 
#                       rotation=90, ha='center', va='bottom', fontsize=8)

################# K3 and K4 pre-anneals ################################### 


for idx in [0, 1]:
    x  = specs[idx].base_wn
    y1 = ispecs[idx].abs_nobase_cm
    y2 = specs[idx].abs_nobase_cm
    axes[idx].plot(x, y1, label='Initial', color='g', linewidth=3,
        linestyle='--')
    axes[idx].plot(x, y2, label='Pre-\nanneal', color='b', linewidth=2)
    axes[idx].fill_between(x, y1, y2, where=y1<=y2, facecolor='b', 
            interpolate=True, alpha=0.2)
    axes[idx].fill_between(x, 0, y1, facecolor='y', alpha=0.1)
    labelpeaks(axes, idx, x, y1, y2, peakpos=[3645, 3617, 3540, 3460, 3443, 3355])

axes[0].legend(loc=1, fontsize=legfont)

# label peaks


################# Kunlun K3 and K4 pre-anneal and after ########################
for idx in np.arange(2, 9, 1):
    x  = specs[idx].base_wn
    y1 = ispecs[idx].abs_nobase_cm
    y2 = specs[idx].abs_nobase_cm
    axes[idx].plot(x, y1, label='Pre-\nanneal', color='b', linewidth=2)
    axes[idx].plot(x, y2, label='Final', color='r', linewidth=1)
    axes[idx].fill_between(x, y1, y2, where=y1>=y2, facecolor='b', 
            interpolate=True, alpha=0.2)
    axes[idx].fill_between(x, 0, y2, facecolor='r', alpha=0.1)
    labelpeaks(axes, idx, x, y1, y2, peakpos=[3645, 3617, 3540, 3460, 3443, 3355])
axes[2].legend(loc=1, fontsize=legfont)


######################### J1 and K5 initial and final ######################## 
for idx in np.arange(9, 11, 1):
    x  = specs[idx].base_wn
    y1 = ispecs[idx].abs_nobase_cm
    y2 = specs[idx].abs_nobase_cm
    axes[idx].plot(x, y1, label='Initial', color='g', linewidth=3,
                linestyle='--')
    axes[idx].plot(x, y2, label='Final', color='r', linewidth=1)
    axes[idx].fill_between(x, 0, y2, facecolor='r', alpha=0.1)
    axes[idx].fill_between(x, y1, y2, where=y1>=y2, facecolor='y', 
            interpolate=True, alpha=0.1)
    labelpeaks(axes, idx, x, y1, y2, peakpos=[3645, 3617, 3540, 3460, 3443, 3355])
axes[9].legend(loc=1, fontsize=legfont)

###################### PMR spectra divided by 5 ##############################
for idx in [11]:
    x  = specs[idx].base_wn
    y1 = ispecs[idx].abs_nobase_cm / 5.
    y2 = specs[idx].abs_nobase_cm / 5.
    axes[idx].plot(x, y1, label='Initial', color='g', linewidth=3,
                linestyle='--')
    axes[idx].plot(x, y2, label='Final', color='r', linewidth=1)
    axes[idx].fill_between(x, 0, y2, facecolor='r', alpha=0.1)
    axes[idx].fill_between(x, y1, y2, where=y1>=y2, facecolor='y', 
            interpolate=True, alpha=0.1)
    labelpeaks(axes, idx, x, y1, y2, peakpos=[3620, 3550, 3460, 3355])


#plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig4.eps', 
#            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig4.tif', dpi=600)
#plt.savefig('Fig4.eps', format='eps', dpi=1000)
