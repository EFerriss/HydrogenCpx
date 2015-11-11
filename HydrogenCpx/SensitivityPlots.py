# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 13:06:23 2015

@author: Ferriss

Plot and play with variations in best-fit diffusivity with initial concentration 
"""
import my_spectra
import json
import numpy as np
import matplotlib.pyplot as plt
import pynams.diffusion as diffusion
from lmfit import minimize
import string

plt.style.use('paper')

func2min = diffusion.diffusion3Dwb_params

wbsK3 = [
#        my_spectra.K3wb_trueInit,
         my_spectra.K3wb_init,
         my_spectra.K3wb_800C_15hr,
         my_spectra.K3wb_6days,
         my_spectra.K4wb_quench,
         my_spectra.K4wb_1hr,
         ]

wbsK4_91 = [
           my_spectra.K4wb_init,
           my_spectra.K3wb_init,
           my_spectra.K3wb_800C_15hr,
           my_spectra.K4wb_quench,
           my_spectra.K4wb_1hr,
           my_spectra.K4wb_91hr
           ]

wbsK4_154 = [
       my_spectra.K4wb_init,
       my_spectra.K3wb_init,
       my_spectra.K3wb_800C_15hr,
       my_spectra.K4wb_quench,
       my_spectra.K4wb_1hr,
       my_spectra.K4wb_91hr,
       my_spectra.K4wb_154hr,
       ]

wbsK5 = [
       my_spectra.K5wb_init,
       my_spectra.K3wb_init,
       my_spectra.K3wb_800C_15hr,
       my_spectra.K4wb_quench,
       my_spectra.K4wb_1hr,
       my_spectra.K5wb_75hr
       ]

wbsJ = [
       my_spectra.J1wb_initial,
       my_spectra.J1wb
       ]

wb_list = wbsK3 + wbsK4_154 + wbsK5 + wbsJ + wbsK4_91

fname_list = ['K4_91', 'K3', 'K5', 'K4_154', 'J']
peak_idx_list = [0, 1, 2, 4, 5]
peak_label_list = ['3645', '3617', '3540', '3460', '3443', '3350']

for wb in wb_list:
    wb.get_baselines()
    wb.get_peakfit()
    wb.setupWB()
    for prof in wb.profiles:
        prof.make_wholeblock(True, True)
        prof.wb_waters = np.array(prof.wb_areas)
        
#%%
def solve_bestd(wb, peak_idx, initial, Dguess=-12, wholeblock=False,
                heights_instead=False):
    """Returns best-fit initial and isotropic D """
    L3 = []
    D3 = []
    for k in range(3):
        L3.append(wb.profiles[k].len_microns)
        D3.append(-12)

    x, y = wb.xy_picker(peak_idx=peak_idx, wholeblock=wholeblock, 
                           heights_instead=heights_instead)    
    for k in range(3):
       if sum(y[k]) < 0.5:
           y[k] = []
           x[k] = []

    params = diffusion.params_setup3D(L3, D3, wb.time_seconds, 
                                     initial=initial, isotropic=True)
    minimize(func2min, params, args=(x, y), 
                  kws={'raypaths' : wb.raypaths, 'show_plot' : False})
                  
    bestD = params['log10Dx'].value
    return bestD

def get_initialD(fname, peak_idx, wholeblock=False):
    """Takes name = 'K4_91', 'K3', 'K5', 'K4_154', or 'J' and peak index
    Returns (1) best-fit initial with (2) best-fit isotropic diffusivity 
    (3) list of all initials held constant and (4) resulting best-fit diffusivities
    Data generated in Sentivity.py    
    """
    if fname not in fname_list:
        print 'fname must be one of the following:'
        print fname_list
        return False, False, False, False
    fi = '../../../CpxPaper/figures/sensitivity_'
    if peak_idx is not None:
        peakbit = peak_label_list[peak_idx]
    else:
        if wholeblock is True:
            peakbit = 'bulk_wb'
        else:
            peakbit = 'bulk'
    filename = ''.join((fi, fname, '_', peakbit, '.txt'))
    f = open(filename, 'r')
    x = json.load(f)
    besti = x[0][0]
    bestD = x[1][0]
    i_list = x[0][1:]
    D_list = x[1][1:]
    return besti, bestD, i_list, D_list

def save_sensitivity_fig(fname, peak_idx, wholeblock=False):
    """Save sensitivity_fname_peaklabel.png"""
    if fname not in fname_list:
        print 'fname must be one of the following:'
        print fname_list
        return False, False        
    fi = '../../../CpxPaper/figures/sensitivity_'
    fend = '.png'
    if peak_idx is not None:
        peakbit = peak_label_list[peak_idx]
    else:
        if wholeblock is True:
            peakbit = 'bulk_wb'
        else:
            peakbit = 'bulk'
    fname = ''.join((fi, fname, '_', peakbit, fend))
    fig.savefig(fname)

def get_max(wb, peak_idx, wholeblock=False, initials=False):
    """return maximum peak area for given whole-block and peak index"""
    ms = []
    if initials is False:
        proflist = wb.profiles
    else:
        proflist = wb.initial_profiles
    for prof in proflist:
        if peak_idx is not None:
            ms.append(max(prof.peak_areas[peak_idx]))
        else:
            if wholeblock is True:
                ms.append(max(prof.wb_areas))
            else:
                ms.append(max(prof.areas_list))
    m = max(ms)
    return m

def get_max_global(wbs, peak_idx, wholeblock=False, initials=False):
    """Returns max peak area for set of whole-blocks"""
    mss = []
    for wb in wbs:
        mss.append(get_max(wb, peak_idx, wholeblock, initials))
    m = max(mss)
    return m

def plot_D(i_list, D_list, wb, peak_idx, besti=None, bestD=None, savefig=True,
           iused=None, show_max=True, wholeblock=False, show_max_init=True,
           maxwbs=None, boundaries=[], lloc=4, makefigure=True, ax=None):
    """Plots initial versus best-fit D. Optional: best initial + best diffusivity. Returns figure, axes."""
    if makefigure is True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        if ax is None:
            print 'Need axis handle; ax='
            return False

    ax.set_xlabel('initial area (cm$^{-2}$)')
    if peak_idx is not None:
        ax.set_ylabel(''.join(('best-fit isotropic log$_{10}$D\nto peak at ',
                          peak_label_list[peak_idx], ' cm$^{-1}$ (m$^2$/s)')))
    else:
        ax.set_ylabel('best-fit isotropic bulk log$_{10}$D$_x$')
                          
    ax.set_title(wb.name, fontsize=14)
    ax.plot(i_list, D_list, 'ok', markersize=6, label='initial held\nconstant')
#    ax.set_ylim(bottom, top)
#    ax.set_xlim(low, high)
#       ax.plot(ax.get_xlim(), [bestD+Dgap, bestD+Dgap], '-k')
#       ax.plot(ax.get_xlim(), [bestD-Dgap, bestD-Dgap], '-k')
    if iused is not None:
        ax.plot([iused, iused], ax.get_ylim(), '-b', linewidth=3,
                label='assumed\ninitial area')

    if (besti is not None) and (bestD is not None):
        bestlabel = ''.join((('Best-fit\ninitial=', '{:.2f}'.format(besti), 
                     '\nlog$_{10}$D=', '{:.2f}'.format(bestD))))                     
        ax.plot(besti, bestD, '*y', markersize=18, label=bestlabel)

    if show_max is True:
        m = get_max(wb, peak_idx, wholeblock)
        ax.plot([m, m], ax.get_ylim(), '--g', linewidth=2,
                label='maximum area\nin profile')
    
    if show_max_init is True:
        m = get_max(wb, peak_idx, wholeblock, initials=True)
        ax.plot([m, m], ax.get_ylim(), ':b', linewidth=2,
                label='maximum area\nin initial profiles')

    if maxwbs is not None:
        ax.plot([maxwbs, maxwbs], ax.get_ylim(), '-r', linewidth=1,
                label='max. area in\nall relevant profiles')

    if len(boundaries) > 0:
        ax.plot(-100., -100., '-.g', linewidth=2, label='assumed boundary')
        
    for i in boundaries:
        ax.plot([i, i], ax.get_ylim(), '-.g', linewidth=2)

    ax.legend(loc=lloc)
    
    if makefigure is True:
        fig.subplots_adjust(left=0.2)
        return fig, ax
    else:
        return ax

#%% K3 3645
peak_idx = 0
wb = my_spectra.K3wb_6days
wbs = wbsK3
fname = 'K3'
wholeblock = True

m = get_max_global(wbs, peak_idx, wholeblock=wholeblock)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx, 
                                          wholeblock=wholeblock)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4,
                 wholeblock=wholeblock)

ax.set_ylim(-30, -10)

idx_mid = 8
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(17, -15.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx, wholeblock=wholeblock)

#%% K3 3617
peak_idx = 1
wb = my_spectra.K3wb_6days
wbs = wbsK3
fname = 'K3'
m = get_max_global(wbs, peak_idx)

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD, 
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)


ax.set_ylim(-30, -10)
ax.set_xlim(0, 60)

idx_min = -7
idx_max = -1

p = np.polyfit(ilist[idx_min:idx_max], Dlist[idx_min:idx_max], 1)
Dmpolyval = np.polyval(p, m)
ax.plot([m, ax.get_xlim()[1]], np.polyval(p, [m, ax.get_xlim()[1]]), '-k', alpha=0.3)

ax.annotate('{:.1f}'.format(Dmpolyval), xy=(m, Dmpolyval), 
            xytext=(10, Dmpolyval),
            arrowprops=dict( arrowstyle='->'))
ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -16), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K3 3540
peak_idx = 2
wb = my_spectra.K3wb_6days
wbs = wbsK3
fname = 'K3'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD, 
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)


ax.set_ylim(-30, -10)
ax.set_xlim(0, 60)

idx_min = -7
idx_max = -1

p = np.polyfit(ilist[idx_min:idx_max], Dlist[idx_min:idx_max], 1)
Dmpolyval = np.polyval(p, m)
ax.plot([m, ax.get_xlim()[1]], np.polyval(p, [m, ax.get_xlim()[1]]), '-k', alpha=0.3)

ax.annotate('{:.1f}'.format(Dmpolyval), xy=(m, Dmpolyval), 
            xytext=(10, Dmpolyval),
            arrowprops=dict( arrowstyle='->'))
ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -16), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)


#%% K3 3443
peak_idx = 4
wb = my_spectra.K3wb_6days
wbs = wbsK3
fname = 'K3'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx,  besti=besti, bestD=bestD, 
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=2)

ax.set_ylim(-30, -10)

idx_min = -6
idx_max = -1

p = np.polyfit(ilist[idx_min:idx_max], Dlist[idx_min:idx_max], 1)
Dmpolyval = np.polyval(p, m)
ax.plot([m, ax.get_xlim()[1]], np.polyval(p, [m, ax.get_xlim()[1]]), '-k', alpha=0.3)

ax.annotate('{:.1f}'.format(Dmpolyval), xy=(m, Dmpolyval), 
            xytext=(24, Dmpolyval),
            arrowprops=dict( arrowstyle='->'))
ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -16), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K3 3350
peak_idx = 5
wb = my_spectra.K3wb_6days
wbs = wbsK3
fname = 'K3'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-30, -10)
#
idx_min = 9
idx_max = -1

p = np.polyfit(ilist[idx_min:idx_max], Dlist[idx_min:idx_max], 1)
Dmpolyval = np.polyval(p, m)
ax.plot([m, ax.get_xlim()[1]], np.polyval(p, [m, ax.get_xlim()[1]]), '-k', alpha=0.3)

ax.annotate('{:.1f}'.format(Dmpolyval), xy=(m, Dmpolyval), 
            xytext=(4, Dmpolyval),
            arrowprops=dict( arrowstyle='->'))
ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -16), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K4 154hr 3645
peak_idx = 0
wb = my_spectra.K4wb_154hr
wbs = wbsK4_154
fname = 'K4_154'
wholeblock = False

m = get_max_global(wbs, peak_idx, wholeblock=wholeblock)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx, 
                                          wholeblock=wholeblock)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4,
                 wholeblock=wholeblock)

ax.set_ylim(-20, -10)

idx_mid = 6
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(5, -14.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx, wholeblock=wholeblock)

#%% K4 3617
peak_idx = 1
wb = my_spectra.K4wb_154hr
wbs = wbsK4_154
fname = 'K4_154'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-20, -10)

idx_max = -1
ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -13.5), arrowprops=dict( arrowstyle='->'))

idx_mid = 13
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(22, -12), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K4 3540
peak_idx = 2
wb = my_spectra.K4wb_154hr
wbs = wbsK4_154
fname = 'K4_154'
maxwbs = get_max_global(wbs, peak_idx)
print maxwbs

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=maxwbs, boundaries=[], lloc=4)

ax.set_ylim(-14, -12)

idx_mid = 11
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(18, -12.5), arrowprops=dict( arrowstyle='->'))

idx_mid = 4
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(12, -12.6), arrowprops=dict( arrowstyle='->'))


ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -12.5), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)


#%% K4 3443
peak_idx = 4
wb = my_spectra.K4wb_154hr
wbs = wbsK4_154
fname = 'K4_154'
maxwbs = get_max_global(wbs, peak_idx)
print maxwbs

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=maxwbs, boundaries=[], lloc=4)

ax.set_ylim(-30, -10)

idx_mid = 15
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(26, -15.5), arrowprops=dict( arrowstyle='->'))

ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -15.5), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)


#%% K4 3350
peak_idx = 5
wb = my_spectra.K4wb_154hr
wbs = wbsK4_154
fname = 'K4_154'
maxwbs = get_max_global(wbs, peak_idx)
print maxwbs

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=maxwbs, boundaries=[], lloc=4)

ax.set_ylim(-20, -10)

idx_mid = 5
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(13, -14), arrowprops=dict( arrowstyle='->'))

idx_mid = 3
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(4, -14.3), arrowprops=dict( arrowstyle='->'))

ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -13.5), arrowprops=dict( arrowstyle='->'))
#
save_sensitivity_fig(fname, peak_idx)

#%% K5 3645
peak_idx = 0
wb = my_spectra.K5wb_75hr
wbs = wbsK5
fname = 'K5'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-20, -10)

ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -13.5), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K5 3617
peak_idx = 1
wb = my_spectra.K5wb_75hr
wbs = wbsK5
fname = 'K5'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-18, -12)

idx_mid = 13
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(29, -13.), arrowprops=dict( arrowstyle='->'))

ax.text(15, -13, '~ measured initials -->', rotation=90, va='top')

idx_mid = 7
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(17, -13.5), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K5 3540
peak_idx = 2
wb = my_spectra.K5wb_75hr
wbs = wbsK5
fname = 'K5'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.annotate('{:.1f}'.format(Dlist[-1]), xy=(ilist[-1], Dlist[-1]), 
            xytext=(40, -13.), arrowprops=dict( arrowstyle='->'))

ax.text(20, -12.5, '~ measured initials -->', rotation=90, va='top')

idx_mid = 9
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(15, -13.), arrowprops=dict( arrowstyle='->'))

idx_mid = 11
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(25, -13.), arrowprops=dict( arrowstyle='->'))


save_sensitivity_fig(fname, peak_idx)

#%% K5 3443
peak_idx = 4
wb = my_spectra.K5wb_75hr
wbs = wbsK5
fname = 'K5'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-20, -10)

idx_mid = 15
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(26, -12), arrowprops=dict( arrowstyle='->'))

idx_mid = 10
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(22, -11), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)


#%% K5 3350
peak_idx = 5
wb = my_spectra.K5wb_75hr
wbs = wbsK5
fname = 'K5'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-14, -12)

idx_mid = 5
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(14, -12.7), arrowprops=dict( arrowstyle='->'))

idx_mid = 3
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(7.5, -13.2), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)


#%% J 3645
peak_idx = 0
wb = my_spectra.J1wb
wbs = wbsJ
fname = 'J'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-12, -9)

idx_mid = 1
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(7, -10.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% J 3617
peak_idx = 1
wb = my_spectra.J1wb
wbs = wbsJ
fname = 'J'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-13, -9)

idx_mid = 3
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(7, -10.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% J 3540
peak_idx = 2
wb = my_spectra.J1wb
wbs = wbsJ
fname = 'J'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

idx_mid = 8
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(20, -11.), arrowprops=dict( arrowstyle='->'))

idx_mid = 5
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(12, -11.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% J 3443
peak_idx = 4
wb = my_spectra.J1wb
wbs = wbsJ

fname = 'J'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

idx_mid = 7
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(20, -11.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% J 3350
peak_idx = 5
wb = my_spectra.J1wb
wbs = wbsJ
fname = 'J'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

idx_mid = -3
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(33, -12.5), arrowprops=dict( arrowstyle='->'))

idx_mid = -10
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(27, -13.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K4 91 hr 3645
peak_idx = 0
wb = my_spectra.K4wb_91hr
wbs = wbsK4_91
fname = 'K4_91'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

idx_mid = -1
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(40, -15.), arrowprops=dict( arrowstyle='->'))

idx_mid = 7
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(17, -15.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K4 91 hr 3617
peak_idx = 1
wb = my_spectra.K4wb_91hr
wbs = wbsK4_91
fname = 'K4_91'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

idx_mid = -1
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(40, -15.), arrowprops=dict( arrowstyle='->'))

idx_mid = 13
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(22, -15.), arrowprops=dict( arrowstyle='->'))

idx_mid = 7
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(16, -15.5), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K4 91 hr 3540
peak_idx = 2
wb = my_spectra.K4wb_91hr
wbs = wbsK4_91
fname = 'K4_91'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

idx_mid = -1
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(40, -13.), arrowprops=dict( arrowstyle='->'))

idx_mid = 11
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(18, -13.), arrowprops=dict( arrowstyle='->'))

idx_mid = 4
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(11, -13), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K4 91 hr 3443
peak_idx = 4
wb = my_spectra.K4wb_91hr
wbs = wbsK4_91
fname = 'K4_91'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-30, -10)

idx_mid = -1
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(40, -15.), arrowprops=dict( arrowstyle='->'))

idx_mid = 9
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(15, -17.), arrowprops=dict( arrowstyle='->'))

idx_mid = 15
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(25, -15), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K3 bulk
peak_idx = None
wb = my_spectra.K3wb_6days
wbs = wbsK3
fname = 'K3'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

ax.set_ylim(-17, -12)
ax.set_xlim(60, 100)

idx_mid = 17
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(70, -15.3), arrowprops=dict( arrowstyle='->'))

idx_mid = 20
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(82, -13.5), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K4 91hr bulk
peak_idx = None
wb = my_spectra.K4wb_91hr
wbs = wbsK4_91
fname = 'K4_91'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

idx_mid = 11
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(50, -12.7), arrowprops=dict( arrowstyle='->'))

idx_mid = 20
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(82, -12.3), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx)

#%% K4 154hr bulk
peak_idx = None
wb = my_spectra.K4wb_154hr
wbs = wbsK4_154
fname = 'K4_154'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

#%% K5 bulk
peak_idx = None
wb = my_spectra.K5wb
wbs = wbsK5
fname = 'K5'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

#%% J bulk
peak_idx = None
wb = my_spectra.J1wb
wbs = wbsJ
fname = 'J'
m = get_max_global(wbs, peak_idx)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4)

#%% K3 bulk WB
peak_idx = None
wb = my_spectra.K3wb_6days
wbs = wbsK3
fname = 'K3'
wholeblock = True

m = get_max_global(wbs, peak_idx, wholeblock=wholeblock)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx, 
                                          wholeblock=wholeblock)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4,
                 wholeblock=wholeblock)

idx_mid = 13
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(1.3, -15.), arrowprops=dict( arrowstyle='->'))

idx_mid = 9
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(1.1, -15.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx, wholeblock=wholeblock)

#%% K4 91hr bulk WB
peak_idx = None
wb = my_spectra.K4wb_91hr
wbs = wbsK4_91
fname = 'K4_91'
wholeblock = True
m = get_max_global(wbs, peak_idx, wholeblock=wholeblock)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx, 
                                          wholeblock=wholeblock)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4,
                 wholeblock=wholeblock)

idx_mid = 13
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(1.3, -15.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx, wholeblock=wholeblock)

#%% K4 154hr bulk WB
peak_idx = None
wb = my_spectra.K4wb_154hr
wbs = wbsK4_154
fname = 'K4_154'
wholeblock = True

m = get_max_global(wbs, peak_idx, wholeblock=wholeblock)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx, 
                                          wholeblock=wholeblock)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4,
                 wholeblock=wholeblock)

idx_mid = 13
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(1.3, -15.), arrowprops=dict( arrowstyle='->'))

idx_mid = 9
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(0.9, -15.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx, wholeblock=wholeblock)

#%% J bulk WB
peak_idx = None
wb = my_spectra.J1wb
wbs = wbsJ
fname = 'J'
wholeblock = True

m = get_max_global(wbs, peak_idx, wholeblock=wholeblock)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx, 
                                          wholeblock=wholeblock)

di = solve_bestd(wb, peak_idx, initial=m, Dguess=-12, wholeblock=True)
print di

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4,
                 wholeblock=wholeblock)
                    
idx_mid = 15
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(1.4, -11.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx, wholeblock=wholeblock)

#%% K5 bulk
peak_idx = None
wb = my_spectra.K5wb_75hr
wbs = wbsK5
fname = 'K5'
wholeblock = True

m = get_max_global(wbs, peak_idx, wholeblock=wholeblock)
print m

besti, bestD, ilist, Dlist = get_initialD(fname=fname, peak_idx=peak_idx, 
                                          wholeblock=wholeblock)

fig, ax = plot_D(ilist, Dlist, wb, peak_idx, besti=besti, bestD=bestD,  
                 iused=None, show_max=True, maxwbs=m, boundaries=[], lloc=4,
                 wholeblock=wholeblock)

idx_mid = 13
ax.annotate('{:.1f}'.format(Dlist[idx_mid]), xy=(ilist[idx_mid], Dlist[idx_mid]), 
            xytext=(1.3, -15.), arrowprops=dict( arrowstyle='->'))

save_sensitivity_fig(fname, peak_idx, wholeblock=wholeblock)

#%% Figure(s) for paper...
nrow = 6
ncol = 5

fig, axarr = plt.subplots(nrow, ncol, sharey='row')
fig.set_size_inches(6.5, 9.)


xlabels = ['K3', 'K4 (91hr)', 'K4 (154hr)', 'K5', 'Jaipur']

for col in range(ncol):
    axarr[0, col].set_title(xlabels[col])

ylabels = ['3645 cm$^{-1}$', '3617 cm$^{-1}$', '3540 cm$^{-1}$', '3443 cm$^{-1}$', 
           '3350 cm$^{-1}$', 'bulk H\n(whole-block)']

peak_idx_list = [0, 1, 2, 4, 5, None]
wholeblock_list = [False]*5 + [True]
fname_list = ['K3', 'K4_91', 'K4_154', 'K5', 'J']

iplots = [my_spectra.i_K3, my_spectra.i_K91, my_spectra.i_K154,
          my_spectra.i_K5, my_spectra.i_J]

xlims = [
         [(0, 40)] * 5,
         [(0, 40)] * 5,
         [(0, 40)] * 5,
         [(0, 40)] * 5,
         [(0, 40)] * 5,
         [(0.5, 2.)] * 5,
        ]

ylims = [
         [(-25, -8)] * 5,
         [(-25, -8)] * 5,
         [(-25, -8)] * 5,
         [(-25, -8)] * 5,
         [(-25, -8)] * 5,
         [(-25, -8)] * 5,
        ]
          
xticksteps = [
              [20]*5,
              [20]*5,
              [20]*5,
              [20]*5,
              [20]*5,
              [0.5]*5,
              ]
yticksteps = [5]*6

idx = 0
for row in range(nrow):
    axarr[row, 0].set_ylabel(ylabels[row])
    for col in range(ncol):
        plt.sca(axarr[row, col])
                
        plt.xticks(np.arange(0., 150., xticksteps[row][col]))
        plt.yticks(np.arange(-50., 0., yticksteps[row]))
        plt.ylim(ylims[row][col])
        plt.xlim(xlims[row][col])

        besti, bestD, ilist, Dlist = get_initialD(fname=fname_list[col], 
                                              peak_idx=peak_idx_list[row],
                                              wholeblock=wholeblock_list[row])
        axarr[row, col].plot(ilist, Dlist, 'or', markersize=3)

        initial = iplots[col][row]        
        axarr[row, col].plot([initial, initial], axarr[row, col].get_ylim(),
                        '-k', linewidth=1.5, color='b')
                        
        if (idx < 26):
            letter = string.ascii_uppercase[idx]
        else:
            letter = string.ascii_uppercase[idx-26]*2
        rax = 0.9*(xlims[row][col][1] - xlims[row][col][0])
        ray = 0.1*(ylims[row][col][1] - ylims[row][col][0])
        letterdot = ''.join((letter, '.'))
        plt.text(xlims[row][col][0] + rax, ylims[row][col][0] + ray, 
                 letterdot, backgroundcolor='w',
                 ha='right')
        idx = idx + 1
#        axarr[row, col].set_title('{:.1f}'.format(initial), fontsize=12)
       

#labels = ['hi']
#for ax in axarr:
#    ax.set_xticks(ticks,)
#    ax.set_xticks(ticks)
#    ax.set_xticks(np.arange(min(x), max(x)+1, 1.0))
#    ax.set_xlim(50, 100)


ax = axarr[nrow/2, 0]
xm = ax.get_xlim()[0]
ax.text(xm-45, 20, 'best-fit log$_{10}$D$_c$ (m$^2$/s)', 
        rotation=90, fontsize=14)
axarr[-1, ncol/2].set_xlabel('initial area (cm$^{-2}$ for peak-specific)', fontsize=14)

#fig.autofmt_xdate()
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1,
                    wspace=0.3, hspace=0.25)


plt.savefig('../../../CpxPaper/figures/Sensitivities.png', dpi=100)

fig.show()
print 'Finished'