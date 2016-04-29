# -*- coding: utf-8 -*-+
"""
Created on Wed Jun 24 17:55:02 2015

@author: Ferriss
Before and after heating comparison of FTIR spectra on the rims of cpx samples
"""

import cpx_spectra
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import pynams.diffusion as diffusion
import string
import json
from lmfit import minimize

#%% DIFFUSIVITIES 
### Saved preferred initial values (i), diffusivities (D), and errors (e) 
### for Kunlun diopside K3, K4 at 91 hours of heating (K91) and 154 hours
### of heating (K154), K5, and Jaipur diopside J1 (J)
### In peak wavenumber order: 3645, 3617, 3540, 3443, 3355, BULK H #####
### The same numbers are used in Fig. 7
### The data and marker styles are set in cpx_spectra.py
i_K3 = cpx_spectra.i_K3
D_K3 = cpx_spectra.D_K3
D_K3 = cpx_spectra.e_K3

i_K91 = cpx_spectra.i_K91
D_K91 = cpx_spectra.D_K91
D_K91 = cpx_spectra.e_K91

i_K154 = cpx_spectra.i_K154
D_K154 = cpx_spectra.D_K154
D_K154 = cpx_spectra.e_K154

i_K5 = cpx_spectra.i_K5
D_K5 = cpx_spectra.D_K5
D_K5 = cpx_spectra.e_K5

i_J = cpx_spectra.i_J
D_J = cpx_spectra.D_J
D_J = cpx_spectra.e_J

#%% Data setup
#iwater_Kunlun = cpx_spectra.K_water_peaks
#iwater_Jaipur = cpx_spectra.J_water_peaks
#iwater_Kunlun_bulk = sum(iwater_Kunlun)
#iwater_Jaipur_bulk = sum(iwater_Jaipur)

# all whole block data associated with Kunlun diopside sample K3
wbsK3 = [cpx_spectra.K3wb_trueInit,
         cpx_spectra.K3wb_init,
         cpx_spectra.K3wb_800C_15hr,
         cpx_spectra.K3wb_6days,
         cpx_spectra.K4wb_quench,
         cpx_spectra.K4wb_1hr,
         ]

# all whole block data relevant to Kunlun diopside sample K4
# only up to heating for 91 hours
wbsK4_91 = [
           cpx_spectra.K4wb_init,
           cpx_spectra.K3wb_init,
           cpx_spectra.K3wb_800C_15hr,
           cpx_spectra.K4wb_quench,
           cpx_spectra.K4wb_1hr,
           cpx_spectra.K4wb_91hr
           ]

# all whole block data relevant to Kunlun diopside sample K4
# including heating for 154 hours
wbsK4_154 = [
       cpx_spectra.K4wb_init,
       cpx_spectra.K3wb_init,
       cpx_spectra.K3wb_800C_15hr,
       cpx_spectra.K4wb_quench,
       cpx_spectra.K4wb_1hr,
       cpx_spectra.K4wb_91hr,
       cpx_spectra.K4wb_154hr,
       ]

# all whole block data relevant to Kunlun diopside sample K5
wbsK5 = [
       cpx_spectra.K5wb_init,
       cpx_spectra.K3wb_init,
       cpx_spectra.K3wb_800C_15hr,
       cpx_spectra.K4wb_quench,
       cpx_spectra.K4wb_1hr,
       cpx_spectra.K5wb
       ]

# all whole block data relevant to Jaipur diopside J1
wbsJ = [
       cpx_spectra.J1wb_initial,
       cpx_spectra.J1wb
       ]

# list of initial whole block data - 1 for each diopside
init_list = [
            cpx_spectra.K3wb_trueInit,
            cpx_spectra.K4wb_init,
            cpx_spectra.K5wb_init,
            cpx_spectra.J1wb_initial]

wb_list = wbsK3 + wbsK4_154 + wbsK5 + wbsJ + wbsK4_91

for wb in wb_list + init_list:
    wb.get_baselines()
    wb.get_peakfit()
    wb.setupWB()
    for prof in wb.profiles:
        prof.make_wholeblock(True, True)
        # Default to all Kunlun scaled water
#        prof.wb_waters = np.array(prof.wb_areas) * iwater_Kunlun_bulk


#%% Plotting details and functions
style_bulk = {'marker' : 'o', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'k', 'markersize' : 7}

style_peak0 = {'marker' : 'x', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'r', 'markersize' : 6, 'label' : '3645 cm$^{-1}$'}

style_peak1 = {'marker' : '+', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'orange', 'markersize' : 6, 'label' : '3617 cm$^{-1}$'}

style_peak2 = {'marker' : 's', 'fillstyle' : 'full', 'linestyle' : 'none',
              'color' : 'k', 'markersize' : 6, 'markerfacecolor' : 'teal',
              'alpha' : 0.5, 'label' : '3540 cm$^{-1}$'}

style_peak3 = {'marker' : '_', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'g', 'markersize' : 6, 'mew' : 2,
              'label' : '3460 cm$^{-1}$'}

style_peak4 = {'marker' : '|', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'b', 'markersize' : 6, 'mew' : 2,
              'label' : '3443 cm$^{-1}$'}

style_peak5 = {'marker' : 'd', 'fillstyle' : 'full', 'linestyle' : 'none',
              'color' : 'k', 'markersize' : 6, 'markerfacecolor' : 'violet',
              'alpha' : 0.5, 'label' : '3355 cm$^{-1}$'}


style_init_K3 =  {'marker' : 'o', 'markeredgecolor' : 'r', 'linestyle' : 'none',
                  'markersize' : 6, 'markerfacecolor' : 'w',
                  'label' : 'initial K3', 'mew' : 1, 'alpha' : 1.}

style_init_K4 =  {'marker' : 's', 'markeredgecolor' : 'b', 'linestyle' : 'none',
                  'markersize' : 6, 'markerfacecolor' : 'w',
                  'label' : 'initial K4', 'mew' : 1, 'alpha' : 0.5}

style_init_K5 =  {'marker' : '^', 'markeredgecolor' : 'g', 'linestyle' : 'none',
                  'markersize' : 6, 'markerfacecolor' : 'w',
                  'label' : 'initial K5', 'mew' : 1, 'alpha' : 0.5}

style_init_J =  {'marker' : 'D', 'markeredgecolor' : 'k', 'linestyle' : 'none',
                 'markersize' : 6, 'markerfacecolor' : 'w',
                 'label' : 'initial J1', 'mew' : 1, 'alpha' : 0.5}

style_K4q =  {'marker' : 'x', 'color' : 'b', 'linestyle' : 'none',
              'markersize' : 6, 'label' : '480 $\degree$C, 0.6hr',
              'mew' : 1, 'alpha' : 0.5}

style_K3q =  {'marker' : '+', 'color' : 'r', 'linestyle' : 'none',
              'markersize' : 6, 'label' : '696 $\degree$C, 2hr', 'mew' : 1,
              'alpha' : 1.}

style_K3_15h =  {'marker' : '3', 'color' : 'r', 'linestyle' : 'none',
                 'markersize' : 6, 'label' : '796 $\degree$C, 15hr', 'mew' : 1,
                 'alpha' : 0.5}

style_K3_6d =  {'marker' : 'o', 'color' : 'r', 'linestyle' : 'none',
                'markersize' : 6, 'label' : '812 $\degree$C, 140hr',
                'alpha' : 0.5, 'mew' : 1}

style_K4_1h =  {'marker' : '|', 'color' : 'c', 'linestyle' : 'none',
                'markersize' : 6, 'label' : '904 $\degree$C, 0.7hr',
                'mew' : 1, 'alpha' : 0.75}

style_K4_91 =  {'marker' : 'd', 'color' : 'c', 'linestyle' : 'none',
                'markersize' : 6, 'label' : '904 $\degree$C, 91hr',
                'alpha' : 0.5}

style_K4_154 =  {'marker' : 's', 'color' : 'b', 'linestyle' : 'none',
                'markersize' : 6, 'label' : '904 $\degree$C, 154hr',
                'alpha' : 0.5}

style_K5 =  {'marker' : '^', 'color' : 'g', 'linestyle' : 'none',
             'markersize' : 6, 'label' : '1000 $\degree$C, 75hr',
             'alpha' : 0.5}

style_J =  {'marker' : 'D', 'color' : 'k', 'linestyle' : 'none',
            'markersize' : 6, 'label' : 'Jaipur\n904 $\degree$C, 0.6 hr',
            'alpha' : 0.5}

style_iline = {'linestyle' : '--', 'color' : 'k',
               'label' : '"initial" line\nfor D model'}

styledict={cpx_spectra.K3wb_trueInit : style_init_K3,
           cpx_spectra.K3wb_init : style_K3q,
           cpx_spectra.K3wb_800C_15hr : style_K3_15h,
           cpx_spectra.K3wb_6days : style_K3_6d,
           cpx_spectra.K4wb_init : style_init_K4,
           cpx_spectra.K4wb_quench : style_K4q,
           cpx_spectra.K4wb_1hr : style_K4_1h,
           cpx_spectra.K4wb_91hr : style_K4_91,
           cpx_spectra.K4wb_154hr : style_K4_154,
           cpx_spectra.K5wb_init : style_init_K5,
           cpx_spectra.K5wb : style_K5,
           cpx_spectra.J1wb_initial : style_init_J,
           cpx_spectra.J1wb : style_J}

# For least squares fitting to obtain diffusivities
func2min = diffusion.diffusion3Dwb_params

def get_best_initialD(fname, peak_idx):
    """Takes name = 'K4_91', 'K3', 'K5', 'K4_154', or 'J' and peak index
    Returns best-fit initial and isotropic diffusivity generated by
    Sentivity.py
    """
    fname_list = ['K4_91', 'K3', 'K5', 'K4_154', 'J']
    if fname not in fname_list:
        print 'fname must be one of the following:'
        print fname_list
        return False, False
    fi = '../../../CpxPaper/figures/sensitivity_'
    peak_label_list = ['3645', '3617', '3540', '3460', '3443', '3355']
    filename = ''.join((fi, fname, '_', peak_label_list[peak_idx], '.txt'))
    f = open(filename, 'r')
    x = json.load(f)
    besti = x[0][0]
    bestD = x[1][0]
    return besti, bestD

def figOutline(top=30., ylab='peak area', ncol=3, nrow=4, figsize=(6.5, 8),
               yTextPos=0.5, tit=None):
    """Make and return figure and axes[12]"""
    fig = plt.figure()
    fig.set_size_inches(figsize)

    gs = gridspec.GridSpec(nrow, ncol)

    axes = []
    for row in range(nrow):
        for col in range(ncol):
            axes.append(plt.subplot(gs[row, col]))

    axes[0].set_title('Profile || a*')
    axes[1].set_title('Profile || b')
    axes[2].set_title('Profile || c')

    if tit is None:
        if ylab == 'peak area':
            fig.text(0.05, yTextPos, 'Peak area (cm$^{-2}$)',
                     ha='center', va='center', rotation='vertical', fontsize=14)
#            fig.text(0.06, yTextPos, 'Bell et al. 1995 calibration: 7.09 +/- 0.32 cm$^{-2}$/ppm H$_2$O',
#                     ha='center', va='center', rotation='vertical', fontsize=12)
        elif ylab == 'K5':
            fig.text(0.04, yTextPos, 'Peak area in Kunlun diopside after 75 hr at 1000 $\degree$C (cm$^{-2}$)',
                     ha='center', va='center', rotation='vertical', fontsize=14)
        else:
            fig.text(0.04, yTextPos, 'Estimated water (ppm H$_2$O) using whole-block data\nscaled up using initial water estimates from polarized data',
                     ha='center', va='center', rotation='vertical', fontsize=14)
    else:
         fig.text(0.04, yTextPos, tit, ha='center', va='center',
                  rotation='vertical', fontsize=14)

    ntwins = len(np.arange(2, len(axes), 3))
    for idx in np.arange(2, len(axes), 3):
        axes.append(axes[idx].twinx())
#        axes[idx].set_ylim(0, 100)

    for ax in axes:
        ax.set_xlim(0., 1.)
        ax.set_ylim(0, top)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)        

    for idx in np.arange(0, 3*nrow, 3):
        plt.setp(axes[idx].get_yticklabels(), visible=True)

    nax = len(axes)
    for idx in [nax-ntwins-3, nax-ntwins-2, nax-ntwins-1]:
        plt.setp(axes[idx].get_xticklabels(), visible=True, rotation=45)

    axes[nax-ntwins-2].set_xlabel('normalized distances (X / Length)', fontsize=14)

    fig.subplots_adjust(bottom=0.1, top=0.95, right=0.88)
    return fig, axes

def plotpeaks(wb_list, wbs_list, ilist, dlist, elist, dlist_slow=None,
              slowb=[False]*6, legidx=13, sidelabels=True, ncol=5,
              wholeblock=False, peak_idx_list = [0, 1, 2, 4, 5],
              show_legend=True, dlabel='Isotropic\nDiffusion curve',
              xtickgrid=[0.2]*6, ytickgrid=[5.]*6):
    """Takes wb = list of main whole block being plotted,
    wbs = list of lists of whole block data being plotted in each panel,
    ilist = list of initials, dlist = list of diffusivities,
    elist = list of errors"""
    #### legend ####
    if show_legend is True:
        styles = []
        for wbs in wbs_list:
            for wbtoplot in wbs:
                styleToAdd = styledict[wbtoplot]
                if styleToAdd not in styles:
                    styles.append(styleToAdd)

        for sty in styles:
            axes[legidx].plot(-100, -100, **sty)
        axes[legidx].plot(-100, -100, '-k', label='"initial"', linewidth=1)
        axes[legidx].plot(-100, -100, '-g', label=dlabel,
                          linewidth=3)
        axes[legidx].plot(-100, -100, '--g', label='+/- error log10 D')
        axes[legidx].legend(fancybox='on', loc=8) # 9 for top

    axesranges = [range(3), range(3,6), range(6,9), range(9,12), range(12,15)]

    ### y axis limits and "initial"
    for r in range(ncol):
        for ax_idx in axesranges[r]:
            ax = axes[ax_idx]
            ax.plot(ax.get_xlim(), [ilist[r], ilist[r]], '-k', linewidth=1)
            ax.set_ylim(0, tops[r])

    xd = [] # x data to plot
    yd = [] # y data to plot
    yf = []
    ys = []

    ax_idx = 0
    for idx in range(ncol): # loop through each row of data

        # get basic data
        peak_idx = peak_idx_list[idx]
        initial = ilist[idx]
        wb = wb_list[idx]
        wbs = wbs_list[idx]

        # diffusion curves
        L3 = []
        D3 = []
        Dfast = []
        Dslow = []       
        for k in range(3):
            L3.append(wb.profiles[k].len_microns)
            D3.append(dlist[idx])
        if slowb[idx] is True:
            D3[1] = dlist_slow[idx]
        for k in range(3):
            Dfast.append(D3[k] + elist[idx])
            Dslow.append(D3[k] - elist[idx])
        params = diffusion.params_setup3D(L3, D3, wb.time_seconds, initial)
        xdiff, ydiff = diffusion.diffusion3Dwb_params(params,
                                  raypaths=wb.raypaths, show_plot=False)
        params = diffusion.params_setup3D(L3, Dfast, wb.time_seconds, initial)
        xdiff, yfast = diffusion.diffusion3Dwb_params(params, 
                                                      raypaths=wb.raypaths, 
                                                      show_plot=False)
        params = diffusion.params_setup3D(L3, Dslow, wb.time_seconds, initial)
        xdiff, yslow = diffusion.diffusion3Dwb_params(params, 
                                                      need_to_center_x_data=False,
                                                      raypaths=wb.raypaths, 
                                                      show_plot=False)
        for k in range(3): 
            m = max(xdiff[k]) 
            xdiff[k] = xdiff[k] / m 
        xd.append(xdiff)
        yd.append(ydiff)
        yf.append(yfast)
        ys.append(yslow)

        ### side labels ####
        ylabelloc = 'left'
        formatter = '{:.2f}'
        if slowb[idx] is True:
            subscript = '$_c$'
        else:
            subscript = ''
    
        if sidelabels is True:
            axes[15].set_ylabel(''.join(('3645\ncm$^{-1}$\n\nlog$_{10}$D', subscript,
                                '\n', formatter.format(peak_D[0]), '\n+/-', str(er[0]))),
                                rotation=0, ha=ylabelloc, va='center')
            axes[16].set_ylabel(''.join(('3617\ncm$^{-1}$\n\nlog$_{10}$D', subscript, '\n',
                                formatter.format(peak_D[1]), '\n+/-', str(er[1]))),
                                rotation=0, ha=ylabelloc, va='center')
            axes[17].set_ylabel(''.join(('3540\ncm$^{-1}$\n\nlog$_{10}$D', subscript, '\n',
                                formatter.format(peak_D[2]), '\n+/-', str(er[2]))),
                                rotation=0, ha=ylabelloc, va='center')
            axes[18].set_ylabel(''.join(('3443\ncm$^{-1}$\n\nlog$_{10}$D', subscript, '\n',
                                formatter.format(peak_D[4]), '\n+/-', str(er[4]))),
                                rotation=0, ha=ylabelloc, va='center')
            axes[19].set_ylabel(''.join(('3355\ncm$^{-1}$\n\nlog$_{10}$D', subscript, '\n',
                                formatter.format(peak_D[5]), '\n+/-', str(er[5]))),
                                rotation=0, ha=ylabelloc, va='center')

        for k in range(3):
            ax = axes[ax_idx]

            # Set tick locations
            xmajorLocator = MultipleLocator(xtickgrid[idx])
            ymajorLocator = MultipleLocator(ytickgrid[idx])
            ax.xaxis.set_major_locator(xmajorLocator)
            ax.yaxis.set_major_locator(ymajorLocator)

            if show_legend is True:
                if ax_idx < legidx:
                    s = ''.join((string.ascii_uppercase[ax_idx],'.'))
                elif ax_idx > legidx:
                    s = ''.join((string.ascii_uppercase[ax_idx-1],'.'))
            else:
                s = ''.join((string.ascii_uppercase[ax_idx],'.'))
            ax.text(0.1, ax.get_ylim()[1] - 0.15*ax.get_ylim()[1], s,
                    backgroundcolor='w')

            # Plot the diffusion curves
            ax.plot(xd[idx][k], yd[idx][k], '-g', linewidth=3)
            ax.plot(xd[idx][k], yf[idx][k], '--g')
            ax.plot(xd[idx][k], ys[idx][k], '--g')

            # Plot all wb data in wbs list of wholeblocks
            for wb_idx in range(len(wbs)):
                if wbs[wb_idx].raypaths[k] == wb.raypaths[k]:
                    xi3, yi3 = wbs[wb_idx].xy_picker(peak_idx=peak_idx,
                                                     wholeblock=wholeblock,
                                                     centered=False)
                    xi = xi3[k]
                    yi = yi3[k]

                    # Normalize x
                    L = wbs[wb_idx].profiles[k].len_microns
                    xi = np.array(xi) / L

                    ax.plot(xi, yi, **styledict[wbs[wb_idx]])

            # label ray paths
            if wb.profiles[k].raypath == 'a':
                R = ''.join(('R || ', wb.profiles[k].raypath,'*'))
            else:
                R = ' '.join(('R ||', wb.profiles[k].raypath))
            axes[ax_idx].text(0.65, ax.get_ylim()[1] - 0.15*ax.get_ylim()[1], R)

            ax_idx = ax_idx + 1

def get_besti(wb, peak_idx, Dguess=-12):
    """Returns best-fit initial and isotropic D """
    x, y = wb.xy_picker(peak_idx=peak_idx, wholeblock=False, heights_instead=False)
    L3 = []
    D3 = []
    for k in range(3):
        L3.append(wb.profiles[k].len_microns)
        D3.append(Dguess)
    # best fit initial and corresponding D
    params = diffusion.params_setup3D(L3, D3, wb.time_seconds, initial=17., isotropic=True, vinit=True)
    minimize(func2min, params, args=(x, y), kws={'raypaths' : wb.raypaths, 'show_plot' : False})
    besti = params['initial_unit_value'].value
    bestD = params['log10Dx'].value
    return besti, bestD

def get_bestd(wb, peak_idx, initial, Dguess=-12):
    """Returns best-fit initial and isotropic D """
    x, y = wb.xy_picker(peak_idx=peak_idx, wholeblock=False, heights_instead=False)
    L3 = []
    D3 = []
    for k in range(3):
        L3.append(wb.profiles[k].len_microns)
        D3.append(Dguess)
    # best fit initial and corresponding D
    params = diffusion.params_setup3D(L3, D3, wb.time_seconds, initial=initial,
                                      isotropic=True, vinit=False)
    minimize(func2min, params, args=(x, y), kws={'raypaths' : wb.raypaths, 'show_plot' : False})
    bestD = params['log10Dx'].value
    return bestD


#%% Figure 5 - peak specific profiles
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5))
ylabelloc = 'left'
axes[15].set_ylabel('Kunlun\ndiopside\n"initial"\n696 $\degree$C\n2 hr', rotation=0, ha=ylabelloc, va='center')
axes[16].set_ylabel('Kunlun\ndiopside\n812 $\degree$C\n6 days', rotation=0, ha=ylabelloc, va='center')
axes[17].set_ylabel('Kunlun\ndiopside\n904 $\degree$C\n154 hr', rotation=0, ha=ylabelloc, va='center')
axes[18].set_ylabel('Kunlun\ndiopside\n1000 $\degree$C\n75 hr', rotation=0, ha=ylabelloc, va='center')
axes[19].set_ylabel('Jaipur\ndiopside\n904 $\degree$C\n0.6 hr', rotation=0, ha=ylabelloc, va='center')

axes_list = range(0, 15)
ax_idx = 0

wbs5 = [cpx_spectra.K3wb_init,
        cpx_spectra.K3wb_6days,
        cpx_spectra.K4wb_154hr,        
        cpx_spectra.K5wb,
        cpx_spectra.J1wb,]
        
for wb in wbs5:
    print wb.name
    for k in range(3):
        L = wb.profiles[k].len_microns
        x = np.array(wb.profiles[k].positions_microns) / L

        if wb.profiles[k].raypath == 'a':
            R = ''.join(('R || ', wb.profiles[k].raypath,'*'))
        else:
            R = ' '.join(('R ||', wb.profiles[k].raypath))
      
        if ax_idx > 0:
            s = ''.join((string.ascii_uppercase[ax_idx-1],'.'))
            axes[axes_list[ax_idx]].text(0.1, 25, s)
            axes[axes_list[ax_idx]].text(0.65, 25, R)
        
        idx = 0
        for spectrum in wb.profiles[k].spectra_list:
            print ''.join((spectrum.fname, ' (x=', '{:.2f}'.format(x[idx]), ')'))
            idx = idx + 1
        print ' '
            
        y0 = wb.profiles[k].peak_areas[0, :]
        axes[axes_list[ax_idx]].plot(x, y0, **style_peak0 )

        y1 = wb.profiles[k].peak_areas[1, :]
        axes[axes_list[ax_idx]].plot(x, y1, **style_peak1 )

        y2 = wb.profiles[k].peak_areas[2, :]
        axes[axes_list[ax_idx]].plot(x, y2, **style_peak2 )

        y3 = wb.profiles[k].peak_areas[3, :]
        axes[axes_list[ax_idx]].plot(x, y3, **style_peak3 )

        y4 = wb.profiles[k].peak_areas[4, :]
        axes[axes_list[ax_idx]].plot(x, y4, **style_peak4 )

        y5 = wb.profiles[k].peak_areas[5, :]
        axes[axes_list[ax_idx]].plot(x, y5, **style_peak5 )

        ax_idx = ax_idx + 1

leg = axes[0].legend(loc=2, ncol=1, fancybox='on', fontsize=8)
plt.savefig('Fig5_CpxProfiles.eps', format='eps', dpi=1000)
fig.savefig('Fig5_CpxProfiles.tif', format='tif', dpi=300)

#%% K3 Diffusivity modeling; Supplemental Figure 1
wb = cpx_spectra.K3wb_6days
fname = 'K3'

iareas = i_K3
peak_D = D_K3
er = e_K3
tops = [17, 35, 30, 40, 15]

style_K3_6d['alpha'] = 1.
style_K3_6d['mew'] = 1.5

wbs = [
       cpx_spectra.K3wb_init,
       cpx_spectra.K3wb_800C_15hr,
       cpx_spectra.K3wb_6days,
       cpx_spectra.K4wb_quench,
       cpx_spectra.K4wb_1hr,
       ]

tit = 'Peak area in Kunlun diopside after 6 days at 816 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)

plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er, sidelabels=False)

peak_idx_list = ['3645', '3617', '3540', '3443', '3355']
facecolors = ['wheat', 'thistle']
alphas = [0.4, 0.4]
xtext = 1.2

for k in range(5):
    ytext_top = axes[k+15].get_ylim()[1]
    ytext_top = ytext_top - 0.03*ytext_top
    ra = axes[k+15].get_ylim()[1] - axes[k+15].get_ylim()[0]
    ytext_gap = ra*0.475
    ytext_bot = ytext_top - ytext_gap 
    
    sidelabel_top = ''.join((peak_idx_list[k], '\ncm$^{-1}$'))
    sidelabel_bot = ''.join(('log$_{10}$D\n', 
                             '{:.2f}'.format(peak_D[k]), 
                             '\n+/-', '{:.1f}'.format(er[0])
                             ))
    
    axes[k+15].text(xtext, ytext_top, 
            sidelabel_top, rotation=0, va='top', ha='center',
            bbox=dict(boxstyle='round', facecolor=facecolors[0],
                      alpha=alphas[0]))
            
    axes[k+15].text(xtext, ytext_bot,
            sidelabel_bot, rotation=0, va='top', ha='center',
            bbox=dict(boxstyle='square', facecolor=facecolors[1],
                      alpha=alphas[1]))
    
fig.show()
print 'Finished'
plt.savefig('Supplement_Fig1_K3.eps', format='eps', dpi=1000)
plt.savefig('Supplement_Fig1_K3.tif', format='tif', dpi=300)

#%% K4 91 hr Diffusivity modeling; Supplemental Fig. 2
style_K4_91['alpha'] = 1.
style_K4_91['mew'] = 1.5

iareas = i_K91
peak_D = D_K91

er = e_K91
tops = [20, 35, 30, 35., 15]

wb = cpx_spectra.K4wb_91hr
wbs = [
       cpx_spectra.K4wb_init,
       cpx_spectra.K3wb_init,
       cpx_spectra.K3wb_800C_15hr,
       cpx_spectra.K4wb_quench,
       cpx_spectra.K4wb_1hr,
       cpx_spectra.K4wb_91hr
       ]

tit = 'Peak area in Kunlun diopside after 91 hr at 904 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)
plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er)
plt.savefig('Supplement_Fig2_K4_91hr.eps', format='eps', dpi=1000)
plt.savefig('Supplement_Fig2_K4_91hr.tif', format='tif', dpi=300)

#%% K4 154 hr Diffusivity modeling; Supplemental Fig. 3
tit = 'Peak area in Kunlun diopside after 154 hr at 904 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)

style_K4_154['alpha'] = 1.
style_K4_154['mew'] = 1.5

wbs = [
       cpx_spectra.K4wb_init,
       cpx_spectra.K3wb_init,
       cpx_spectra.K3wb_800C_15hr,
       cpx_spectra.K4wb_quench,
       cpx_spectra.K4wb_1hr,
       cpx_spectra.K4wb_91hr,
       cpx_spectra.K4wb_154hr,
       ]

iareas = i_K154
peak_D = D_K154

er = e_K154
tops = [21, 36, 31, 35, 16]
wb = cpx_spectra.K4wb_154hr

plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er)
plt.savefig('Supplement_Fig3_K4_154hr.eps', format='eps', dpi=1000)
plt.savefig('Supplement_Fig3_K4_154hr.tif', format='tif', dpi=300)

#%% K5 Diffusivity modeling; Supplemental Fig. 4
tit = 'Peak area in Kunlun diopside after 75 hr at 1000 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)

wb = cpx_spectra.K5wb

style_K5['alpha'] = 1.
style_K5['mew'] = 1.5

wbs = [
       cpx_spectra.K5wb_init,
       cpx_spectra.K3wb_init,
       cpx_spectra.K3wb_800C_15hr,
       cpx_spectra.K4wb_quench,
       cpx_spectra.K4wb_1hr,
       cpx_spectra.K5wb
       ]

iareas = i_K5
peak_D = D_K5
er = e_K5
tops = [23, 35, 30, 35, 15]

plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er)
plt.savefig('Supplement_Fig4_K5.eps', format='eps', dpi=1000)
plt.savefig('Supplement_Fig4_K5.tif', format='tif', dpi=300)

#%% J diffusivity modeling; Supplemental Figure 5
wb = cpx_spectra.J1wb

style_J['alpha'] = 0.5
style_J['mew'] = 1.5

wbs = [
       cpx_spectra.J1wb_initial,
       cpx_spectra.J1wb
       ]

iareas = i_J
peak_D = D_J
er = e_J
tops = [4, 8, 23, 20, 50]
peak_D_slow = np.array(peak_D) - 1.
ytickgrid = [1, 2, 5, 5, 10]

tit = 'Peak area in Jaipur diopside after 0.6 hr at 904 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)
plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er, slowb=[True]*5, dlist_slow=peak_D_slow, legidx=4,
          dlabel='Diffusion curve', ytickgrid=ytickgrid)

plt.savefig('Supplement_Fig5_J1.eps', format='eps', dpi=1000)
plt.savefig('Supplement_Fig5_J1.tif', format='tif', dpi=300)

#%% Bulk WB areas with diffusivity estimates; Fig. 6
wb_list = [cpx_spectra.K3wb_6days,
           cpx_spectra.K4wb_91hr,
           cpx_spectra.K4wb_154hr,
           cpx_spectra.J1wb,
           cpx_spectra.K5wb]

typelabels = ['Kunlun\ndiopside\n812 $\degree$C\n6 days',
              'Kunlun\ndiopside\n904 $\degree$C\n91 hr',
              'Kunlun\ndiopside\n904 $\degree$C\n154 hr',
              'Jaipur\ndiopside\n904 $\degree$C\n0.6 hr',
              'Kunlun\ndiopside\n1000 $\degree$C\n75 hr']

my_ilist = [i_K3, i_K91, i_K154,
            i_J, i_K5]
my_dlist = [D_K3, D_K91, D_K154,
            D_J, D_K5]
my_elist = [e_K3, e_K91, e_K154,
            e_J, e_K5]

tops = [ 1.6, 1.6, 1.6, 1.8, 1.6,]
slowb_list = [False, False, False, True, False, False]
wbs_list = [[wb_list[0]], [wb_list[1]], [wb_list[2]], [wb_list[3]], [wb_list[4]]]

### Moving Jaipur to the end for all of them
for li in [wb_list, typelabels, my_ilist, my_dlist, my_elist,
           tops, slowb_list, wbs_list]:
    li.insert(4, li.pop(3))

iareas = np.ones(5)
Dxz = np.ones(5)
er = np.ones(5)
for idx in range(5):
    iareas[idx] = my_ilist[idx][-1]
    Dxz[idx] = my_dlist[idx][-1]
    er[idx] = my_elist[idx][-1]
peak_D_slow = np.array(Dxz) - 1.

tit = 'Bulk H (Total area / Initial total area)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5),
                       top=1.6, tit=tit)
plotpeaks(wb_list=wb_list, wbs_list=wbs_list, ilist=iareas, dlist=Dxz,
          elist=er, slowb=slowb_list, dlist_slow=peak_D_slow,
          show_legend=False, sidelabels=False,
          wholeblock=True, peak_idx_list=[None]*5, ytickgrid=[0.5]*5)

peak_idx_list = ['3645', '3617', '3540', '3443', '3355']
facecolors = ['wheat', 'thistle']
alphas = [0.4, 0.4]
xtext = 1.25
textsize = 9

for k in range(5):
    ytext_top = axes[k+15].get_ylim()[1] / 2.
    ytext_bot = axes[k+15].get_ylim()[0]
    sidelabel_top = typelabels[k]
    
    axes[k+15].text(xtext, ytext_top, 
            sidelabel_top, rotation=0, va='center', ha='center', fontsize=textsize,
            bbox=dict(boxstyle='round', facecolor=facecolors[0],
                      alpha=alphas[0]))

idx = 0
for k in [2, 5, 8, 11, 14]:
    sidelabel_bot = ''.join(('log$_{10}$D$_c$\n', 
                             '{:.2f}'.format(my_dlist[idx][-1]), 
                             '+/-', '{:.2f}'.format(my_elist[idx][-1])
                             ))
    if k == 14:
        yloc = 0.
    else:
        yloc = 0.15

    axes[k].text(0.5, yloc, sidelabel_bot, va='bottom', ha='center',
                fontsize=textsize)
    idx = idx + 1

plt.savefig('Fig6_DiffusivityFitting.eps', format='eps', dpi=1000)
plt.savefig('Fig6_DiffusivityFitting.tif', format='tif', dpi=300)