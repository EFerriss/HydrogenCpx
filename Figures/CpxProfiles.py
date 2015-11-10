# -*- coding: utf-8 -*-+
"""
Created on Wed Jun 24 17:55:02 2015

@author: Ferriss
Before and after heating comparison of FTIR spectra on the rims of cpx samples
"""

import my_spectra
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import pynams.diffusion as diffusion
import string
import json
from lmfit import minimize

func2min = diffusion.diffusion3Dwb_params

#%% Data setup
iwater_Kunlun = my_spectra.K_water_peaks
iwater_Jaipur = my_spectra.J_water_peaks
iwater_Kunlun_bulk = sum(iwater_Kunlun)
iwater_Jaipur_bulk = sum(iwater_Jaipur)

wbsK3 = [my_spectra.K3wb_trueInit,
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

init_list = [
            my_spectra.K3wb_trueInit,
            my_spectra.K4wb_init,
            my_spectra.K5wb_init,
            my_spectra.J1wb_initial]

wb_list = wbsK3 + wbsK4_154 + wbsK5 + wbsJ + wbsK4_91

for wb in wb_list + init_list:
    wb.get_baselines()
    wb.get_peakfit()
    wb.setupWB()
    for prof in wb.profiles:
        prof.make_wholeblock(True, True)
        # Default to all Kunlun scaled water
        prof.wb_waters = np.array(prof.wb_areas) * iwater_Kunlun_bulk


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
              'alpha' : 0.5, 'label' : '3350 cm$^{-1}$'}


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


styledict={my_spectra.K3wb_trueInit : style_init_K3,
           my_spectra.K3wb_init : style_K3q,
           my_spectra.K3wb_800C_15hr : style_K3_15h,
           my_spectra.K3wb_6days : style_K3_6d,
           my_spectra.K4wb_init : style_init_K4,
           my_spectra.K4wb_quench : style_K4q,
           my_spectra.K4wb_1hr : style_K4_1h,
           my_spectra.K4wb_91hr : style_K4_91,
           my_spectra.K4wb_154hr : style_K4_154,
           my_spectra.K5wb_init : style_init_K5,
           my_spectra.K5wb_75hr : style_K5,
           my_spectra.J1wb_initial : style_init_J,
           my_spectra.J1wb : style_J}

# Area figure setup and labeling
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
    peak_label_list = ['3645', '3617', '3540', '3460', '3443', '3350']
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

    xd = []
    yd = []
    yf = []
    ys = []

    ax_idx = 0
    for idx in range(ncol):
        peak_idx = peak_idx_list[idx]
        initial = ilist[idx]
        wb = wb_list[idx]
        wbs = wbs_list[idx]

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
#            axes[15].set_ylabel(r"This is \textbf{line 1}",
#                                rotation=0, ha=ylabelloc, va='center')
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
            axes[19].set_ylabel(''.join(('3350\ncm$^{-1}$\n\nlog$_{10}$D', subscript, '\n',
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

            ax.plot(xd[idx][k], yd[idx][k], '-g', linewidth=3)
            ax.plot(xd[idx][k], yf[idx][k], '--g')
            ax.plot(xd[idx][k], ys[idx][k], '--g')

            # Plot all wb data in wbs list of wholeblocks\
            for wb_idx in range(len(wbs)):
                if wbs[wb_idx].raypaths[k] == wb.raypaths[k]:
                    xi3, yi3 = wbs[wb_idx].xy_picker(peak_idx=peak_idx,
                                            wholeblock=wholeblock)
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

#            # set wb diffusivities
#            wb.profiles[k].D_peakarea_wb_error[peak_idx] = elist[idx]
#            wb.profiles[k].D_peakarea_wb[peak_idx] = dlist[idx]

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


#%% Main profiles figure
fig, axes = figOutline()
fig.set_figsize_inches = (6.5, 7.5)
ylabelloc = 'left'
axes[12].set_ylabel('Kunlun\ndiopside\n812 $\degree$C\n6 days', rotation=0, ha=ylabelloc, va='center')
axes[13].set_ylabel('Kunlun\ndiopside\n904 $\degree$C\n154 hr', rotation=0, ha=ylabelloc, va='center')
axes[14].set_ylabel('Jaipur\ndiopside\n904 $\degree$C\n0.6 hr', rotation=0, ha=ylabelloc, va='center')
axes[15].set_ylabel('Kunlun\ndiopside\n1000 $\degree$C\n75 hr', rotation=0, ha=ylabelloc, va='center')

axes_list = range(0, 12)
ax_idx = 0
for wb in wb_list[-4:]:
#    print wb.name
    for k in range(3):
        L = wb.profiles[k].len_microns
        x = np.array(wb.profiles[k].positions_microns) / L

        if wb.profiles[k].raypath == 'a':
            R = ''.join(('R || ', wb.profiles[k].raypath,'*'))
        else:
            R = ' '.join(('R ||', wb.profiles[k].raypath))
        axes[axes_list[ax_idx]].text(0.65, 25, R)

        if ax_idx == 0:
            s = ''.join((string.ascii_uppercase[ax_idx],'.'))
        else:
            s = ''.join((string.ascii_uppercase[ax_idx-1],'.'))
        axes[axes_list[ax_idx]].text(0.1, 25, s)

    #    y = wb.profiles[k].areas_list
    #    axes[axes_list[k]].plot(x, y, **style_bulk)

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

axes[1].legend()
fig.savefig('../../../CpxPaper/figures/CpxProfiles.png', dpi=100)

#%% K3 profiles figure
fig, axes = figOutline()
ylabelloc = 'left'
axes[12].set_ylabel('Kunlun\ndiopside\n696 $\degree$C\n2 hr', rotation=0, ha=ylabelloc, va='center')
axes[13].set_ylabel('Kunlun\ndiopside\n696 $\degree$C\n19 hr', rotation=0, ha=ylabelloc, va='center')
axes[14].set_ylabel('Kunlun\ndiopside\n696 $\degree$C\n35 hr', rotation=0, ha=ylabelloc, va='center')
axes[15].set_ylabel('Kunlun\ndiopside\n796 $\degree$C\n15.7 hr', rotation=0, ha=ylabelloc, va='center')

axes_list = range(0, 12)
ax_idx = 0
for wb in wb_list[0:4]:
    print wb.name
    for k in range(3):
        L = wb.profiles[k].len_microns
        x = np.array(wb.profiles[k].positions_microns) / L

    #    y = wb.profiles[k].areas_list
    #    axes[axes_list[k]].plot(x, y, **style_bulk)

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


#%% K4 profiles figure
fig, axes = figOutline()
ylabelloc = 'left'
axes[12].set_ylabel('', rotation=0, ha=ylabelloc, va='center')
axes[13].set_ylabel('Kunlun\ndiopside K4\n480 $\degree$C\n0.6 hr', rotation=0, ha=ylabelloc, va='center')
axes[14].set_ylabel('Kunlun\ndiopside K4\n904 $\degree$C\n0.7 hr', rotation=0, ha=ylabelloc, va='center')
axes[15].set_ylabel('Kunlun\ndiopside K4\n904 $\degree$C\n91 hr', rotation=0, ha=ylabelloc, va='center')

axes_list = range(0, 12)
ax_idx = 0
for wb in [wb_list[1]] + wb_list[5:8]:
    print wb.name
    for k in range(3):
        L = wb.profiles[k].len_microns
        x = np.array(wb.profiles[k].positions_microns) / L

    #    y = wb.profiles[k].areas_list
    #    axes[axes_list[k]].plot(x, y, **style_bulk)

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

axes[1].legend()
fig.savefig('../../../CpxPaper/figures/.png', dpi=100)

#%% water profiles figure
fig, axes = figOutline(ylab='water', top=20)
ylabelloc = 'left'
axes[12].set_ylabel('Kunlun\ndiopside\n812 $\degree$C\n6 days', rotation=0, ha=ylabelloc, va='center')
axes[13].set_ylabel('Kunlun\ndiopside\n904 $\degree$C\n154 hr', rotation=0, ha=ylabelloc, va='center')
axes[14].set_ylabel('Jaipur\ndiopside\n904 $\degree$C\n0.6 hr', rotation=0, ha=ylabelloc, va='center')
axes[15].set_ylabel('Kunlun\ndiopside\n1000 $\degree$C\n75 hr', rotation=0, ha=ylabelloc, va='center')

axes_list = range(0, 12)
ax_idx = 0
for wb in wb_list[-4:]:
#    print wb.name
    for k in range(3):
        L = wb.profiles[k].len_microns
        x = np.array(wb.profiles[k].positions_microns) / L

        if wb.profiles[k].raypath == 'a':
            R = ''.join(('R || ', wb.profiles[k].raypath,'*'))
        else:
            R = ' '.join(('R ||', wb.profiles[k].raypath))
        axes[axes_list[ax_idx]].text(0.65, 17, R)

        if ax_idx == 0:
            s = ''.join((string.ascii_uppercase[ax_idx],'.'))
        else:
            s = ''.join((string.ascii_uppercase[ax_idx-1],'.'))
        axes[axes_list[ax_idx]].text(0.1, 17, s)

        if 'Jaipur' in wb.name:
            iwater = iwater_Jaipur
#            print 'Jaipur! 3350', iwater[5]
        else:
            iwater = iwater_Kunlun
#            print 'Kunlun! 3350', iwater[5]
    #    y = wb.profiles[k].areas_list
    #    axes[axes_list[k]].plot(x, y, **style_bulk)

        y0 = wb.profiles[k].peak_wb_areas[0] * iwater[0]
        axes[axes_list[ax_idx]].plot(x, y0, **style_peak0 )

        y1 = wb.profiles[k].peak_wb_areas[1] * iwater[1]
        axes[axes_list[ax_idx]].plot(x, y1, **style_peak1 )

        y2 = wb.profiles[k].peak_wb_areas[2] * iwater[2]
        axes[axes_list[ax_idx]].plot(x, y2, **style_peak2 )

        y3 = wb.profiles[k].peak_wb_areas[3] * iwater[3]
        axes[axes_list[ax_idx]].plot(x, y3, **style_peak3 )

        y4 = wb.profiles[k].peak_wb_areas[4] * iwater[4]
        axes[axes_list[ax_idx]].plot(x, y4, **style_peak4 )

        y5 = wb.profiles[k].peak_wb_areas[5] * iwater[5]
        axes[axes_list[ax_idx]].plot(x, y5, **style_peak5 )

        ax_idx = ax_idx + 1

axes[1].legend()
fig.savefig('../../../CpxPaper/figures/CpxProfilesWater.png', dpi=100)

#%% Main profiles figure + initial
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5))
ylabelloc = 'left'
axes[15].set_ylabel('Kunlun\ndiopside\n"initial"\n696 $\degree$C\n2 hr', rotation=0, ha=ylabelloc, va='center')
axes[16].set_ylabel('Kunlun\ndiopside\n812 $\degree$C\n6 days', rotation=0, ha=ylabelloc, va='center')
axes[17].set_ylabel('Kunlun\ndiopside\n904 $\degree$C\n154 hr', rotation=0, ha=ylabelloc, va='center')
axes[18].set_ylabel('Kunlun\ndiopside\n1000 $\degree$C\n75 hr', rotation=0, ha=ylabelloc, va='center')
axes[19].set_ylabel('Jaipur\ndiopside\n904 $\degree$C\n0.6 hr', rotation=0, ha=ylabelloc, va='center')

axes_list = range(0, 15)
ax_idx = 0

wbs5 = [my_spectra.K3wb_init,
        my_spectra.K3wb_6days,
        my_spectra.K4wb_154hr,        
        my_spectra.K5wb,
        my_spectra.J1wb,]
        
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
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig5.eps', 
            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig5', dpi=300)

#%% Most things by peak rather than sample
#tit = 'Peak area in Kunlun diopside after 6 days at 816 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 9), tit=None)

ylabelloc = 'left'
axes[15].set_ylabel('3645\ncm$^{-1}$', rotation=0, ha=ylabelloc, va='center')
axes[16].set_ylabel('3617\ncm$^{-1}$', rotation=0, ha=ylabelloc, va='center')
axes[17].set_ylabel('3540\ncm$^{-1}$', rotation=0, ha=ylabelloc, va='center')
axes[18].set_ylabel('3443\ncm$^{-1}$', rotation=0, ha=ylabelloc, va='center')
axes[19].set_ylabel('3350\ncm$^{-1}$', rotation=0, ha=ylabelloc, va='center')

wbs = [
       my_spectra.K3wb_trueInit,
       my_spectra.K3wb_init,
       my_spectra.K3wb_800C_15hr,
       my_spectra.K3wb_6days,

       my_spectra.K4wb_init,
       my_spectra.K4wb_quench,
       my_spectra.K4wb_1hr,
       my_spectra.K4wb_154hr,

       my_spectra.K5wb_init,
       my_spectra.K5wb_75hr,
       ]

styles = [style_init_K3,
          style_K3q,
          style_K3_15h,
          style_K3_6d,

          style_init_K4,
          style_K4q,
          style_K4_1h,
          style_K4_154,

          style_init_K5,
          style_K5,
          ]

#### main plotting ####
ax_idx = 0
for peak_idx in [0, 1, 2, 4, 5]:
    for k in range(3):
        ax = axes[ax_idx]
        L = wb.profiles[k].len_microns
        x = np.array(wb.profiles[k].positions_microns) / L
        print len(x)

        s = ''.join((string.ascii_uppercase[ax_idx],'.'))
        ax.text(0.1, ax.get_ylim()[1] - 0.2*ax.get_ylim()[1], s)

        for wb_idx in range(len(wbs)):
            xi = np.array(wbs[wb_idx].profiles[k].positions_microns) / L
            yi = wbs[wb_idx].profiles[k].peak_areas[peak_idx, :]
            ax.plot(xi, yi, **styles[wb_idx])

        ax_idx = ax_idx + 1

fig.savefig('../../../CpxPaper/figures/CpxProfiles_AllPeaks.png', dpi=200)

#%% K3 Diffusivity modeling
wb = my_spectra.K3wb_6days
fname = 'K3'

iareas = my_spectra.i_K3
peak_D = my_spectra.D_K3
er = my_spectra.e_K3
tops = [17, 35, 30, 40, 15]

style_K3_6d['alpha'] = 1.
style_K3_6d['mew'] = 1.5

wbs = [
       my_spectra.K3wb_init,
       my_spectra.K3wb_800C_15hr,
       my_spectra.K3wb_6days,
       my_spectra.K4wb_quench,
       my_spectra.K4wb_1hr,
       ]

tit = 'Peak area in Kunlun diopside after 6 days at 816 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)

plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er, sidelabels=False)

peak_idx_list = ['3645', '3617', '3540', '3443', '3350']
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
#fig.savefig('../../../CpxPaper/figures/CpxProfiles_K3.png', dpi=200)

#%% K4 91 hr Diffusivity modeling
style_K4_91['alpha'] = 1.
style_K4_91['mew'] = 1.5

iareas = my_spectra.i_K91
peak_D = my_spectra.D_K91

er = my_spectra.e_K91
tops = [20, 35, 30, 35., 15]

wb = my_spectra.K4wb_91hr
wbs = [
       my_spectra.K4wb_init,
       my_spectra.K3wb_init,
       my_spectra.K3wb_800C_15hr,
       my_spectra.K4wb_quench,
       my_spectra.K4wb_1hr,
       my_spectra.K4wb_91hr
       ]

tit = 'Peak area in Kunlun diopside after 91 hr at 904 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)
plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er)
fig.savefig('../../../CpxPaper/figures/CpxProfiles_K4_91.png', dpi=200)

#%% K4 154 hr Diffusivity modeling
tit = 'Peak area in Kunlun diopside after 154 hr at 904 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)

style_K4_154['alpha'] = 1.
style_K4_154['mew'] = 1.5

wbs = [
       my_spectra.K4wb_init,
       my_spectra.K3wb_init,
       my_spectra.K3wb_800C_15hr,
       my_spectra.K4wb_quench,
       my_spectra.K4wb_1hr,
       my_spectra.K4wb_91hr,
       my_spectra.K4wb_154hr,
       ]

iareas = my_spectra.i_K154
peak_D = my_spectra.D_K154

er = my_spectra.e_K154
tops = [21, 36, 31, 35, 16]
wb = my_spectra.K4wb_154hr

plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er)
fig.savefig('../../../CpxPaper/figures/CpxProfiles_K4_154.png', dpi=200)

#%% K5
tit = 'Peak area in Kunlun diopside after 75 hr at 1000 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)

wb = my_spectra.K5wb_75hr

style_K5['alpha'] = 1.
style_K5['mew'] = 1.5

wbs = [
       my_spectra.K5wb_init,
       my_spectra.K3wb_init,
       my_spectra.K3wb_800C_15hr,
       my_spectra.K4wb_quench,
       my_spectra.K4wb_1hr,
       my_spectra.K5wb_75hr
       ]

iareas = my_spectra.i_K5
peak_D = my_spectra.D_K5
er = my_spectra.e_K5
tops = [23, 35, 30, 35, 15]

plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er)
fig.savefig('../../../CpxPaper/figures/CpxProfiles_K5.png', dpi=200)

#%% J diffusivity modeling

wb = my_spectra.J1wb

style_J['alpha'] = 0.5
style_J['mew'] = 1.5

wbs = [
       my_spectra.J1wb_initial,
       my_spectra.J1wb
       ]

iareas = my_spectra.i_J
peak_D = my_spectra.D_J
er = my_spectra.e_J
tops = [4, 8, 23, 20, 50]
peak_D_slow = np.array(peak_D) - 1.
ytickgrid = [1, 2, 5, 5, 10]

tit = 'Peak area in Jaipur diopside after 0.6 hr at 904 $\degree$C (cm$^{-2}$)'
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5), tit=tit)
plotpeaks(wb_list=[wb]*5, wbs_list=[wbs]*5, ilist=iareas, dlist=peak_D,
          elist=er, slowb=[True]*5, dlist_slow=peak_D_slow, legidx=4,
          dlabel='Diffusion curve', ytickgrid=ytickgrid)
fig.savefig('../../../CpxPaper/figures/CpxProfiles_J.png', dpi=200)

#%% Bulk WB areas with diffusivity estimates
wb_list = [my_spectra.K3wb_6days,
           my_spectra.K4wb_91hr,
           my_spectra.K4wb_154hr,
           my_spectra.J1wb,
           my_spectra.K5wb_75hr]

typelabels = ['Kunlun\ndiopside\n812 $\degree$C\n6 days',
              'Kunlun\ndiopside\n904 $\degree$C\n91 hr',
              'Kunlun\ndiopside\n904 $\degree$C\n154 hr',
              'Jaipur\ndiopside\n904 $\degree$C\n0.6 hr',
              'Kunlun\ndiopside\n1000 $\degree$C\n75 hr']

my_ilist = [my_spectra.i_K3, my_spectra.i_K91, my_spectra.i_K154,
            my_spectra.i_J, my_spectra.i_K5]
my_dlist = [my_spectra.D_K3, my_spectra.D_K91, my_spectra.D_K154,
            my_spectra.D_J, my_spectra.D_K5]
my_elist = [my_spectra.e_K3, my_spectra.e_K91, my_spectra.e_K154,
            my_spectra.e_J, my_spectra.e_K5]

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

peak_idx_list = ['3645', '3617', '3540', '3443', '3350']
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


plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig6.eps', 
            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig6', dpi=300)
#plt.savefig('Fig6.eps', format='eps', dpi=1000)

#%% Bulk WB areas with diffusivity estimates + additional symbols
wb_list = [my_spectra.K3wb_6days,
           my_spectra.K4wb_91hr,
           my_spectra.K4wb_154hr,
           my_spectra.J1wb,
           my_spectra.K5wb_75hr]

typelabels = ['Kunlun\ndiopside\n812 $\degree$C\n6 days',
              'Kunlun\ndiopside\n904 $\degree$C\n91 hr',
              'Kunlun\ndiopside\n904 $\degree$C\n154 hr',
              'Jaipur\ndiopside\n904 $\degree$C\n0.6 hr',
              'Kunlun\ndiopside\n1000 $\degree$C\n75 hr']

my_ilist = [my_spectra.i_K3, my_spectra.i_K91, my_spectra.i_K154,
            my_spectra.i_J, my_spectra.i_K5]
my_dlist = [my_spectra.D_K3, my_spectra.D_K91, my_spectra.D_K154,
            my_spectra.D_J, my_spectra.D_K5]
my_elist = [my_spectra.e_K3, my_spectra.e_K91, my_spectra.e_K154,
            my_spectra.e_J, my_spectra.e_K5]

tops = [ 1.6, 1.6, 1.6, 1.8, 1.6,]
slowb_list = [False, False, False, True, False, False]
wbs_list = [wbsK3, wbsK4_91, wbsK4_154, wbsJ, wbsK5]

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

## side labels
ylabelloc = 'left'
formatter = '{:.2f}'
subscript = '$_x$'
for k in range(5):
    axes[15+k].set_ylabel(''.join((typelabels[k], '\nlog$_{10}$D$_c$\n',
                        formatter.format(Dxz[k]), '\n+/-',
                        formatter.format(er[k]))),
                        rotation=0, ha=ylabelloc, va='center', fontsize=10)

fig.savefig('../../../CpxPaper/figures/CpxProfiles_bulkH_WB.png', dpi=200)
fig.savefig('../../../CpxPaper/figures/CpxProfiles_bulkH_WB.svg')

#%% Bulk WB areas with diffusivity estimates
wb_list = [my_spectra.K3wb_6days,
           my_spectra.K4wb_91hr,
           my_spectra.K4wb_154hr,
           my_spectra.J1wb,
           my_spectra.K5wb_75hr]

typelabels = ['Kunlun\ndiopside\n812 $\degree$C\n6 days',
              'Kunlun\ndiopside\n904 $\degree$C\n91 hr',
              'Kunlun\ndiopside\n904 $\degree$C\n154 hr',
              'Jaipur\ndiopside\n904 $\degree$C\n0.6 hr',
              'Kunlun\ndiopside\n1000 $\degree$C\n75 hr']

my_ilist = [my_spectra.i_K3, my_spectra.i_K91, my_spectra.i_K154,
            my_spectra.i_J, my_spectra.i_K5]
my_dlist = [my_spectra.D_K3, my_spectra.D_K91, my_spectra.D_K154,
            my_spectra.D_J, my_spectra.D_K5]
my_elist = [my_spectra.e_K3, my_spectra.e_K91, my_spectra.e_K154,
            my_spectra.e_J, my_spectra.e_K5]

tops = [ 1.6, 1.6, 1.6, 1.8, 1.6,]
slowb_list = [False, False, False, True, False, False]
wbs_list = []
for wb in wb_list:
    wbs_list.append([wb])

### Moving Jaipur to the end for all of them
#for li in [wb_list, typelabels, my_ilist, my_dlist, my_elist,
#           tops, slowb_list, wbs_list]:
#    li.insert(4, li.pop(3))

iareas = np.ones(5)
Dxz = np.ones(5)
er = np.ones(5)
for idx in range(5):
    iareas[idx] = my_ilist[idx][-1]
    Dxz[idx] = my_dlist[idx][-1]
    er[idx] = my_elist[idx][-1]
peak_D_slow = np.array(Dxz) - 1.


tit = ''.join(('Bulk H (Total area / Initial total area)',
               '\n1=~29 ppm H$_2$O in Kunlun and ~30 ppm H$_2$O in Jaipur diopside'))
fig, axes = figOutline(nrow=5, ncol=3, figsize=(6.5, 7.5),
                       top=1.6, tit=tit)
plotpeaks(wb_list=wb_list, wbs_list=wbs_list, ilist=iareas, dlist=Dxz,
          elist=er, slowb=slowb_list, dlist_slow=peak_D_slow,
          show_legend=False, sidelabels=False,
          wholeblock=True, peak_idx_list=[None]*5, ytickgrid=[0.5]*5)

## side labels
ylabelloc = 'left'
formatter = '{:.2f}'
subscript = '$_x$'
for k in range(5):
    axes[15+k].set_ylabel(''.join((typelabels[k], '\nlog$_{10}$D$_c$\n',
                        formatter.format(Dxz[k]), '\n+/-',
                        formatter.format(er[k]))),
                        rotation=0, ha=ylabelloc, va='center', fontsize=10)

fig.savefig('../../../CpxPaper/figures/CpxProfiles_bulkH_poster.svg')
