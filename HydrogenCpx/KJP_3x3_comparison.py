# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 13:59:16 2015

@author: Ferriss
H diffusion in cpx project

Create 3x3 subplots of initial polarized FTIR spectra for Kunlun, Jaipur, PMR
"""
import my_spectra as sp
from pynams import pynams as nams
from pynams import styles as styles
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import matplotlib.lines as mlines
import csv

plt.close('all')
plt.style.use('paper')
savedpi = 300
numformat = '{:.0f}'

# wavenumber range to plot
high = 3800
low = 3200

# %% First make water estimates for all, 
# saving all baselines as -3baselines.CSV
Kunlun_list = [sp.ave_K6_Ea, sp.ave_K6_Eb, sp.ave_K6_Ec]
K_area, K_water = nams.water_from_spectra(Kunlun_list, proper3=True, 
                                        savebaselines=False, show_plots=False)


Jaipur_list = [sp.J_Ea, sp.J_Eb, sp.J_Ec]
J_area, J_water = nams.water_from_spectra(Jaipur_list, proper3=True, 
                                        savebaselines=False, show_plots=False)


PMR_list = [sp.PMR_Ea, sp.PMR_Eb, sp.PMR_Ec]
PMR_area, PMR_water = nams.water_from_spectra(PMR_list, proper3=True,
                                            savebaselines=False, 
                                            show_plots=False)

spec_list = Kunlun_list + Jaipur_list + PMR_list

# Save main baselines in -baseline.CSV for use in Matlab peakfitting
for x in spec_list:
    # Without this, it just saves the most recent (i.e., too low) baseline
    x.make_baseline(2)
    x.save_baseline(2)

for x in spec_list:
    print x.fname
    
elabels = ['E || a*', 'E || b $(\\beta)$', 'E || c',
           'E || a', 'E || b $(\\beta)$', 'E || c*',
           'E || $\\alpha$', 'E || b $(\\beta)$', 'E || $\gamma$']

# %%  Figure showing 3 baselines
ax = nams.plotsetup_3x3()

PMRfac = 15

abs_shift = 0.12

for k in range(9):
    bdata = spec_list[k].get_3baselines()
    if k > 5:
        ax[k].plot(spec_list[k].wn, spec_list[k].abs_cm / PMRfac + abs_shift, 
                        **styles.style_spectrum)
        for x in range(1, 4):
            ax[k].plot(bdata[:,0], bdata[:,x]/PMRfac + abs_shift,
                        **styles.style_baseline)
    else:
        ax[k].plot(spec_list[k].wn, spec_list[k].abs_cm + abs_shift,
                        **styles.style_spectrum)
        for x in range(1, 4):
            ax[k].plot(bdata[:,0], bdata[:,x] + abs_shift, 
                        **styles.style_baseline)
    ax[k].text(3030, 0.75, elabels[k], horizontalalignment='right',
                        backgroundcolor='w', fontsize=11)

ax[0].set_ylabel('Kunlun diopside\n%s ppm H$_2$O' % numformat.format(K_water))
ax[3].set_ylabel('absorbance (cm$^{-1}$)\n' +
                    'Jaipur diopside\n%s ppm H$_2$O' % 
                    numformat.format(J_water))
ax[6].set_ylabel('augite PMR-53 / %i\n%s ppm H$_2$O' % 
                  (PMRfac, numformat.format(PMR_water)))

#savefile = sp.default_savefolder + 'KJP_3x3_baselines'
#plt.savefig(savefile, dpi=savedpi)

plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig2.eps', 
            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig2', dpi=300)
plt.savefig('Fig2.eps', format='eps', dpi=1000)

# %% Figure showing baseline-subtracted spectra and fit peaks
ax = nams.plotsetup_3x3(xhi=high, xlo=low, xtickgrid=200, yhi=1.5,
                        fig_size_inches=(6.5, 6.3))

peakpos_diopside = [3645, 3617, 3540, 3460, 3350]
peakpos_augite = [3620, 3550, 3460, 3355]
shrinker = 0.1 # for squeezing in legend at the bottom
PMRfac = 6 # Divide PMR spectra

for k in range(9):
    bdata = spec_list[k].get_3baselines()   
    gaussian, summed_spectrum = spec_list[k].get_peakfit()
        
    if k > 5:
        y = bdata[:,5] / PMRfac
        ysum = summed_spectrum / PMRfac
        peakpos = peakpos_augite
    else:
        y = bdata[:,5]
        ysum = summed_spectrum        
        peakpos = peakpos_diopside

    ax[k].plot(bdata[:,0], y, **styles.style_spectrum)
    ax[k].plot(bdata[:,0], ysum, **styles.style_summed)

    for x in range(len(gaussian)):
        if k > 5:
            ax[k].plot(bdata[:,0], gaussian[x]/PMRfac, **styles.style_fitpeak)
        else:            
            ax[k].plot(bdata[:,0], gaussian[x], **styles.style_fitpeak)
            
    # label ray path
    ax[k].text(3230, 1.2, elabels[k], horizontalalignment='right',
                backgroundcolor='w', fontsize=11)

    # label peaks
    ax[k].grid(b=False)
    for peakwn in peakpos:
        idx = (np.abs(bdata[:,0]-peakwn)).argmin()
        ypeak = y[idx]
        xpeak = peakwn
        if peakwn == 3617:
            xtext = 3600
        else:
            xtext = xpeak
        ax[k].annotate(str(peakwn), xy=(xpeak, ypeak), 
                       xytext=(xtext, ypeak+0.2), fontsize=8,
                        rotation=90, ha='center', va='bottom',
                        arrowprops=dict(facecolor='black', arrowstyle='->',
                                        linewidth=1))

#        ax[k].text(xpeak, ypeak, ''.join(('$\leftarrow$', str(peakwn))), 
#                    rotation=90, ha='center', va='bottom', fontsize=8)

    # Shrink to make room for legend at the bottom
    box = ax[k].get_position()
    ax[k].set_position([box.x0, box.y0 + box.height*shrinker, 
                     box.width, box.height*(1.0-shrinker)])

ax[0].set_ylabel('Kunlun diopside\n%s ppm H$_2$O' % numformat.format(K_water))
ax[3].set_ylabel('absorbance (cm$^{-1}$)\n' +
                    'Jaipur diopside\n%s ppm H$_2$O' % 
                    numformat.format(J_water))
ax[6].set_ylabel('augite PMR-53 / %i\n%s ppm H$_2$O' % 
                  (PMRfac, numformat.format(PMR_water)))

add_marker1 = mlines.Line2D([], [], label='observed spectrum', 
                            **styles.style_spectrum)
add_marker2 = mlines.Line2D([], [], label='Gaussians', **styles.style_fitpeak)
add_marker3 = mlines.Line2D([], [], label='sum of Gaussians', 
                            **styles.style_summed)
leg_handles = [add_marker1, add_marker2, add_marker3]


main_legend = plt.legend(handles=leg_handles, ncol=3,
#                          bbox_to_anchor=(1, 0, 0, -0.5), 
                          bbox_to_anchor=(250, 0, 1500, -0.75),
                          bbox_transform=ax[6].transData, 
                          frameon=True)

#savefile = sp.default_savefolder + 'KJP_3x3_peakfit'
#plt.savefig(savefile, dpi=savedpi)

plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig3.eps', 
            format='eps', dpi=1000)
plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig3', dpi=300)
#plt.savefig('Fig3.eps', format='eps', dpi=1000)
# %% Now sum up the peaks for each one

ax = nams.plotsetup_3stacked(yhi=2, ytickgrid=0.5)

cpx_calib_Bell95 = 1.0 / ufloat(7.08, 0.32)
xh = 3800
xl = 3200
dx = xh-xl
xplot = np.linspace(xl, xh, 150)

def sum_and_plot_peaks(rlist, axplot, fac=1):
    npeaks = len(spec_list[rlist[0]].peakpos)
    ssum = np.zeros((npeaks, len(xplot)))
    bigsum = np.zeros_like(xplot)
    sum_areas = np.zeros(npeaks)
    for k in rlist:
        for m in range(npeaks):
            position = spec_list[k].peakpos[m]
            height = spec_list[k].peak_heights[m]
            width = spec_list[k].peak_widths[m]
            gauss = nams.make_gaussian(position, height, width, x=xplot)
            ssum[m] = ssum[m] + gauss
    for m in range(npeaks):
        ax[axplot].plot(xplot, ssum[m]/fac, **styles.style_fitpeak)
        bigsum = bigsum + ssum[m]
        dy = np.mean(ssum[m])
        sum_areas[m] = dx*dy
    ax[axplot].plot(xplot, bigsum/fac, **styles.style_summed)
    return sum_areas

### version 3
#ax[0].set_ylabel('Kunlun diopside\n%s ppm H$_2$O' % numformat.format(K_water))
#ax[1].set_ylabel('Total polarized FTIR absorbance ' + 
#                    'over 3 directions (cm$^{-1})$\n' +
#                    'Jaipur diopside\n%s ppm H$_2$O' % 
#                    numformat.format(J_water))
#ax[2].set_ylabel('augite PMR-53 / %i\n%s ppm H$_2$O' % 
#                  (PMRfac, numformat.format(PMR_water)))

ax[0].set_ylabel('absorbance (cm$^{-1}$)')
ax[1].set_ylabel('absorbance (cm$^{-1}$)')
ax[2].set_ylabel(''.join(('absorbance / ', '{:.0f}'.format(PMRfac), 
                    ' (cm$^{-1}$)')))

Kunlun_peak_areas = sum_and_plot_peaks(range(3), axplot=0)
Jaipur_peak_areas = sum_and_plot_peaks(range(3, 6), axplot=1)
PMR_peak_areas = sum_and_plot_peaks(range(6, 9), axplot=2, fac=PMRfac)

print np.sum(Kunlun_peak_areas*cpx_calib_Bell95)
print np.sum(Jaipur_peak_areas*cpx_calib_Bell95)
print np.sum(PMR_peak_areas*cpx_calib_Bell95)

savefile = sp.default_savefolder + 'KJP_3stack.svg'
plt.savefig(savefile, dpi=300)

savefile = sp.default_savefolder + 'KJP_3stack.png'
plt.savefig(savefile, dpi=300)


# %% print and save results of summing
def print_peak_water(spec, arealist, water):
    """Takes spectrum with peak positions, list of areas, and total bulk 
    water, and prints out peak area, fraction of bulk, and water estimate.
    """
    npeaks = len(spec.peakpos)
    print "'True' bulk water:   ", water, 'ppm H2O'
    for m in range(npeaks):
        frac = arealist[m] / np.sum(arealist)
        w = frac*water
#        print (numformat.format(spec.peakpos[m]), 
#               numformat.format(arealist[m]), numformat.format(frac), w)
        print spec.peakpos[m], 'cm-1,', (
                '{:.0f}'.format(arealist[m])), 'cm-2 area,', (
                '{:.0f}'.format(frac*100.0)), '% bulk,', w, 'ppm H2O'

def csv_peaks(spec, arealist, water, csvfi):
    """Takes spectrum with peak position information, peak-specific areas, 
    and initial water and saves to .CSV file"""
    npeaks = len(spec.peakpos)
    for m in range(npeaks):
        frac = arealist[m] / np.sum(arealist)
        w = frac*water
        csvfi.writerow([spec.peakpos[m], spec.peak_widths[m], spec.peak_heights[m],
                       arealist[m], frac*100.0, w])

# 'True' water contents
#K_water_known = ufloat(29, 4)
J_water_known = J_water
#PMR_water_known = ufloat(268, 8)

# print to screen            
print 'Kunlun'
print 'Total water estimate:', K_water, 'ppm H2O'
print_peak_water(spec = Kunlun_list[0], 
                 arealist = Kunlun_peak_areas, water=sp.K_water_known)
print ' '
print 'Jaipur'
print 'Total water estimate:', J_water, 'ppm H2O'
print_peak_water(spec = Jaipur_list[0], 
                 arealist = Jaipur_peak_areas, water=J_water_known)
print ' '
print 'PMR'
print 'Total water estimate:', PMR_water, 'ppm H2O'
print_peak_water(spec = PMR_list[0], 
                 arealist=PMR_peak_areas, water=sp.PMR_water_known)

# save to .CSV
savefile = 'C:\Users\Ferriss\Documents\CpxPaper\Table1_cpx_peakfit_intial.CSV'
with open(savefile, 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    spamwriter.writerow(['Kunlun diopside'])
    spamwriter.writerow(['Bulk water', sp.K_water_known, 'ppm H2O'])
    spamwriter.writerow(['peak position (cm-1)',  'width (cm-1)', 
                         'height (cm-1)', 'area (cm-2)', 
                         '% bulk water', 'water (ppm H2O)'])
    csv_peaks(spec=Kunlun_list[0], arealist=Kunlun_peak_areas, 
              water=sp.K_water_known, csvfi=spamwriter)

    spamwriter.writerow(['Jaipur diopside'])
    spamwriter.writerow(['Bulk water', J_water_known, 'ppm H2O'])
    csv_peaks(spec=Jaipur_list[0], arealist=Jaipur_peak_areas,
              water=J_water_known, csvfi=spamwriter)
             
    spamwriter.writerow(['augite PMR-53'])
    spamwriter.writerow(['Bulk water', sp.PMR_water_known, 'ppm H2O'])
    csv_peaks(spec=PMR_list[0], arealist=PMR_peak_areas,
              water=sp.PMR_water_known, csvfi=spamwriter)

#%% relative polarization
npeaks = 4
biglist = PMR_list
polar = np.zeros([3, npeaks])

for k in range(3):
    spec = biglist[k]
    gaussian, summed_spectrum = spec.get_peakfit()
    for peak in range(npeaks):
        polar[k][peak] = spec.peak_areas[peak]
        
peaksum = polar.sum(axis=0)

print '\npolarization - check with sample in program'
for peak_idx in range(npeaks):
    a = polar[0][peak_idx] / peaksum[peak_idx]
    b = polar[1][peak_idx] / peaksum[peak_idx]
    c = polar[2][peak_idx] / peaksum[peak_idx]
    print a, b, c