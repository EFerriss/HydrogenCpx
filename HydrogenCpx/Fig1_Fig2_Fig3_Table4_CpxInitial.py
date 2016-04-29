# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 13:59:16 2015

@author: Ferriss
H diffusion in cpx project

Create 3x3 subplots of initial polarized FTIR spectra for Kunlun, Jaipur, PMR
"""
import cpx_spectra as sp
from pynams import pynams as nams
from pynams import styles as styles
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import matplotlib.lines as mlines
import csv

plt.close('all')
plt.style.use('paper')
numformat = '{:.0f}'

# spectra of interest
Kunlun_list = [sp.ave_K6_Ea, sp.ave_K6_Eb, sp.ave_K6_Ec]
Jaipur_list = [sp.J_Ea, sp.J_Eb, sp.J_Ec] # These are from Woods et al.
J2_list = [sp.J2_Ea, sp.J2_Eb, sp.J2_Ec] # My Jaipur diopside sample J2
PMR_list = [sp.PMR_Ea, sp.PMR_Eb, sp.PMR_Ec]
spec_list = Kunlun_list + J2_list + PMR_list

#%% Determine baselines, areas and bulk water estimates for each sample

# Ultimately make three estimates and take the average and stdev.
water_estimates = []

def make_water_estimates():
    """ Take current baselines for the 9 spectra in order
    Kunlun, Jaipur, PMR and produce water estimates
    for each sample by adding up all the areas and applying Bell's
    calibration for cpx
    """
    areas = np.zeros(9)
    for idx, spec in enumerate(spec_list):
        areas[idx] = spec.area_under_curve(require_saved_baseline=False, 
                                           show_plot=False, printout=False)
    
    Kunlun_areas = areas[0:3]
    Jaipur_areas = areas[3:6]
    PMR_areas = areas[6:9]
    
    Kunlun_area_total = np.sum(Kunlun_areas)
    Jaipur_area_total = np.sum(Jaipur_areas)
    PMR_area_total = np.sum(PMR_areas)
    
    Kunlun_water = nams.area2water(Kunlun_area_total, phase='cpx', calibration='Bell')
    Jaipur_water = nams.area2water(Jaipur_area_total, phase='cpx', calibration='Bell')
    PMR_water = nams.area2water(PMR_area_total, phase='cpx', calibration='Bell')
    
    return Kunlun_water.n, Jaipur_water.n, PMR_water.n

# First make 3 baselines for each spectrum
# the wavenumber ranges remain constant throughout
wn_high_list = [3775]*6 + [3700, 3775., 3700.]
wn_low_list = [3150., 3150., 3150., 3100., 3150., 3100., 3100., 3100., 3100.]

# the extent of the shift for each quadratic changes
shift_list_main = [0.05, 0.08, 0.03, 0.13, 0.055, 0.1,  1.5, 0.5, 0.4]
shift_list_low =  [0.07, 0.10, 0.05, 0.2,  0.08,  0.2,  2.,  0.8, 0.6]
shift_list_high = [0.03, 0.06, 0.01, 0.06, 0.11,  0.01, 1.1, 0.2, 0.2]

sp.PMR_Ea.base_mid_wn=3500.

# make central baselines and water estimates
for idx, spec in enumerate(spec_list):
    spec.bline = []
    spec.bline.append(spec.make_baseline(show_plot=False, 
                                         linetype='quadratic', 
                                         wn_low=wn_low_list[idx], 
                                         wn_high=wn_high_list[idx], 
                                         shiftline=shift_list_main[idx]))
    spec.save_baseline() # save main baselines to file for peak fitting
water_estimates.append(make_water_estimates())

# make low baselines and water estimates
for idx, spec in enumerate(spec_list):
    spec.bline.append(spec.make_baseline(show_plot=False, linetype='quadratic', 
                                         wn_low=wn_low_list[idx], 
                                         wn_high=wn_high_list[idx], 
                                         shiftline=shift_list_low[idx]))
water_estimates.append(make_water_estimates())
                       
# make high baselines and water estimates
for idx, spec in enumerate(spec_list):
    spec.bline.append(spec.make_baseline(show_plot=False, linetype='quadratic', 
                                         wn_low=wn_low_list[idx], 
                                         wn_high=wn_high_list[idx], 
                                         shiftline=shift_list_high[idx]))
water_estimates.append(make_water_estimates())

# Get mean and standard deviation for each
K_water_n = np.mean(np.array(water_estimates)[:, 0])
J2_water_n = np.mean(np.array(water_estimates)[:, 1])
PMR_water_n = np.mean(np.array(water_estimates)[:, 2])

K_water_s = np.std(np.array(water_estimates)[:, 0])
J2_water_s = np.std(np.array(water_estimates)[:, 1])
PMR_water_s = np.std(np.array(water_estimates)[:, 2])

K_water = ufloat(K_water_n, K_water_s)
J2_water = ufloat(J2_water_n, J2_water_s)
PMR_water = ufloat(PMR_water_n, PMR_water_s)

### Figure showing 3 baselines
ax = nams.plotsetup_3x3(yhi=1.5)

# wavenumber range to plot
high = 3800
low = 3200

elabels = ['E || a*', 'E || b $(\\beta)$', 'E || c',
           'E || a*', 'E || b $(\\beta)$', 'E || c',
           'E || $\\alpha$', 'E || b $(\\beta)$', 'E || $\gamma$']

PMRfac = 10 # how much to divide PMR spectra absorbances by

abs_shift = 0.15 # how much to shift spectra so baselines are clear

for k, spec in enumerate(spec_list):
    if k > 5: # PMR
        ax[k].plot(spec.wn_full, spec.abs_full_cm / PMRfac + abs_shift, 
                   **styles.style_spectrum)
        for x in xrange(3): # loop through each baseline
            ax[k].plot(spec.base_wn, spec.bline[x]/PMRfac + abs_shift,
                        **styles.style_baseline)
    else:
        ax[k].plot(spec.wn_full, spec.abs_full_cm + abs_shift,
                   **styles.style_spectrum)
        for x in xrange(3):
            ax[k].plot(spec.base_wn, spec.bline[x] + abs_shift, 
                        **styles.style_baseline)
                        
    ax[k].text(3030, 1.2, elabels[k], horizontalalignment='right',
                        backgroundcolor='w', fontsize=11)
 
ax[0].set_ylabel('Kunlun diopside\n%s ppm H$_2$O' % numformat.format(K_water))
ax[3].set_ylabel('absorbance (cm$^{-1}$)\n' +
                    'Jaipur diopside\n%s ppm H$_2$O' % 
                    numformat.format(J2_water))
ax[6].set_ylabel('augite PMR-53 / %i\n%s ppm H$_2$O' % 
                  (PMRfac, numformat.format(PMR_water)))

plt.savefig('Fig2.png', format='png', dpi=300)
plt.savefig('Fig2.eps', format='eps')

# %% Figure showing baseline-subtracted spectra and fit peaks
ax = nams.plotsetup_3x3(xhi=high, xlo=low, xtickgrid=200, yhi=1.5,
                        fig_size_inches=(6.5, 6.3))

peakpos_diopside = [3645, 3617, 3540, 3443, 3355]
peakpos_J2 = [3675, 3645, 3617, 3540, 3460, 3355]
peakpos_augite = [3620, 3550, 3460, 3355]
shrinker = 0.1 # for squeezing in legend at the bottom
PMRfac = 6 # Divide PMR spectra

for k, spec in enumerate(spec_list):
    spec.get_baseline()
    xdata = spec.base_wn
    bdata = spec.abs_nobase_cm
    gaussian, summed_spectrum = spec.get_peakfit()
        
    if k > 5:
        y = bdata / PMRfac
        ysum = summed_spectrum / PMRfac
        peakpos = peakpos_augite
    elif k > 2:
        y = bdata
        ysum = summed_spectrum        
        peakpos = peakpos_J2        
    else:
        y = bdata 
        ysum = summed_spectrum        
        peakpos = peakpos_diopside

    # plot 
    ax[k].plot(xdata, y, **styles.style_spectrum)
    ax[k].plot(xdata, ysum, **styles.style_summed)

    for x in range(len(gaussian)):
        if k > 5:
            ax[k].plot(xdata, gaussian[x]/PMRfac, **styles.style_fitpeak)
        else:            
            ax[k].plot(xdata, gaussian[x], **styles.style_fitpeak)
            
    # label ray path
    ax[k].text(3230, 1.2, elabels[k], horizontalalignment='right',
                backgroundcolor='w', fontsize=11)

    # label peaks
    ax[k].grid(b=False)
    for peakwn in peakpos:
        idx = (np.abs(xdata-peakwn)).argmin()
        ytext = y[idx]
#        if peakwn == 3617:
#            xtext = 3600
#        else:
        xtext = peakwn

        ax[k].text(xtext, ytext, ''.join(('$\leftarrow$', str(peakwn))), 
                   rotation=90, ha='center', va='bottom', fontsize=8)

    # Shrink to make room for legend at the bottom
    box = ax[k].get_position()
    ax[k].set_position([box.x0, box.y0 + box.height*shrinker, 
                     box.width, box.height*(1.0-shrinker)])

ax[0].set_ylabel('Kunlun diopside\n%s ppm H$_2$O' % numformat.format(K_water))
ax[3].set_ylabel('absorbance (cm$^{-1}$)\n' +
                    'Jaipur diopside\n%s ppm H$_2$O' % 
                    numformat.format(J2_water))
ax[6].set_ylabel('augite PMR-53 / %i\n%s ppm H$_2$O' % 
                  (PMRfac, numformat.format(PMR_water)))

add_marker1 = mlines.Line2D([], [], label='observed spectrum', 
                            **styles.style_spectrum)
add_marker2 = mlines.Line2D([], [], label='Gaussians', **styles.style_fitpeak)
add_marker3 = mlines.Line2D([], [], label='sum of Gaussians', 
                            **styles.style_summed)
leg_handles = [add_marker1, add_marker2, add_marker3]


main_legend = plt.legend(handles=leg_handles, ncol=3,
                          bbox_to_anchor=(250, 0, 1500, -0.75),
                          bbox_transform=ax[6].transData, 
                          frameon=True)

#plt.savefig('C:\\Users\\Ferriss\\Documents\\CpxPaper\\Fig3.tif')
plt.savefig('Fig3.tif', format='tif', dpi=300)
plt.savefig('Fig3.eps', format='eps')

# %% Now sum up the peaks for each one
ax = nams.plotsetup_3stacked(yhi=2, ytickgrid=0.5)

cpx_calib_Bell95 = 1.0 / ufloat(7.08, 0.32)
xh = 3800
xl = 3200
dx = xh-xl
xplot = np.linspace(xl, xh, 150)

def sum_and_plot_peaks(rlist, axplot, fac=1, specs=spec_list):
    """ Requires that each direction have the same number of peaks for each sample.
    """
    if len(rlist) != 3:
        print 'rlist is a list of exactly 3 indices in list of spectra in specs'
        return
        
    npeaks = len(specs[rlist[0]].peakpos) # number of peaks
    
    ssum = np.zeros((npeaks, len(xplot))) 
    bigsum = np.zeros_like(xplot)
    sum_areas = np.zeros(npeaks)
    
    # for determining polarization of peaks
    peakareas = np.zeros((3, npeaks))
    peakareas_total = np.zeros(npeaks)
    polarization = np.zeros((npeaks, 3))

    for idx, k in enumerate(rlist): # loop through three spectra with indices given in rlist
        for m in range(npeaks): # loop through each peak
            position = specs[k].peakpos[m]
            height = specs[k].peak_heights[m]
            width = specs[k].peak_widths[m]
            peakareas[idx,m] = specs[k].peak_areas[m]
            gauss = nams.make_gaussian(position, height, width, x=xplot)
            ssum[m] = ssum[m] + gauss # sum of all gaussians
            
    for m in range(npeaks):
        ax[axplot].plot(xplot, ssum[m]/fac, **styles.style_fitpeak)
        bigsum = bigsum + ssum[m]
        dy = np.mean(ssum[m])
        sum_areas[m] = dx*dy
        
        # polarization
        peakareas_total[m] = np.sum(peakareas[:,m])
        
        for idx in xrange(3):
            polarization[m, idx] = peakareas[idx,m] / peakareas_total[m]

    ax[axplot].plot(xplot, bigsum/fac, **styles.style_summed)
    return sum_areas, polarization

ax[0].set_ylabel('absorbance (cm$^{-1}$)')
ax[1].set_ylabel('absorbance (cm$^{-1}$)')
ax[2].set_ylabel(''.join(('absorbance / ', '{:.0f}'.format(PMRfac), 
                    ' (cm$^{-1}$)')))

Kunlun_peak_areas, Kpolar = sum_and_plot_peaks(range(3), axplot=0)
Jaipur_peak_areas, Jpolar = sum_and_plot_peaks(range(3, 6), axplot=1)
PMR_peak_areas, Ppolar = sum_and_plot_peaks(range(6, 9), axplot=2, fac=PMRfac)

print np.sum(Kunlun_peak_areas*cpx_calib_Bell95)
print np.sum(Jaipur_peak_areas*cpx_calib_Bell95)
print np.sum(PMR_peak_areas*cpx_calib_Bell95)

plt.savefig('Fig1_KJP_3stack.svg', dpi=300)

#%% print and save results of summing
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

def csv_peaks(spec, arealist, water, csvfi, polar):
    """Takes spectrum with peak position information, peak-specific areas, 
    initial water, and polarization information and saves to .CSV file"""
    npeaks = len(spec.peakpos)
    for m in range(npeaks):
        frac = arealist[m] / np.sum(arealist)
        w = frac*water
        csvfi.writerow([spec.peakpos[m], spec.peak_widths[m], spec.peak_heights[m],
                       arealist[m], frac*100.0, w, 
                       polar[m][0], polar[m][1], polar[m][2]])

# print to screen            
print 'Kunlun'
print 'Total water estimate:', K_water, 'ppm H2O'
print_peak_water(spec = Kunlun_list[0], 
                 arealist = Kunlun_peak_areas, water=K_water)
print ' '
print 'Jaipur'
print 'Total water estimate:', J2_water, 'ppm H2O'
print_peak_water(spec = J2_list[0], 
                 arealist = Jaipur_peak_areas, water=J2_water)
print ' '
print 'PMR'
print 'Total water estimate:', PMR_water, 'ppm H2O'
print_peak_water(spec = PMR_list[0], 
                 arealist=PMR_peak_areas, water=PMR_water)

# save to .CSV
savefile = 'Table3_cpx_peakfit_intial.CSV'
with open(savefile, 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    
    # Kunlun
    spamwriter.writerow(['Kunlun diopside', '', '', '', '', '', 
                         'relative polarization'])
    spamwriter.writerow(['Bulk water', K_water, 'ppm H2O'])
    spamwriter.writerow(['peak position (cm-1)',  'width (cm-1)', 
                         'height (cm-1)', 'area (cm-2)', 
                         '% bulk water', 'water (ppm H2O)', 
                         'alpha', 'beta', 'gamma'])
    csv_peaks(spec=Kunlun_list[0], arealist=Kunlun_peak_areas, 
              water=K_water, csvfi=spamwriter, 
              polar=Kpolar)

    # Jaipur
    spamwriter.writerow(['Jaipur diopside'])
    spamwriter.writerow(['Bulk water', J2_water, 'ppm H2O'])
    csv_peaks(spec=J2_list[0], arealist=Jaipur_peak_areas,
              water=J2_water, csvfi=spamwriter, polar=Jpolar)

    # PMR              
    spamwriter.writerow(['augite PMR-53'])
    spamwriter.writerow(['Bulk water', PMR_water, 'ppm H2O'])
    csv_peaks(spec=PMR_list[0], arealist=PMR_peak_areas,
              water=PMR_water, csvfi=spamwriter, polar=Ppolar)
