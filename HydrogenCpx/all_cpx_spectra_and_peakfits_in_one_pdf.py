# -*- coding: utf-8 -*-
"""
Created on Thu May 12 16:31:48 2016

@author: Ferriss

Create a pdf with plots showing all baselines and resolved peaks for all
spectra used to determine site-specific hydrogen diffusivities 
in the paper 'Site-specific hydrogen diffusion during
clinopyroxene dehyration' by Ferriss et al. 2016
"""

import cpx_spectra
from matplotlib.backends.backend_pdf import PdfPages


#%%
pp = PdfPages('all_cpx_spectra_and_peakfits.pdf')

wb_list = [cpx_spectra.K3wb_init,
           cpx_spectra.K3wb_700C_19hr,
           cpx_spectra.K3wb_700C_35hr,
           cpx_spectra.K3wb_800C_15hr,
           cpx_spectra.K3wb_6days,
           cpx_spectra.K4wb_init,
           cpx_spectra.K4wb_quench,
           cpx_spectra.K4wb_91hr,
           cpx_spectra.K4wb_154hr,
           cpx_spectra.K5wb_init,
           cpx_spectra.K5wb,
           cpx_spectra.J1wb_initial,
           cpx_spectra.J1wb
           ]

for wb in wb_list:
    for profile in wb.profiles:   
        for idx, spectrum in enumerate(profile.spectra_list):
            spectrum.get_baseline()
            spectrum.get_peakfit()            
            fig, ax = spectrum.plot_peakfit_and_baseline()
            fig.set_size_inches(6, 6)
            title = ''.join((profile.profile_name, '\n ', 
                             '{:.1f}'.format(profile.positions_microns[idx]),
                             ' $\mu$m || ', profile.direction, ', ray path || ',
                             profile.raypath))
            top = ax.get_ylim()[1]
            ax.set_ylim(-0.1, top)                             
            ax.set_title(title)                         
            pp.savefig()
            fig.clf()

time_series = cpx_spectra.PMR_unpol
for idx, spectrum in enumerate(time_series.spectra_list):
    spectrum.get_baseline()
    spectrum.get_peakfit()
    fig, ax = spectrum.plot_peakfit_and_baseline()
    fig.set_size_inches(6, 6)
    title = ''.join(('augite PMR dehydrated at 800$\degree\C\n ', 
                     '{:.1f}'.format(time_series.times_hours[idx]), ' hours'))
    ax.set_title(title)                         
    pp.savefig()
    fig.clf()
       
pp.close()        
        
    

