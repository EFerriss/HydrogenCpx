# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 11:18:32 2016

@author: Ferriss

Convert .CSV files of polarized FTIR measurements of untreated Kunlun and 
Jaipur diopside to tab-delimited txt file appropriate for upload to the 
PULI database, http://puli.mfgi.hu
First column: wavenumber
Second column: absorbances
Normalized to cm

Wavenumber range 2700-4000 cm-1 following the example on the PULI website.

"""

import cpx_spectra as sp
import numpy as np
import csv
#import matplotlib.pyplot as plt
import pynams.pynams 

reload(pynams)
reload(sp)

Kunlun_list = [sp.ave_K6_Ea, sp.ave_K6_Eb, sp.ave_K6_Ec] 
J2_list = [sp.J2_Ea, sp.J2_Eb, sp.J2_Ec] # My Jaipur diopside sample J2

spec_list = J2_list + Kunlun_list
name_list = ['Jaipur_diopside_Ea', 'Jaipur_diopside_Eb', 'Jaipur_diopside_Ec',
             'Kunlun_diopside_Ea', 'Kunlun_diopside_Eb', 'Kunlun_diopside_Ec', ]

wn_high_puli = 4000.
wn_low_puli = 2500.

for idx, spec in enumerate(spec_list):
    spec.divide_by_thickness()
    idx_hi = (np.abs(spec.wn_full-wn_high_puli)).argmin()
    idx_lo = (np.abs(spec.wn_full-wn_low_puli)).argmin()

    wn_puli = spec.wn_full[idx_lo:idx_hi]
    abs_puli = spec.abs_full_cm[idx_lo:idx_hi]
    
    fig, ax = spec.plot_spectrum(wn_xlim_left=wn_high_puli,
                                 wn_xlim_right=wn_low_puli)
    
    filename = 'Ferriss_' + name_list[idx] + '_puli.txt'
    
    with open(filename, 'wb') as pulifile:
        spamwriter = csv.writer(pulifile, dialect='excel-tab')
        for idx_wn, wn in enumerate(wn_puli):
            spamwriter.writerow([wn, abs_puli[idx_wn]])
