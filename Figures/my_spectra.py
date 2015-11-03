# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 14:20:16 2015

@author: Ferriss

My data: samples, FTIR spectra, and experiments
Class and function definitions are in pynams.py 

See Ferriss et al. 2015 American Mineralogist for description.
"""
import pynams.pynams as nams
import numpy as np
from uncertainties import ufloat
import matplotlib.pyplot as plt

plt.style.use('paper')
default_folder = 'C:\\Users\\Ferriss\\Documents\\FTIR\\'
default_savefolder = 'C:\\Users\\Ferriss\\Documents\\CpxPaper\\Figures\\'
style_profile_default = {'markeredgecolor' : 'blue', 'linestyle' : 'none', 
                         'marker' : 's', 'fillstyle' : 'none', 
                         'markersize' : 10, 'alpha' : 0.5}

# cpx paper initial concentrations for Kunlun, Jaipur, and PMR-53
# See KJP_3x3_comparison.py for details and peak-specific values
K_water_known = nams.ufloat(29, 4)
J_water_known = nams.ufloat(30, 2)
PMR_water_known = nams.ufloat(268, 8)

# peak distributions
K_init_peak_fraction = np.array([0.12, 0.32, 0.18, 0.19, 0.19])
J_init_peak_fraction = np.array([0.06, 0.09, 0.13, 0.34, 0.38])

# peak-specific water content - leaving out 3460 cm-1
K_water_peaks = K_water_known * K_init_peak_fraction
J_water_peaks = J_water_known * J_init_peak_fraction
#
### Saved preferred initial values (i), diffusivities (D), and errors (e) ####
### In order: 3645, 3617, 3540, 3443, 3350, BULK H ###
i_K3 = [12., 25., 23.3, 23., 10., 1.24]
D_K3 = [-13.9, -13.4, -13., -13.7, -13.5, -13.6]
e_K3 = [0.5, 0.3, 0.3, 0.4, 0.4, 0.3]

i_K91 = [15.18, 20., 23.3, 23., 10., 1.24]
D_K91 = [-13.5, -12.9, -12.4, -13.5, -13.25, -13.35]
e_K91 = [0.50, 0.40, 0.15, 0.40, 0.25, 0.25]

i_K154 = [15.18, 20., 23.3, 23., 10., 1.24]
D_K154 = [-13.60, -13.05, -12.6, -13.40, -13.25, -13.35]
e_K154 = [0.50, 0.45, 0.15, 0.30, 0.25, 0.3]

i_K5 = [16.12, 20., 23.53, 23., 10., 1.24]
D_K5 = [-13.30, -12.65, -12.10, -13.00, -12.75, -12.9]
e_K5 = [0.50, 0.25, 0.10, 0.20, 0.25, 0.3]

i_J = [2.95, 6.32, 16.71, 14.37, 30., 1.32]
D_J = [-10.00, -9.70, -9.60, -9.80, -9.95, -9.9]
e_J = [0.30, 0.30, 0.10, 0.20, 0.25, 0.4]

#my_ilist = [i_K3, i_K91, i_K154, i_J, i_K5]
K_dlist = [D_K3, D_K91, D_K154, D_K5]
K_elist = [e_K3, e_K91, e_K154,e_K5]

K_bulk_D = np.ones(4)
K_3645_D = np.ones(4)
K_3617_D = np.ones(4)
K_3540_D = np.ones(4)
K_3443_D = np.ones(4)
K_3350_D = np.ones(4)
K_bulk_e = np.ones(4)
K_3645_e = np.ones(4)
K_3617_e = np.ones(4)
K_3540_e = np.ones(4)
K_3443_e = np.ones(4)
K_3350_e = np.ones(4)

for idx in range(4):
    K_bulk_D[idx] = K_dlist[idx][-1]
    K_3645_D[idx] = K_dlist[idx][0]
    K_3617_D[idx] = K_dlist[idx][1]
    K_3540_D[idx] = K_dlist[idx][2]
    K_3443_D[idx] = K_dlist[idx][3]
    K_3350_D[idx] = K_dlist[idx][4]
    K_bulk_e[idx] = K_elist[idx][-1]
    K_3645_e[idx] = K_elist[idx][0]
    K_3617_e[idx] = K_elist[idx][1]
    K_3540_e[idx] = K_elist[idx][2]
    K_3443_e[idx] = K_elist[idx][3]
    K_3350_e[idx] = K_elist[idx][4]

# Kunlun diopside samples and spectra
# 
class Sample_Kunlun(nams.Sample):
    mineral_name = 'diopside'
    IGSN = 'IEFERKUN0'
    source = 'AMNH, Kunlun Mts China'
    initial_water = K_water_known

K3 = Sample_Kunlun()
K3.twoA_list = [1998, 1988, 2002, 1998, 1997]
K3.twoB_list = [1481, 1484, 1479, 1482, 1482]
K3.twoC_list = [1803, 1801, 1800, 1796, 1804]
K3.sample_thick_microns = nams.get_3thick(K3)

K4 = Sample_Kunlun()
K4.twoA_list = 7000
K4.twoB_list = [2185, 2190, 2188, 2185, 2188]
K4.twoC_list = [1546, 1551, 1536, 1548, 1548]
K4.sample_thick_microns = nams.get_3thick(K4)

K5 = Sample_Kunlun()
K5.twoA_list = 3450
K5.twoB_list = [1614, 1609, 1612, 1602, 1607]
K5.twoC_list = [1756, 1759, 1748, 1763, 1759]
K5.sample_thick_microns = nams.get_3thick(K5)

K5slice = Sample_Kunlun()
K5slice.twoA_list = [955, 955, 954, 957, 953]
K5slice.twoB_list = K5.twoB_list
K5slice.twoC_list = K5.twoC_list

K6 = Sample_Kunlun()
K6.twoA_list = [6912, 6913, 6917, 6913, 6917]
K6.twoB_list = [2731, 2739, 2741, 2723, 2705]
K6.twoC_list = [1524, 1511, 1500, 1488, 1517]
K6.sample_thick_microns = nams.get_3thick(K6)
            
           
class SpecKunlun_RaypathC(nams.Spectrum):
    base_low_wn = 3500

class SpecKunlun_RaypathA(nams.Spectrum):
    base_low_wn = 3500
    
class KunlunProfile(nams.Profile):
    pass
#    def __init__(self):
#        self.style_base = {'markeredgecolor' : 'b', 'linestyle' : 'none', 
#                      'marker' : '^', 'markerfacecolor' : 'b', 
#                      'fillstyle' : 'none', 'markersize' : 10}

#
#%% Kunlun diopside FTIR spectra
# K3: 
# K4: 904 C
# K5: 1000 C

#%% K3: Kunlun diopside heated at 700 C for 2 hrs, then 817 C

### K3 true initial ### 
profile_K3_trueInit_raypathA = KunlunProfile()
profile_K3_trueInit_raypathA.sample = K3
profile_K3_trueInit_raypathA.profile_name = 'K3 initial R || a*'
profile_K3_trueInit_raypathA.direction = 'a'
profile_K3_trueInit_raypathA.raypath = 'b'
profile_K3_trueInit_raypathA.initial_profile = profile_K3_trueInit_raypathA
leng = profile_K3_trueInit_raypathA.set_len()
profile_K3_trueInit_raypathA.fname_list = ['K3_sp13_Bda', 'K3_sp14_Cda']
profile_K3_trueInit_raypathA.positions_microns = [leng/2., leng/2.]

profile_K3_trueInit_raypathB = KunlunProfile()
profile_K3_trueInit_raypathB.sample = K3
profile_K3_trueInit_raypathB.profile_name = 'K3 initial R || b'
profile_K3_trueInit_raypathB.direction = 'b'
profile_K3_trueInit_raypathB.raypath = 'c'
profile_K3_trueInit_raypathB.initial_profile = profile_K3_trueInit_raypathB
leng = profile_K3_trueInit_raypathB.set_len()
profile_K3_trueInit_raypathB.fname_list = ['K3_sp1_Adb', 'K3_sp3_Cdb']
profile_K3_trueInit_raypathB.positions_microns = [leng/2., leng/2.]

profile_K3_trueInit_raypathC = KunlunProfile()
profile_K3_trueInit_raypathC.sample = K3
profile_K3_trueInit_raypathC.profile_name = 'K3 initial R || c'
profile_K3_trueInit_raypathC.direction = 'c'
profile_K3_trueInit_raypathC.raypath = 'a'
profile_K3_trueInit_raypathC.initial_profile = profile_K3_trueInit_raypathC
profile_K3_trueInit_raypathC.spectrum_class_name = SpecKunlun_RaypathC
leng = profile_K3_trueInit_raypathC.set_len()
profile_K3_trueInit_raypathC.fname_list = ['K3_sp11_Adc', 'K3_sp12_Bdc']
profile_K3_trueInit_raypathC.positions_microns = [leng/2., leng/2.]

K3wb_trueInit = nams.WholeBlock()
K3wb_trueInit.name = (
    'Initial measurements for Kunlun diopside K3 along 3 ray paths')
K3wb_trueInit.profiles = [profile_K3_trueInit_raypathA, 
                          profile_K3_trueInit_raypathB, 
                          profile_K3_trueInit_raypathC]
K3wb_trueInit.time_seconds = 10.
K3wb_trueInit.worksheetname = 'K3 initial'

### K3 heated at 696 C for 2 hours 7/26/12 "Initial" ###
profile_K3_init_A = KunlunProfile()
profile_K3_init_A.sample = K3
profile_K3_init_A.profile_name = 'K3 heated at 696 C for 2hr || a*'
profile_K3_init_A.initial_profile = profile_K3_trueInit_raypathB
#profile_K3_init_A.initial_profile = profile_K3_init_A
profile_K3_init_A.direction = 'a'
profile_K3_init_A.raypath = 'b'
leng = profile_K3_init_A.set_len()
profile_K3_init_A.fname_list = ['K3_cdb06']
profile_K3_init_A.positions_microns = [leng/2.]

profile_K3_init_B = KunlunProfile()
profile_K3_init_B.sample = K3
profile_K3_init_B.profile_name = 'K3 heated at 696 C for 2hr || b, OFF CENTER'
profile_K3_init_B.initial_profile = profile_K3_trueInit_raypathC
#profile_K3_init_B.initial_profile = profile_K3_init_B
profile_K3_init_B.direction = 'b'
profile_K3_init_B.raypath = 'c'
profile_K3_init_B.spectrum_class_name = SpecKunlun_RaypathC
leng = profile_K3_init_B.set_len()
profile_K3_init_B.fname_list = ['K3_bdc01', 'K3_bdc02', 'K3_bdc03', 
                                'K3_bdc04', 'K3_bdc05', 'K3_bdc06', 
                                'K3_bdc07', 'K3_bdc08', 'K3_bdc09', 
                                'K3_bdc10']
profile_K3_init_B.positions_microns = []
for k in range(10):
    profile_K3_init_B.positions_microns.append((k+1)*leng/11)

profile_K3_init_C = KunlunProfile()
profile_K3_init_C.sample = K3
profile_K3_init_C.profile_name = 'K3 heated at 696 C for 2hr || c'
profile_K3_init_C.initial_profile = profile_K3_trueInit_raypathB
#profile_K3_init_C.initial_profile = profile_K3_init_C
profile_K3_init_C.direction = 'c'
profile_K3_init_C.raypath = 'b'
leng = profile_K3_init_C.set_len()
profile_K3_init_C.fname_list = ['K3_cdb01', 'K3_cdb02', 'K3_cdb03', 
                                'K3_cdb04', 'K3_cdb05', 'K3_cdb06', 
                                'K3_cdb07', 'K3_cdb08', 'K3_cdb09', 
                                'K3_cdb10', 'K3_cdb11', 'K3_cdb12', 
                                'K3_cdb13']
profile_K3_init_C.positions_microns = np.linspace(50, leng-50, 
                                          len(profile_K3_init_C.fname_list))

K3wb_init = nams.WholeBlock()
K3wb_init.name = (
'"Initial" (pre-anneal at 696 C for 2 hrs) profiles for Kunlun diopside K3')
K3wb_init.profiles = [profile_K3_init_A, profile_K3_init_B, profile_K3_init_C]
K3wb_init.time_seconds = 2.*3600.
K3wb_init.worksheetname = 'K3 696C 2hr'

### K3 heated at 696 C for an additional 17hr 15m = 19hr 15m total
profile_K3_696C_19hr_A = KunlunProfile()
profile_K3_696C_19hr_A.sample = K3
profile_K3_696C_19hr_A.profile_name = 'K3 696 C for 19hr 15m || a'
profile_K3_696C_19hr_A.initial_profile = profile_K3_init_A
profile_K3_696C_19hr_A.direction = 'a'
profile_K3_696C_19hr_A.raypath = 'b'
leng = profile_K3_696C_19hr_A.set_len()
profile_K3_696C_19hr_A.fname_list = ['K3h_cdb02']
profile_K3_696C_19hr_A.positions_microns = [leng/2.]

profile_K3_696C_19hr_B = KunlunProfile()
profile_K3_696C_19hr_B.sample = K3
profile_K3_696C_19hr_B.profile_name = 'K3 696 C for 19hr 15m || b OFF-CENTER'
profile_K3_696C_19hr_B.initial_profile = profile_K3_init_B
profile_K3_696C_19hr_B.direction = 'b'
profile_K3_696C_19hr_B.raypath = 'c'
profile_K3_696C_19hr_B.spectrum_class_name = SpecKunlun_RaypathC
leng = profile_K3_696C_19hr_B.set_len()
profile_K3_696C_19hr_B.fname_list = ['K3h_bdc01', 'K3h_bdc02', 'K3h_bdc03',
                                     'K3h_bdc04', 'K3h_bdc05']
profile_K3_696C_19hr_B.positions_microns = [leng/5., 2.*leng/5., leng/2.,
                                            3.*leng/5., 4.*leng/5.]

profile_K3_696C_19hr_C = KunlunProfile()
profile_K3_696C_19hr_C.sample = K3
profile_K3_696C_19hr_C.profile_name = 'K3 696C for 19hr 15m || c'
profile_K3_696C_19hr_C.initial_profile = profile_K3_init_C
profile_K3_696C_19hr_C.direction = 'c'
profile_K3_696C_19hr_C.raypath = 'b'
leng = profile_K3_696C_19hr_C.set_len()
profile_K3_696C_19hr_C.fname_list = ['K3h_cdb01', 'K3h_cdb02', 'K3h_cdb03']
profile_K3_696C_19hr_C.positions_microns = [leng/4., leng/2., 3.*leng/4.]

K3wb_700C_19hr = nams.WholeBlock()
K3wb_700C_19hr.name = 'Kunlun diopside K3 at 696 C for 19hr 15 min'
K3wb_700C_19hr.profiles = [profile_K3_696C_19hr_A, profile_K3_696C_19hr_B,
                           profile_K3_696C_19hr_C]
K3wb_700C_19hr.time_seconds = (19.*3600.) + (15.*60.)
K3wb_700C_19hr.worksheetname = 'K3 696C 19hr'

### K3 heated at 696 C for an additional 16hr = 35hr 15m total
profile_K3_696C_35hr_A = KunlunProfile()
profile_K3_696C_35hr_A.sample = K3
profile_K3_696C_35hr_A.profile_name = 'K3 696C for 35hr 15m || a'
profile_K3_696C_35hr_A.initial_profile = profile_K3_init_A
profile_K3_696C_35hr_A.direction = 'a'
profile_K3_696C_35hr_A.raypath = 'b'
leng = profile_K3_696C_35hr_A.set_len()
profile_K3_696C_35hr_A.fname_list = ['K3h700_cdb01']
profile_K3_696C_35hr_A.positions_microns = [leng/2.]

profile_K3_696C_35hr_B = KunlunProfile()
profile_K3_696C_35hr_B.sample = K3
profile_K3_696C_35hr_B.profile_name = 'K3 696C for 35hr 15m || b, OFF CENTER'
profile_K3_696C_35hr_B.initial_profile = profile_K3_init_B
profile_K3_696C_35hr_B.direction = 'b'
profile_K3_696C_35hr_B.raypath = 'c'
profile_K3_696C_35hr_B.spectrum_class_name = SpecKunlun_RaypathC
leng = profile_K3_696C_35hr_B.set_len()
profile_K3_696C_35hr_B.fname_list = ['K3h700_bdc01', 'K3h700_bdc02', 
                                     'K3h700_bdc03']
profile_K3_696C_35hr_B.positions_microns = [leng/4., leng/2., 3.*leng/4.]

profile_K3_696C_35hr_C = KunlunProfile()
profile_K3_696C_35hr_C.sample = K3
profile_K3_696C_35hr_C.profile_name = 'K3 696C for 35hr 15m || c'
profile_K3_696C_35hr_C.initial_profile = profile_K3_init_C
profile_K3_696C_35hr_C.direction = 'c'
profile_K3_696C_35hr_C.raypath = 'b'
leng = profile_K3_696C_35hr_C.set_len()
profile_K3_696C_35hr_C.fname_list = ['K3h700_cdb01', 'K3h700_cdb02',
                                     'K3h700_cdb03']
profile_K3_696C_35hr_C.positions_microns = [leng/4., leng/2., 3.*leng/4.]

K3wb_700C_35hr = nams.WholeBlock()
K3wb_700C_35hr.name = 'Kunlun diopside K3 at 696 C for 35hr 15 min'
K3wb_700C_35hr.profiles = [profile_K3_696C_35hr_A, profile_K3_696C_35hr_B,
                           profile_K3_696C_35hr_C]
K3wb_700C_35hr.time_seconds = (19.*3600.) + (15.*60.) + (16.*3600.)
K3wb_700C_35hr.worksheetname = 'K3 696C 35hr'

### K3 heated at 796 C for 15 hr 40 m
profile_K3_800C_15hr40m_A = KunlunProfile()
profile_K3_800C_15hr40m_A.initial_profile = profile_K3_init_A
profile_K3_800C_15hr40m_A.sample = K3
profile_K3_800C_15hr40m_A.profile_name = (
                                'Kunlun diopside K3 at 796C for 15hr 40m || a')
profile_K3_800C_15hr40m_A.direction = 'a'
profile_K3_800C_15hr40m_A.raypath = 'b'
profile_K3_800C_15hr40m_A.fname_list = ['K3h800_cdb02']
leng = profile_K3_800C_15hr40m_A.set_len()
profile_K3_800C_15hr40m_A.positions_microns = [leng/2.]

profile_K3_800C_15hr40m_B = KunlunProfile()
profile_K3_800C_15hr40m_B.initial_profile = profile_K3_init_B
profile_K3_800C_15hr40m_B.sample = K3
profile_K3_800C_15hr40m_B.profile_name = (
                                'Kunlun diopside K3 at 796C for 15hr 40m || b')
profile_K3_800C_15hr40m_B.direction = 'b'
profile_K3_800C_15hr40m_B.raypath = 'c'
profile_K3_800C_15hr40m_B.spectrum_class_name = SpecKunlun_RaypathC
profile_K3_800C_15hr40m_B.fname_list = ['K3h800_bdc01', 'K3h800_bdc02', 
                                        'K3h800_bdcMID']
leng = profile_K3_800C_15hr40m_B.set_len()
profile_K3_800C_15hr40m_B.positions_microns = [50., 150., leng/2.]

profile_K3_800C_15hr40m_C = KunlunProfile()
profile_K3_800C_15hr40m_C.initial_profile = profile_K3_init_C
profile_K3_800C_15hr40m_C.sample = K3
profile_K3_800C_15hr40m_C.profile_name = (
                                'Kunlun diopside K3 at 796C for 15hr 40m || c')
profile_K3_800C_15hr40m_C.direction = 'c'
profile_K3_800C_15hr40m_C.raypath = 'b'
profile_K3_800C_15hr40m_C.fname_list = ['K3h800_cdb01', 'K3h800_cdb02', 
                                        'K3h800_cdbMID']
leng = profile_K3_800C_15hr40m_C.set_len()
profile_K3_800C_15hr40m_C.positions_microns = [50., 150., leng/2.]

K3wb_800C_15hr = nams.WholeBlock()
K3wb_800C_15hr.name = 'Kunlun diopside K3 at 796 C for 15hr 40min'
K3wb_800C_15hr.profiles = [profile_K3_800C_15hr40m_A, 
                           profile_K3_800C_15hr40m_B,
                           profile_K3_800C_15hr40m_C]
K3wb_800C_15hr.time_seconds = (15.*3600) + 40.
K3wb_800C_15hr.worksheetname = 'K3 800C 15hr'

### K3 heated at 817 C for 6 days
profile_K3_817C_6days_a = KunlunProfile()
profile_K3_817C_6days_a.time_seconds = 6.*24.*3600.
profile_K3_817C_6days_a.sample = K3
profile_K3_817C_6days_a.profile_name = 'K3 initial || a*'
profile_K3_817C_6days_a.short_name = 'K3g_adb'
profile_K3_817C_6days_a.initial_profile = profile_K3_init_A
profile_K3_817C_6days_a.direction = 'a' 
profile_K3_817C_6days_a.raypath = 'b'
leng = profile_K3_817C_6days_a.set_len()
profile_K3_817C_6days_a.fname_list = ['K3g_adb01', 'K3g_adb02', 'K3g_adb05', 
                                    'K3g_adb07', 'K3g_adb10', 'K3g_adb15', 
                                    'K3g_adb19', 'K3g_adb20']
profile_K3_817C_6days_a.positions_microns = [50., 150., 450., 
                                           650., 950., 1450., 
                                           leng-150., leng-50.]

profile_K3_817C_6days_b = KunlunProfile()
profile_K3_817C_6days_b.time_seconds = 6.*24.*3600
profile_K3_817C_6days_b.sample = K3
profile_K3_817C_6days_b.profile_name = 'K3 initial || b'
profile_K3_817C_6days_b.short_name = 'K3g_bdc'
profile_K3_817C_6days_b.initial_profile = profile_K3_init_B
profile_K3_817C_6days_b.direction = 'b' 
profile_K3_817C_6days_b.raypath = 'c'
profile_K3_817C_6days_b.spectrum_class_name = SpecKunlun_RaypathC
leng = profile_K3_817C_6days_b.set_len()
profile_K3_817C_6days_b.fname_list = ['K3g_bdc01', 'K3g_bdc02', 'K3g_bdc06', 
                                      'K3g_bdc10', 'K3g_bdc12']
profile_K3_817C_6days_b.positions_microns = [1150., 950., 550.,
                                             150., 50.]
# Confirm orientations...
#                                            50., 150., 550., 
#                                             950., 1150.]

profile_K3_817C_6days_c = KunlunProfile()
profile_K3_817C_6days_c.time_seconds = 6.*24.*3600
profile_K3_817C_6days_c.sample = K3
profile_K3_817C_6days_c.profile_name = 'K3 initial || c'
profile_K3_817C_6days_c.short_name = 'K3g_cdb'
profile_K3_817C_6days_c.initial_profile = profile_K3_init_C
profile_K3_817C_6days_c.direction = 'c' 
profile_K3_817C_6days_c.raypath = 'b'
leng = profile_K3_817C_6days_c.set_len()
profile_K3_817C_6days_c.fname_list = ['K3g_cdb01', 'K3g_cdb02', 'K3g_cdb03', 
                                      'K3g_cdb05', 'K3g_cdb06', 'K3g_cdb08', 
                                      'K3g_adb10', 'K3g_cdb12', 'K3g_cdb14', 
                                      'K3g_cdb17', 'K3g_cdb18']
# CHECK POSITION ORIENTATION HERE
profile_K3_817C_6days_c.positions_microns = [50., 150., 250., 
                                             450., 550., 750., 
                                             950., 1150., 1350., 
                                             leng-150., leng-50.]

K3wb_6days = nams.WholeBlock()
K3wb_6days.name = 'Kunlun diopside K3 heated at 817 C for 6 days'
K3wb_6days.profiles = [profile_K3_817C_6days_a, profile_K3_817C_6days_b, 
                       profile_K3_817C_6days_c]
K3wb_6days.time_seconds = 6.*24.*3600.
K3wb_6days.temperature_celsius = 817.
K3wb_6days.worksheetname = 'K3 817C 6days'

#%% K4: Kunlun diopside heated at 904 C

### Initial profiles ###
# initial // a (K4adcI)
profile_K4_init_A = KunlunProfile()
profile_K4_init_A.sample = K4
profile_K4_init_A.fname_list = ["K4_adcIL", "K4_dcImid", "K4_adcIr4"]
profile_K4_init_A.spectrum_class_name = SpecKunlun_RaypathC
profile_K4_init_A.direction = 'a'
profile_K4_init_A.raypath = 'c'
leng = profile_K4_init_A.set_len()
profile_K4_init_A.positions_microns = np.array([1000, leng/2.0, leng-1000])
profile_K4_init_A.profile_name = 'K4 initial || a'
profile_K4_init_A.initial_profile = profile_K4_init_A

# initial // b (K4bdcI)
profile_K4_init_B = KunlunProfile()
profile_K4_init_B.sample = K4
profile_K4_init_B.direction = 'b'
profile_K4_init_B.raypath = 'c'
leng = profile_K4_init_B.set_len()
profile_K4_init_B.spectrum_class_name = SpecKunlun_RaypathC
profile_K4_init_B.profile_name = 'K4 initial || b'
profile_K4_init_B.initial_profile = profile_K4_init_B
profile_K4_init_B.positions_microns = np.array([leng/2.0, leng/6.0, 
                                                leng*(2./3.), leng*(5./6.)])
profile_K4_init_B.fname_list = ['K4_dcImid', 'K4_dcIleft', 'K4_dcIr2', 
                                'K4_dcIright']

# initial // c (K4cdbI)
profile_K4_init_C = KunlunProfile()
profile_K4_init_C.sample = K4
profile_K4_init_C.profile_name = 'K4 initial || c'
profile_K4_init_C.direction = 'c'
profile_K4_init_C.raypath = 'b'
profile_K4_init_C.initial_profile = profile_K4_init_C
leng = profile_K4_init_C.set_len()
profile_K4_init_C.fname_list = ['K4_cdbIleft', 'K4_cdbImid', 'K4_cdbIright']
profile_K4_init_C.positions_microns = [leng/6.0, leng/2.0, leng*(5.0/6.0)]

K4wb_init = nams.WholeBlock()
K4wb_init.name = 'Initial profiles for Kunlun diopside K4'
K4wb_init.profiles = [profile_K4_init_A, profile_K4_init_B, profile_K4_init_C]
K4wb_init.time_seconds = 10.
K4wb_init.worksheetname = 'Kunlun 904C initial'

### K4 0-time experiement (quenched immediately by falling to bottom) ###
### Used as initial for other esperiments ###

# quench // a (K4adcQ)
profile_K4_quench_A = KunlunProfile()
profile_K4_quench_A.sample = K4
profile_K4_quench_A.profile_name = 'Kunlun 0-time experiment || a'
profile_K4_quench_A.direction = 'a'
profile_K4_quench_A.raypath = 'c'
profile_K4_quench_A.initial_profile = profile_K4_init_A
leng = profile_K4_quench_A.set_len()
profile_K4_quench_A.fname_list = ['K4q_adc05', 'K4q_bdcMID', 'K4q_adc65']
profile_K4_quench_A.positions_microns = [525., leng/2., 6525.]
profile_K4_quench_A.spectrum_class_name = SpecKunlun_RaypathC

profile_K4_quench_B = KunlunProfile()
profile_K4_quench_B.sample = K4
profile_K4_quench_B.profile_name = 'Kunlun 0-time experiment || b'
profile_K4_quench_B.direction = 'b'
profile_K4_quench_B.raypath = 'c'
profile_K4_quench_B.initial_profile = profile_K4_init_B
leng = profile_K4_quench_B.set_len()
profile_K4_quench_B.fname_list = ['K4q_bdc01', 'K4q_bdcMID', 'K4q_bdc02']
profile_K4_quench_B.positions_microns = [120., leng/2., leng-120.]
profile_K4_quench_B.spectrum_class_name = SpecKunlun_RaypathC

profile_K4_quench_C = KunlunProfile()
profile_K4_quench_C.sample = K4
profile_K4_quench_C.profile_name = 'Kunlun 0-time experiment || c'
profile_K4_quench_C.direction = 'c'
profile_K4_quench_C.raypath = 'b'
profile_K4_quench_C.initial_profile = profile_K4_init_C
leng = profile_K4_quench_C.set_len()
profile_K4_quench_C.fname_list = ['K4q_cdb01', 'K4q_cdbMID', 'K4q_cdb02',
                                  'K4q_cdb03']
profile_K4_quench_C.positions_microns = [100., leng/2., leng-110., leng-220.]

K4wb_quench = nams.WholeBlock()
K4wb_quench.name = '0-time (37m at 480C) profiles at 904 C for Kunlun diopside K4'
K4wb_quench.profiles = [profile_K4_quench_A, profile_K4_quench_B, 
                        profile_K4_quench_C]
K4wb_quench.time_seconds = 30.
K4wb_quench.worksheetname = 'Kunlun 904C 0-time'

### K4 heated 1 hour ###
# Kunlun K4adcH
profile_K4_904C_1hr_A = KunlunProfile()
profile_K4_904C_1hr_A.sample = K4
profile_K4_904C_1hr_A.direction = 'a'
profile_K4_904C_1hr_A.raypath = 'c'
profile_K4_904C_1hr_A.profile_name = 'Kunlun heated 1 hr at 904 C || a'
profile_K4_904C_1hr_A.fname_list = ['K4h_adc01i', 'K4h_adc02', 'K4h_adc15',
                                    'K4h_adc22', 'K4h_adc29', 'K4h_dcMID',
                                    'K4h_adc42', 'K4h_adc56', 'K4h_adc66',
                                    'K4h_adc68i', 'K4h_adc68e']
profile_K4_904C_1hr_A.spectrum_class_name = SpecKunlun_RaypathC
leng = profile_K4_904C_1hr_A.set_len()
profile_K4_904C_1hr_A.positions_microns = [100, 200, 1500, 2200, 2900, 
                                           leng/2.0, 4200, 5600, 6600, 
                                           leng-150, leng-100]
# Why are these not included? Did I lose the files?
del profile_K4_904C_1hr_A.positions_microns[-2:]
del profile_K4_904C_1hr_A.fname_list[-2:]
profile_K4_904C_1hr_A.initial_profile = profile_K4_quench_A

# Kunlun K4bdcH
profile_K4_904C_1hr_B = KunlunProfile()
profile_K4_904C_1hr_B.sample = K4
profile_K4_904C_1hr_B.direction = 'b'
profile_K4_904C_1hr_B.raypath = 'c'
leng = profile_K4_904C_1hr_B.set_len()
profile_K4_904C_1hr_B.profile_name = 'Kunlun heated 1 hr at 904 C || b'
profile_K4_904C_1hr_B.spectrum_class_name = SpecKunlun_RaypathC
profile_K4_904C_1hr_B.initial_profile = profile_K4_quench_B
profile_K4_904C_1hr_B.fname_list = ['K4h_bdc01e', 'K4h_bdc01i', 'K4h_bdc02', 
                                    'K4h_bdc04', 'K4h_bdc06', 'K4h_dcMID', 
                                    'K4h_bdc13', 'K4h_bdc15', 'K4h_bdc17', 
                                    'K4h_bdc20', 'K4h_bdc21i', 'K4h_bdc21e']
profile_K4_904C_1hr_B.positions_microns = [50., 100., 200., 
                                           400., 600., leng/2., 
                                           1300., 1500., 1700., 
                                           leng-200, leng-100., leng-50.]

# Kunlun K4cdbH
profile_K4_904C_1hr_C = KunlunProfile()
profile_K4_904C_1hr_C.sample = K4
profile_K4_904C_1hr_C.direction = 'c'
profile_K4_904C_1hr_C.raypath = 'b'
leng = profile_K4_904C_1hr_C.set_len()
profile_K4_904C_1hr_C.profile_name = 'Kunlun heated 1 hr at 904 C || c'
profile_K4_904C_1hr_C.initial_profile = profile_K4_quench_C
profile_K4_904C_1hr_C.fname_list = ['K4h_cdb01', 'K4h_cdb02', 'K4h_cdb03', 
                                    'K4h_cdb04', 'K4h_cdb08', 'K4h_cdb10', 
                                    'K4h_cdb12', 'K4h_cdb14', 'K4h_cdb15', 
                                    'K4h_cdb16']
profile_K4_904C_1hr_C.positions_microns = [115., 215., 315., 
                                           415., leng/2., 1015., 
                                           1215., 1415., 1515., 
                                           leng-115]

K4wb_1hr = nams.WholeBlock()
K4wb_1hr.name = 'Kunlun diopside K4 heated at 904 C for 1 hour'
K4wb_1hr.profiles = [profile_K4_904C_1hr_A, profile_K4_904C_1hr_B, 
                     profile_K4_904C_1hr_C]
K4wb_1hr.time_seconds = 3600.
K4wb_1hr.worksheetname = 'Kunlun 904C 1 hour'

### K4 heated 91 hours (K4P) ###
profile_K4_904C_91hr_A = KunlunProfile()
profile_K4_904C_91hr_A.time_seconds = 91.*3600.
profile_K4_904C_91hr_A.sample = K4
profile_K4_904C_91hr_A.profile_name = 'Kunlun K4 heated 904 C for 91 hr || a'
profile_K4_904C_91hr_A.short_name = 'K4adc_P'
profile_K4_904C_91hr_A.direction = 'a'
profile_K4_904C_91hr_A.raypath = 'c'
profile_K4_904C_91hr_A.spectrum_class_name = SpecKunlun_RaypathC
profile_K4_904C_91hr_A.initial_profile = profile_K4_quench_A
#profile_K4_904C_91hr_A.initial_profile = profile_K4_904C_1hr_A
leng = profile_K4_904C_91hr_A.set_len()
profile_K4_904C_91hr_A.fname_list = ['K4p_adc01', 'K4p_adc02', 'K4p_adc08', 
                                     'K4p_adc14', 'K4p_adc22', 'K4p_adc29', 
                                     'K4p_adc35', 'K4p_adc42', 'K4p_adc49', 
                                     'K4p_adc56', 'K4p_adc61', 'K4p_adc67']
profile_K4_904C_91hr_A.positions_microns = [100., 200., 800., 
                                            1400., 2200., 2900., 
                                            3500., 4200., 4900., 
                                            5600., 6100., 6700.]

profile_K4_904C_91hr_B = KunlunProfile()
profile_K4_904C_91hr_B.time_seconds = 91.*3600.
profile_K4_904C_91hr_B.sample = K4
profile_K4_904C_91hr_B.profile_name = 'Kunlun K4 heated 904 C for 91 hr || b'
profile_K4_904C_91hr_B.short_name = 'K4bdc_P'
profile_K4_904C_91hr_B.direction = 'b'
profile_K4_904C_91hr_B.raypath = 'c'
profile_K4_904C_91hr_B.spectrum_class_name = SpecKunlun_RaypathC
profile_K4_904C_91hr_B.initial_profile = profile_K4_quench_B
#profile_K4_904C_91hr_B.initial_profile = profile_K4_904C_1hr_B
leng = profile_K4_904C_91hr_B.set_len()
profile_K4_904C_91hr_B.fname_list = ['K4p_bdc01', 'K4p_bdc02', 'K4p_bdc04', 
                                     'K4p_bdc06', 'K4p_bdc09', 'K4p_adc35', 
                                     'K4p_bdc13', 'K4p_bdc15', 'K4p_bdc17', 
                                     'K4p_bdc19', 'K4p_bdc21']

profile_K4_904C_91hr_B.positions_microns = [100., 200., 400., 
                                            600., 900., 1200., 
                                            1400., 1500., 1700., 
                                            1900., 2100.]

profile_K4_904C_91hr_C = KunlunProfile()
profile_K4_904C_91hr_C.time_seconds = 91.*3600.
profile_K4_904C_91hr_C.sample = K4
profile_K4_904C_91hr_C.profile_name = 'Kunlun K4 heated 904 C for 91 hr || c'
profile_K4_904C_91hr_C.short_name = 'K4cdb_P'
profile_K4_904C_91hr_C.direction = 'c'
profile_K4_904C_91hr_C.raypath = 'b'
profile_K4_904C_91hr_C.initial_profile = profile_K4_quench_C
#profile_K4_904C_91hr_C.initial_profile = profile_K4_904C_1hr_C
leng = profile_K4_904C_91hr_C.set_len()
profile_K4_904C_91hr_C.fname_list = ['K4p_cdb01', 'K4p_cdb02', 'K4p_cdb04', 
                                     'K4p_cdb06', 'K4p_cdb09', 'K4p_cdb11', 
                                     'K4p_cdb13', 'K4p_cdb15', 'K4p_cdb16']
profile_K4_904C_91hr_C.positions_microns = [50., 150., 400., 
                                            600., 900., 1100., 
                                            1300., leng-150., leng-50.]

K4wb_91hr = nams.WholeBlock()
K4wb_91hr.name = 'Kunlun diopside K4 heated at 904 C for 91 hours'
K4wb_91hr.worksheetname = 'Kunlun 904C 91 hours'
K4wb_91hr.temperature_celsius = 904.
K4wb_91hr.profiles = [profile_K4_904C_91hr_A, 
                      profile_K4_904C_91hr_B, 
                      profile_K4_904C_91hr_C]
K4wb_91hr.time_seconds = 91.*3600.

### K4 heated 154 hours ###
# Kunlun K4adcF
profile_K4_904C_154hr_A = KunlunProfile()
profile_K4_904C_154hr_A.sample = K4
profile_K4_904C_154hr_A.direction = 'a'
profile_K4_904C_154hr_A.raypath = 'c'
profile_K4_904C_154hr_A.spectrum_class_name = SpecKunlun_RaypathC
profile_K4_904C_154hr_A.initial_profile = profile_K4_quench_A
profile_K4_904C_154hr_A.profile_name = (
                            'Kunlun heated 154 hours at 904C || a*')
profile_K4_904C_154hr_A.short_name = 'K4adcF'
profile_K4_904C_154hr_A.time_seconds = 154.*3600.
profile_K4_904C_154hr_A.fname_list = ['K4f_adc01', 'K4f_adc02', 'K4f_adc04',
                                      'K4f_adc06', 'K4f_adc08', 'K4f_adc11',
                                      'K4f_adc14', 'K4f_adc18', 'K4f_adc22',
                                      'K4f_adc25', 'K4f_adc29', 'K4f_adc35',
                                      'K4f_adc42', 'K4f_adc49', 'K4f_adc56',
                                      'K4f_adc61', 'K4f_adc67', 'K4f_adc68']
profile_K4_904C_154hr_A.positions_microns = [100, 200, 400, 
                                             600, 800, 1100,
                                             1400, 1800, 2200, 
                                             2500, 2900, 3500, 
                                             4200, 4900, 5600, 
                                             6100, 6700, 6800]

profile_K4_904C_154hr_B = KunlunProfile()
profile_K4_904C_154hr_B.sample = K4
profile_K4_904C_154hr_B.direction = 'b'
profile_K4_904C_154hr_B.raypath = 'c'
profile_K4_904C_154hr_B.time_seconds = 154.*3600.
profile_K4_904C_154hr_B.initial_profile = profile_K4_quench_B
profile_K4_904C_154hr_B.spectrum_class_name = SpecKunlun_RaypathC
profile_K4_904C_154hr_B.profile_name = (
                            'Kunlun heated 154 hours at 904C || b')
profile_K4_904C_154hr_B.short_name = 'K4bdcF'
leng = profile_K4_904C_154hr_B.set_len()
profile_K4_904C_154hr_B.fname_list = ['K4f_bdc01', 'K4f_bdc02', 'K4f_bdc03', 
                                      'K4f_bdc04', 'K4f_bdc05', 'K4f_bdc06', 
                                      'K4f_bdc07', 'K4f_bdc08', 'K4f_bdc09', 
                                      'K4f_bdc10', 'K4f_adc35', 'K4f_bdc15', 
                                      'K4f_bdc17', 'K4f_bdc19', 'K4f_bdc21', 
                                      'K4f_bdc22']
profile_K4_904C_154hr_B.positions_microns = [100., 200., 300., 
                                             400., 300., 600., 
                                             700., 800., 900., 
                                             1000., 1200., 1500., 
                                             1700., 1900., leng-150., 
                                             leng-50.]
                                             
profile_K4_904C_154hr_C = KunlunProfile()
profile_K4_904C_154hr_C.sample = K4
profile_K4_904C_154hr_C.direction = 'c'
profile_K4_904C_154hr_C.raypath = 'b'
profile_K4_904C_154hr_C.initial_profile = profile_K4_quench_C
profile_K4_904C_154hr_C.profile_name = (
                            'Kunlun heated 154 hours at 904C || c')
profile_K4_904C_154hr_C.short_name = 'K4cdbF'
profile_K4_904C_154hr_C.fname_list = ['K4f_cdb01', 'K4f_cdb02', 'K4f_cdb03',
                                      'K4f_cdb04', 'K4f_cdb05', 'K4f_cdb06',
                                      'K4f_cdb07', 'K4f_cdb09', 'K4f_cdb11',
                                      'K4f_cdb12', 'K4f_cdb13', 'K4f_cdb14',
                                      'K4f_cdb16']
leng = profile_K4_904C_154hr_C.set_len()
profile_K4_904C_154hr_C.positions_microns = [50, 120, 300, 400, 500, 600, 700,
                                             900, 1100, 1200, 1300, 1400, 
                                             leng-50]

K4wb_154hr = nams.WholeBlock()
K4wb_154hr.name = 'Kunlun diopside K4 heated at 904 C for 154 hours'
K4wb_154hr.worksheetname = 'Kunlun 904C 154 hours'
K4wb_154hr.temperature_celsius = 904.
K4wb_154hr.profiles = [profile_K4_904C_154hr_A, 
                      profile_K4_904C_154hr_B, 
                      profile_K4_904C_154hr_C]
K4wb_154hr.time_seconds = 154.*3600.

#%% K5 - bulk data published in whole-block paper - 1000C

### Initial profiles
# Kunlun K5adbI
profile_K5_initial_A = KunlunProfile()
profile_K5_initial_A.sample = K5
profile_K5_initial_A.direction = 'a'
profile_K5_initial_A.raypath = 'b'
profile_K5_initial_A.profile_name = 'Kunlun K5 initial || a'
profile_K5_initial_A.fname_list = ['K5i_adb01', 'K5i_dbMID', 'K5i_adb32']
leng = profile_K5_initial_A.set_len()
profile_K5_initial_A.positions_microns = [100., leng/2., leng-200]
profile_K5_initial_A.wb_initial_profile = profile_K5_initial_A

# Kunlun K5bdcI
profile_K5_initial_B = KunlunProfile()
profile_K5_initial_B.sample = K5
profile_K5_initial_B.direction = 'b'
profile_K5_initial_B.raypath = 'c'
profile_K5_initial_B.fname_list = ['K5i_bdc01', 'K5i_bdcMID', 'K5i_bdc16']
profile_K5_initial_B.profile_name = 'Kunlun K5 initial || b'
leng = profile_K5_initial_B.set_len()
profile_K5_initial_B.positions_microns = [50., leng/2., leng-50]
profile_K5_initial_B.spectrum_class_name = SpecKunlun_RaypathC

# Kunlun K5cdbI
profile_K5_initial_C = KunlunProfile()
profile_K5_initial_C.sample = K5
profile_K5_initial_C.direction = 'c'
profile_K5_initial_C.raypath = 'b'
profile_K5_initial_C.profile_name = 'Kunlun K5 initial || c'
profile_K5_initial_C.fname_list = ['K5i_cdb01', 'K5i_dbMID', 'K5i_cdb16']
leng = profile_K5_initial_C.set_len()
profile_K5_initial_C.positions_microns = [100., leng/2., 1600.]

K5wb_init = nams.WholeBlock()
K5wb_init.name = 'Kunlun diopside K5 initial'
K5wb_init.worksheetname = 'Kunlun 1000C initial'
K5wb_init.profiles = [profile_K5_initial_A,
                      profile_K5_initial_B, 
                      profile_K5_initial_C]
K5wb_init.time_seconds = 10.

### Final profiles
# Kunlun K5adb
profile_K5_1000C_75hr_A = KunlunProfile()
profile_K5_1000C_75hr_A.sample = K5
profile_K5_1000C_75hr_A.direction = 'a'
profile_K5_1000C_75hr_A.raypath = 'b'
profile_K5_1000C_75hr_A.time_seconds = 75.*3600.
profile_K5_1000C_75hr_A.profile_name = 'Kunlun heated at 1000C for 75hr || a*'
profile_K5_1000C_75hr_A.short_name = 'K5adb'
profile_K5_1000C_75hr_A.initial_profile = profile_K5_initial_A
profile_K5_1000C_75hr_A.diffusivity_log10m2s = -13.04
profile_K5_1000C_75hr_A.diff_error = 0.12
profile_K5_1000C_75hr_A.peak_diffusivities = [0., -13.27, -13.08, 0., 0., 0.]
profile_K5_1000C_75hr_A.peak_diff_error = [0., 0.29, 0.25, 0., 0., 0.]
profile_K5_1000C_75hr_A.fname_list = ['K5_adb01', 'K5_adb02', 'K5_adb03', 
                                      'K5_adb04', 'K5_adb05', 'K5_adb06',
                                      'K5_adb07', 'K5_adb09', 'K5_adb11',
                                      'K5_adb14', 'K5_cdb09', 'K5_adb21',
                                      'K5_adb24', 'K5_adb25', 'K5_adb26',
                                      'K5_adb27', 'K5_adb28', 'K5_adb29',
                                      'K5_adb30', 'K5_adb31']
profile_K5_1000C_75hr_A.positions_microns = [50., 150., 300., 
                                             400., 500., 600.,
                                             700., 900., 1100., 
                                             1400., 1800., 2100., 
                                             2400., 2500., 2600., 
                                             2700., 2800., 2900.,
                                             3000., 3100.]                                      

# Kunlun K5bdc
profile_K5_1000C_75hr_B = KunlunProfile()
profile_K5_1000C_75hr_B.sample = K5
profile_K5_1000C_75hr_B.direction = 'b'
profile_K5_1000C_75hr_B.raypath = 'c'
profile_K5_1000C_75hr_B.time_seconds = 75.*3600.
profile_K5_1000C_75hr_B.profile_name = 'Kunlun heated at 1000C for 75hr || b'
profile_K5_1000C_75hr_B.short_name = 'K5bdc'
profile_K5_1000C_75hr_B.initial_profile = profile_K5_initial_B
profile_K5_1000C_75hr_B.diffusivity_log10m2s = -13.41
profile_K5_1000C_75hr_B.diff_error = 0.08
profile_K5_1000C_75hr_B.peak_diffusivities = [0., -12.38, -12.07, 0., 0., 0.]
profile_K5_1000C_75hr_B.peak_diff_error = [0., 0.03, 0.02, 0., 0., 0.]
profile_K5_1000C_75hr_B.fname_list = ['K5_bdc01', 'K5_bdc02', 'K5_bdc03',
                                      'K5_bdc04', 'K5_bdc05', 'K5_bdc06',
                                      'K5_bdc07', 'K5_bdc08', 'K5_bdc09',
                                      'K5_bdc10', 'K5_bdc11', 'K5_bdc12',
                                      'K5_bdc13', 'K5_bdc14', 'K5_bdc15',
                                      'K5_bdc16']
leng = profile_K5_1000C_75hr_B.set_len()
profile_K5_1000C_75hr_B.positions_microns = [50., 150., 300., 
                                             400., 500., 600.,
                                             700., 800., 900., 
                                             1000., 1100., 1200., 
                                             1300., 1400., leng-150., 
                                            leng-50.]
profile_K5_1000C_75hr_B.spectrum_class_name = SpecKunlun_RaypathC

# Kunlun K5cdb
profile_K5_1000C_75hr_C = KunlunProfile()
profile_K5_1000C_75hr_C.sample = K5
profile_K5_1000C_75hr_C.direction = 'c'
profile_K5_1000C_75hr_C.raypath = 'b'
profile_K5_1000C_75hr_C.time_seconds = 75.*3600.
profile_K5_1000C_75hr_C.profile_name = 'Kunlun heated at 1000C for 75hr || c'
profile_K5_1000C_75hr_C.short_name = 'K5cdb'
profile_K5_1000C_75hr_C.diffusivity_log10m2s = -13.58
profile_K5_1000C_75hr_C.diff_error = 0.12
profile_K5_1000C_75hr_C.peak_diffusivities = [0., -13.36, -13.08, 0., 0., 0.]
profile_K5_1000C_75hr_C.peak_diff_error = [0., 0.18, 0.17, 0., 0., 0.]
profile_K5_1000C_75hr_C.initial_profile = profile_K5_initial_C
profile_K5_1000C_75hr_C.fname_list = ['K5_cdb01', 'K5_cdb02', 'K5_cdb03', 
                                      'K5_cdb04', 'K5_cdb05', 'K5_cdb06', 
                                      'K5_cdb07', 'K5_cdb08', 'K5_cdb09', 
                                      'K5_cdb10', 'K5_cdb11', 'K5_cdb12', 
                                      'K5_cdb13', 'K5_cdb14', 'K5_cdb15', 
                                      'K5_cdb16']
profile_K5_1000C_75hr_C.positions_microns = [50., 150., 300., 
                                             400., 500., 600., 
                                             700., 800., 900., 
                                             1000., 1100., 1200., 
                                             1300., 1400., 1500., 
                                             1600]

K5wb_75hr = nams.WholeBlock()
K5wb_75hr.name = 'Kunlun diopside K5 heated at 1000 C for 75 hours'
K5wb_75hr.worksheetname = 'Kunlun 1000C 75hr'
K5wb_75hr.temperature_celsius = 1000.
K5wb_75hr.profiles = [profile_K5_1000C_75hr_A,
                      profile_K5_1000C_75hr_B, 
                      profile_K5_1000C_75hr_C]
K5wb_75hr.time_seconds = 75.*3600.


#profile K5 SLICE K5sbda
#initial state Use K3_sp13_Bda.CSV with thickness transforms
profile_K5_slice_1000C_75hr_B = KunlunProfile()
profile_K5_slice_1000C_75hr_B.sample = K5slice
profile_K5_slice_1000C_75hr_B.direction = 'b'
profile_K5_slice_1000C_75hr_B.raypath = 'a'
profile_K5_slice_1000C_75hr_B.profile_name = 'Kunlun slice data || b'
#profile_K5_slice_1000C_75hr_B.spectrum_class_name = SpecKunlun_RaypathA
# Notes in matlab also list K5c_bda05 at 490 microns, _bda07 at 690,
# _bda09 at 790, bda11 at 1090, and bda12 at 1190
profile_K5_slice_1000C_75hr_B.fname_list = [
                        'K5c_bda01', 'K5c_bda02', 'K5c_bda03', 'K5c_bda04', 
                                     'K5c_bda06',              'K5c_cda09', 
                                     'K5c_bda10',              
                        'K5c_bda13', 'K5c_bda14', 'K5c_bda15', 'K5c_bda16']
profile_K5_slice_1000C_75hr_B.positions_microns = [ 50., 140., 290., 390., 
                                                         590.,       790., 
                                                         990.,
                                                   1290., 1390., 1490., 1538.]

# profile K5 SLICE K5ccda
profile_K5_slice_1000C_75hr_C = KunlunProfile()
profile_K5_slice_1000C_75hr_C.sample = K5slice
profile_K5_slice_1000C_75hr_C.direction = 'c'
profile_K5_slice_1000C_75hr_C.raypath = 'a'
profile_K5_slice_1000C_75hr_C.profile_name = 'Kunlun slice data || c'
profile_K5_slice_1000C_75hr_C.fname_list = [
                    'K5c_cda01', 'K5c_cda02', 'K5c_cda03', 'K5c_cda04', 
                    'K5c_cda05', 'K5c_cda06', 'K5c_cda07', 'K5c_cda08', 
                    'K5c_cda09', 'K5c_cda10', 'K5c_cda11', 'K5c_cda12', 
                    'K5c_cda13', 'K5c_cda14', 'K5c_cda15', 'K5c_cda16', 
                    'K5c_cda17']
leng = profile_K5_slice_1000C_75hr_C.set_len()
profile_K5_slice_1000C_75hr_C.positions_microns = [
                                                49.1, 149.1, 249.1, 349.1, 
                                                449.1, 549.1, 649.1, 749.1, 
                                                849.1, 949.1, 1049.1, 1149.1, 
                                                1249.1, 1349.1, 1446.7, 1546.7, 
                                                leng-50]

K5wb = nams.WholeBlock()
K5wb.name = 'Kunlun diopside heated at 1000 C for 75 hours'
K5wb.profiles = [profile_K5_1000C_75hr_A,
                 profile_K5_1000C_75hr_B,
                 profile_K5_1000C_75hr_C]
K5wb.style_base = style_profile_default
K5wb.time_seconds = 75.*3600
#K5wb.diffusivities_log10_m2s = [-13.2, -13.2, -13.2]


#%% K6: I have polarized measurements for

class SpecK6(nams.Spectrum):
    sample = K6


class SpecK6_db(SpecK6):
    raypath = 'b'
    thick_microns = K6.sample_thick_microns[1]


class SpecK6_dc(SpecK6):
    raypath = 'c'
    thick_microns = K6.sample_thick_microns[2]

K6_db_Ea = SpecK6_db()
K6_db_Ea.fname = 'K6_db_Ea'
K6_db_Ea.polar = 'E || a'

K6_db_Ea_2 = SpecK6_db()
K6_db_Ea_2.fname = 'K6_db_Ea_2'
K6_db_Ea_2.polar = 'E || a'

K6_dc_Ea = SpecK6_dc()
K6_dc_Ea.fname = 'K6_dc_Ea'
K6_dc_Ea.polar = 'E || a'

K6_dc_Ea_2 = SpecK6_dc()
K6_dc_Ea_2.fname = 'K6_dc_Ea_2'
K6_dc_Ea_2.polar = 'E || a'

K6_dc_Eb = SpecK6_dc()
K6_dc_Eb.fname = 'K6_dc_Eb'
K6_dc_Eb.polar = 'E || b'

K6_dc_Eb_2 = SpecK6_dc()
K6_dc_Eb_2.fname = 'K6_dc_Eb_2'
K6_dc_Eb_2.polar = 'E || b'

K6_db_Ec_2 = SpecK6_db()
K6_db_Ec_2.fname = 'K6_db_Ec_2'
K6_db_Ec_2.polar = 'E || c'

nams.make_filenames()

ave_K6_Ea = SpecK6()
ave_K6_Ea.other_name = 'Average of 4 initial K6 spectra, E || a'
ave_K6_Ea.make_average_spectra([K6_db_Ea, K6_db_Ea_2, K6_dc_Ea, K6_dc_Ea_2])
ave_K6_Ea.polar = 'E || a'
ave_K6_Ea.fname = 'ave_K6_Ea'
ave_K6_Ea.base_low_wn = 3200
ave_K6_Ea.base_high_wn = 3770
ave_K6_Ea.base_mid_wn = 3500
ave_K6_Ea.base_mid_yshift = 0.09
ave_K6_Ea.base_w_small = 0.035
ave_K6_Ea.base_w_large = 0.035

ave_K6_Eb = SpecK6()
ave_K6_Eb.other_name = 'Average of 2 initial K6 spectra, E || b'
ave_K6_Eb.make_average_spectra([K6_dc_Eb, K6_dc_Eb_2])
ave_K6_Eb.polar = 'E || b'
ave_K6_Eb.fname = 'ave_K6_Eb'
ave_K6_Eb.base_low_wn = 3200
ave_K6_Eb.base_high_wn = 3700
ave_K6_Eb.base_mid_wn = 3450
ave_K6_Eb.base_mid_yshift = 0.04
ave_K6_Eb.base_w_small = 0.035
ave_K6_Eb.base_w_large = 0.035

ave_K6_Ec = SpecK6()
ave_K6_Ec.other_name = 'Initial K6 spectrum with E || c'
ave_K6_Ec.make_average_spectra([K6_db_Ec_2])
ave_K6_Ec.polar = 'E || c'
ave_K6_Ec.fname = 'ave_K6_Ec'
ave_K6_Ec.base_low_wn = 3200
ave_K6_Ec.base_high_wn = 3750
ave_K6_Ec.base_mid_wn = 3450
ave_K6_Ec.base_mid_yshift = 0.05
ave_K6_Ec.base_w_small = 0.035
ave_K6_Ec.base_w_large = 0.035

#%%
class Sample_Jaipur(nams.Sample):
    mineral_name = 'diopside'
    source = 'Stephen Mackwell; north of Jaipur, India'
    initial_water = 26

J_CurveSnap = Sample_Jaipur()
J_CurveSnap.sample_thick_microns = 0.8e6
    
class SpecJ_CurveSnap(nams.Spectrum):
    sample = J_CurveSnap
    instrument = 'CurveSnap software - Woods et al. 2000'

# polarized spectra in three dimensions from Woods et al. 2000
J_Ea = SpecJ_CurveSnap()
J_Ea.other_name = 'J_Ea'
J_Ea.fname = 'J_Ea'
J_Ea.polar = 'a'
J_Ea.thick_microns = J_Ea.sample.sample_thick_microns
J_Ea.base_low_wn = 3200
J_Ea.base_high_wn = 3675
J_Ea.base_mid_wn = 3600
J_Ea.base_mid_yshift = 0.085
J_Ea.base_w_small = 0.01
J_Ea.base_w_large = 0.01

J_Eb = SpecJ_CurveSnap()
J_Eb.other_name = 'J_Eb'
J_Eb.fname = 'J_Eb'
J_Eb.polar = 'b'
J_Eb.thick_microns = J_Eb.sample.sample_thick_microns
J_Eb.base_low_wn = 3200
J_Eb.base_high_wn = 3700
J_Eb.base_mid_wn = 3400
J_Eb.base_mid_yshift = 0.05
J_Eb.base_w_small = 0.04
J_Eb.base_w_large = 0.04

J_Ec = SpecJ_CurveSnap()
J_Ec.other_name = 'J_Ec'
J_Ec.fname = 'J_Ec'
J_Ec.polar = 'c'
J_Ec.thick_microns = J_Ec.sample.sample_thick_microns
J_Ec.base_low_wn = 3200
J_Ec.base_high_wn = 3700
J_Ec.base_mid_wn = 3600
J_Ec.base_mid_yshift = 0.04
J_Ec.base_w_small = 0.025
J_Ec.base_w_large = 0.025

# J1: 
# a -> c
# b -> a
# c -> b
J1 = Sample_Jaipur()
J1.twoC_list = [2465, 2460, 2467, 2470, 2452] # originally twoA
J1.twoA_list = [4367, 4367, 4366, 4365, 4364] # originally twoB
J1.twoB_list = [3218, 3236, 3190, 3232, 3231] # originally twoC
J1.sample_thick_microns = nams.get_3thick(J1)

class ProfileJ1(nams.Profile):
    sample = J1
        
class SpecJ1(nams.Spectrum):
    sample = J1

class SpecJaipur_RaypathA(SpecJ1):
    base_high_wn = 3600

# J1 initial profiles
profile_J1_init_A = ProfileJ1()
profile_J1_init_A.fname_list = ['J1_dcI_t', 'J1_sp2', 'J1_dcI_bot']
profile_J1_init_A.direction = 'a'
profile_J1_init_A.raypath = 'b'
leng = profile_J1_init_A.set_len()
profile_J1_init_A.positions_microns = [leng*0.25, leng/2.0, leng*0.75]
profile_J1_init_A.profile_name = 'Jaipur diopside J1 initial || a*'
profile_J1_init_A.initial_profile = profile_J1_init_A

profile_J1_init_B = ProfileJ1()
profile_J1_init_B.fname_list = ['J1_sp3']
profile_J1_init_B.direction = 'b'
profile_J1_init_B.raypath = 'a'
leng = profile_J1_init_B.set_len()
profile_J1_init_B.positions_microns = [leng/2.0]
profile_J1_init_B.spectrum_class_name = SpecJaipur_RaypathA
profile_J1_init_B.profile_name = 'Jaipur diopside J1 initial || b'
profile_J1_init_B.initial_profile = profile_J1_init_B

profile_J1_init_C = ProfileJ1()
profile_J1_init_C.fname_list = ['J1_dcI_L', 'J1_sp2', 'J1_dcI_R']
profile_J1_init_C.direction = 'c'
profile_J1_init_C.raypath = 'b'
leng = profile_J1_init_C.set_len()
profile_J1_init_C.positions_microns = [leng*0.25, leng/2.0, leng*0.75]
profile_J1_init_C.profile_name = 'Jaipur diopside J1 initial || c'
profile_J1_init_C.initial_profile = profile_J1_init_C


# J1 after heating 30 minutes at 904 C       
profile_J1_904C_30m_A = ProfileJ1()
profile_J1_904C_30m_A.profile_name = 'Jaipur diopside at 904 C for 30 m || a*'
profile_J1_904C_30m_A.short_name = 'J1a'
profile_J1_904C_30m_A.direction = 'a'
profile_J1_904C_30m_A.raypath = 'b'
profile_J1_904C_30m_A.time_seconds = 30.*60.
leng = profile_J1_904C_30m_A.set_len()
profile_J1_904C_30m_A.initial_profile = profile_J1_init_A
profile_J1_904C_30m_A.fname_list = ['J1_bdc01', 'J1_bdc02', 'J1_bdc06', 
                                    'J1_bdc10', 'J1_bdc14', 'J1_bdc18', 
                                    'J1_dcMID', 'J1_bdc25', 'J1_bdc29', 
                                    'J1_bdc34', 'J1_bdc39', 'J1_bdc40', 
                                    'J1_bdc41', 'J1h_bdc05', 'J1h_bdc06', 
                                    'J1h_bdc07', 'J1h_bdc16', 'J1h_bdc17', 
                                    'J1h_bdc18', 'J1h_bdc19', 'J1h_bdc21', 
                                    'J1h_bdc22', 'J1h_bdc24', 'J1h_bdc26', 
                                    'J1h_bdc28',               'J1h_bdc33', 
                                    'J1h_bdc35', 'J1h_bdc37', 'J1h_bdc40', 
                                    'J1h_bdc41', 'J1h_bdc42']
profile_J1_904C_30m_A.positions_microns = [100., 200., 600., 
                                           1000., 1400., 1800., # crack at 1400
                                           # amphibole at 2900.
                                           leng/2., 2500., 2900., 
                                           3400., leng-400., leng-300., 
                                           leng-200., 500., 600., 
                                           700., 1600., 1700., 
                                           1800., 1900., 2100., 
                                           2200., 2400., 2600., 
                                           2800.,        3300., 
                                           3500., 3700., 4000., 
                                           4100., 4200.]

profile_J1_904C_30m_B = ProfileJ1()
profile_J1_904C_30m_B.profile_name = 'Jaipur diopside at 904 C for 30 m || b'
profile_J1_904C_30m_B.short_name = 'J1b'
profile_J1_904C_30m_B.direction = 'b'
profile_J1_904C_30m_B.raypath = 'a'
leng = profile_J1_904C_30m_B.set_len()
profile_J1_904C_30m_B.fname_list = ['J1_cdb01', 'J1_cdb02', 'J1_cdb05', 
                                    'J1_cdb08', 'J1_cdb12', 'J1_cdb15', 
                                    'J1_cdb18', 'J1_cdb22', 'J1_cdb25', 
                                    'J1_cdb28']
profile_J1_904C_30m_B.positions_microns = [100., 200., 500., 
                                           800., 1200., leng/2., 
                                           1800., 2200., leng-600., 
                                           2800.]
profile_J1_904C_30m_B.spectrum_class_name = SpecJaipur_RaypathA
profile_J1_904C_30m_B.initial_profile = profile_J1_init_B

# J1adc_orig.m
profile_J1_904C_30m_C = ProfileJ1()
profile_J1_904C_30m_C.profile_name = 'Jaipur diopside at 904 C for 30 m || c'
profile_J1_904C_30m_C.short_name = 'J1c'
profile_J1_904C_30m_C.direction = 'c'
profile_J1_904C_30m_C.raypath = 'b'
profile_J1_904C_30m_C.fname_list = ['J1_adc01', 'J1_adc02', 'J1_adc05', 
                                    'J1_adc10', 'J1_dcMID', 'J1_adc15', 
                                    'J1_adc20', 'J1_adc24', 'J1_adc25', 
                                    'J1h_adc09', 'J1h_adc10', 'J1h_adc11', 
                                    'J1h_adc12', 'J1h_adc13', 'J1h_adc15', 
                                    'J1h_adc16', 'J1h_adc17', 'J1h_adc18',
                                    'J1h_adc19', 'J1h_adc20', 'J1h_adc21', 
                                    'J1h_adc22', 'J1h_adc23', 'J1h_adc24', 
                                    'J1h_adc25']                                  
leng = profile_J1_904C_30m_C.set_len()
profile_J1_904C_30m_C.positions_microns = [75., 75.+100., 75+400., 
                                           75.+900., leng/2., (leng/2.) + 200., 
                                            (leng/2)+700., leng-175., leng-75., 
                                            900., 1000., 1100., 
                                            1200., 1300., 1500., 
                                            1600., 1700., 1800., 
                                            1900., 2000., 2100., 
                                            2200., 2300., 2400., 
                                            leng-50]
profile_J1_904C_30m_C.initial_profile = profile_J1_init_C


nams.make_filenames()

J1wb_initial = nams.WholeBlock()
J1wb_initial.name = 'Jaipur diopside initial'
J1wb_initial.time_seconds = 10.
J1wb_initial.worksheetname = 'Jaipur 904C initial'
J1wb_initial.profiles = [profile_J1_init_A, 
                         profile_J1_init_B, 
                         profile_J1_init_C]

J1wb = nams.WholeBlock()
J1wb.name = 'Jaipur diopside heated at 904 C for 30m'
J1wb.time_seconds = 30.*60.
J1wb.worksheetname = 'Jaipur 904C 30m'
J1wb.temperature_celsius = 904.
J1wb.profiles = [profile_J1_904C_30m_A, 
                 profile_J1_904C_30m_B, 
                 profile_J1_904C_30m_C]


class PMR(nams.Sample):
    initial_water = 268
    mineral_name = 'augite'
    source = 'David Bell'

PMR_CurveSnap = PMR()
PMR_CurveSnap.sample_thick_microns = 1e4

class SpecP_CurveSnap(nams.Spectrum):
    sample = PMR_CurveSnap
    instrument = 'CurveSnap software - Bell et al. 1995(?)'

PMR_Ea = SpecP_CurveSnap()
PMR_Ea.other_name = 'PMR_Ea'
PMR_Ea.fname = 'PMR_Ea'
PMR_Ea.polar = 'E || alpha'
PMR_Ea.thick_microns = PMR_Ea.sample.sample_thick_microns
PMR_Ea.base_low_wn = 3200
PMR_Ea.base_high_wn = 3700
PMR_Ea.base_mid_wn = 3500
PMR_Ea.base_mid_yshift = 1.65
PMR_Ea.base_w_small = 0.3
PMR_Ea.base_w_large = 0.3

PMR_Eb = SpecP_CurveSnap()
PMR_Eb.other_name = 'PMR_Eb'
PMR_Eb.fname = 'PMR_Eb'
PMR_Eb.polar = 'E || beta'
PMR_Eb.thick_microns = PMR_Eb.sample.sample_thick_microns
PMR_Eb.base_low_wn = 3200
PMR_Eb.base_high_wn = 3750
PMR_Eb.base_mid_wn = 3400
PMR_Eb.base_mid_yshift = 0.3
PMR_Eb.base_w_small = 0.3
PMR_Eb.base_w_large = 0.35

PMR_Ec = SpecP_CurveSnap()
PMR_Ec.other_name = 'PMR_Ec'
PMR_Ec.fname = 'PMR_Ec'
PMR_Ec.polar = 'E || gamma'
PMR_Ec.thick_microns = PMR_Ec.sample.sample_thick_microns
PMR_Ec.base_low_wn = 3200
PMR_Ec.base_high_wn = 3700
PMR_Ec.base_mid_wn = 3400
PMR_Ec.base_mid_yshift = 0.3
PMR_Ec.base_w_small = 0.3
PMR_Ec.base_w_large = 0.35


PULI_PMR53_alpha = nams.Spectrum()
PULI_PMR53_alpha.fname = "PULI_PMR53_alpha"
PULI_PMR53_alpha.thick_microns = 1e3

PULI_PMR53_alpha_baseline = nams.Spectrum()
PULI_PMR53_alpha_baseline.fname = 'PULI_PMR53_alpha_baseline'
PULI_PMR53_alpha_baseline.thick_microns = 1e3

PULI_PMR53_beta = nams.Spectrum()
PULI_PMR53_beta.fname = "PULI_PMR53_beta"
PULI_PMR53_beta.thick_microns = 1e3

PULI_PMR53_beta_baseline = nams.Spectrum()
PULI_PMR53_beta_baseline.fname = 'PULI_PMR53_beta_baseline'
PULI_PMR53_beta_baseline.thick_microns = 1e3

PULI_PMR53_gamma = nams.Spectrum()
PULI_PMR53_gamma.fname = "PULI_PMR53_gamma"
PULI_PMR53_gamma.thick_microns = 1e3

PULI_PMR53_gamma_baseline = nams.Spectrum()
PULI_PMR53_gamma_baseline.fname = 'PULI_PMR53_gamma_baseline'
PULI_PMR53_gamma_baseline.thick_microns = 1e3

PMR_alpha_both = nams.Profile()
PMR_beta_both = nams.Profile()
PMR_gamma_both = nams.Profile()
PMR_alpha_both.spectra_list = [PULI_PMR53_alpha, PULI_PMR53_alpha_baseline]
PMR_beta_both.spectra_list = [PULI_PMR53_beta, PULI_PMR53_beta_baseline]
PMR_gamma_both.spectra_list = [PULI_PMR53_gamma, PULI_PMR53_gamma_baseline]

nams.make_filenames()

###### PMR-53 augite that I dehydrated ######################################
PMR_dehydrated = nams.Sample()
PMR_dehydrated.sample_thick_microns = np.mean([885., 867., 873., 878., 879.])

class PMR_spectra(nams.Spectrum):
    sample = PMR_dehydrated
    thick_microns = PMR_dehydrated.sample_thick_microns

PMR_0 = PMR_spectra()
PMR_0.fname = 'P_0_unpol'

PMR_6 = PMR_spectra()
PMR_6.fname = 'P_6_unpol'

#############################################################################

ANU = nams.Sample()
ANU.mineral_name = 'clinopyroxene'
ANU.sample_thick_microns = 1235


class SpecANU(nams.Spectrum):
    sample = ANU
    thick_microns = ANU.sample_thick_microns
    polar = 'unpolarized'
    instrument = 'ANU (Penny King)'

ANU_pt57 = SpecANU()
ANU_pt57.fname = 'ANU-pt57'
ANU_pt57.other_name = 'CIP98-62-5'

nams.make_filenames()

#%% Groups of whole-blocks with the same diffusivties
heights_instead = False

D_PMR = nams.diffusion.Diffusivities()
D_PMR.description = 'augite PMR-53 bulk H'
D_PMR.celsius_all = [800.]
D_PMR.logD_unoriented = [-11.]
D_PMR.basestyle = {'color' : 'orange', 'marker' : 'D', 
                   'markersize' :  14, 'markeredgewidth' : 0,
                   'linestyle' : 'none', 'alpha' : 1.}



class K(nams.diffusion.Diffusivities):
    def __init__(self, diffusivities, description, color, marker, 
                 markersize=12, mew=1., alpha=1., error=[]):
        self.logDx = diffusivities
        self.logDy = diffusivities
        self.logDz = diffusivities
        self.logDx_error = error
        self.celsius_all = [816., 904., 904., 1000.]
        self.description = description
        self.basestyle = self.basestyle.copy()
        self.basestyle['color'] = color
        self.basestyle['marker'] = marker
        self.basestyle['markersize'] = markersize
        self.basestyle['mew'] = mew
        self.basestyle['alpha'] = alpha

#D_Kunlun_bulk = K(K_bulk_D, 'Bulk H', 'black', 'o', 
#           alpha=0.5,markersize=10, 
#           error=K_bulk_e)
#D_Kunlun_bulkH.basestyle = {'color' : 'crimson', 'marker' : 's', 
#                            'markersize' :  10, 'markeredgewidth' : 0,
#                            'linestyle' : 'none', 'alpha' : 0.5}
                            
K_bulk = K(K_bulk_D, 'Bulk H', 'black', 'o', 
           alpha=0.5,markersize=10, 
           error=K_bulk_e)
K_3645 = K(K_3645_D, '3645 cm$^{-1}$', 
           'red', 'x', mew=1.5, markersize=10, 
           error=K_3645_e)
K_3617 = K(K_3617_D, '3617 cm$^{-1}$', 
           'orange', '+', mew=2, markersize=10, 
           error=K_3617_e)
K_3540 = K(K_3540_D, '3540 cm$^{-1}$', 
           'teal', 's', markersize=8, 
           error=K_3540_e)
K_3443 = K(K_3443_D, '3443 cm$^{-1}$', 
           'blue', '3', mew=2, markersize=10, 
           error=K_3443_e)
K_3350 = K(K_3350_D, '3350 cm$^{-1}$', 
           'indigo', 'd', markersize=8, alpha=0.5, 
           error=K_3350_e)

#D_Kunlun_peak0 = nams.diffusion.Diffusivities()
#D_Kunlun_peak0.description = 'Kunlun 3645 cm$^{-1}$'
#D_Kunlun_peak0.wholeblocks = [K3wb_6days, K4wb_91hr, K4wb_154hr, K5wb_75hr]
##D_Kunlun_peak0.wholeblocks = [K4wb_154hr]
#D_Kunlun_peak0.basestyle = {'color' : 'peachpuff', 'marker' : 'o', 
#                            'markersize' :  10, 'markeredgewidth' : 0.5,
#                            'linestyle' : 'none', 'alpha' : 1,}
#D_Kunlun_peak0.get_from_wholeblock(peak_idx=0,
#                                   heights_instead=heights_instead)
#
#D_Kunlun_peak1 = nams.diffusion.Diffusivities()
#D_Kunlun_peak1.description = 'Kunlun 3617 cm$^{-1}$'
#D_Kunlun_peak1.wholeblocks = [K3wb_6days, K4wb_91hr, K4wb_154hr, K5wb_75hr]
#D_Kunlun_peak1.basestyle = {'color' : 'lightsalmon', 'marker' : '^', 
#                            'markersize' : 10, 'markeredgewidth' : 1.,
#                            'linestyle' : 'none', 'alpha' : 0.5,}
#D_Kunlun_peak1.get_from_wholeblock(peak_idx=1,
#                                   heights_instead=heights_instead)
#
#D_Kunlun_peak2 = nams.diffusion.Diffusivities()
#D_Kunlun_peak2.description = 'Kunlun 3540 cm$^{-1}$'
#D_Kunlun_peak2.wholeblocks = [K3wb_6days, K4wb_91hr, K4wb_154hr, K5wb_75hr]
#D_Kunlun_peak2.basestyle = {'color' : 'magenta', 'marker' : 'd', 
#                            'markersize' :  10, 'markeredgewidth' : 0,
#                            'linestyle' : 'none', 'alpha' : 1}
#D_Kunlun_peak2.get_from_wholeblock(peak_idx=2,
#                                   heights_instead=heights_instead)
#
#D_Kunlun_peak4 = nams.diffusion.Diffusivities()
#D_Kunlun_peak4.description = 'Kunlun 3443 cm$^{-1}$'
#D_Kunlun_peak4.wholeblocks = [K3wb_6days, K4wb_91hr, K4wb_154hr, K5wb_75hr]
#D_Kunlun_peak4.basestyle = {'color' : 'mistyrose', 'marker' : '*', 
#                            'markersize' :  17, 'markeredgewidth' : 1,
#                            'markeredgecolor' : 'black',
#                            'linestyle' : 'none', 'alpha' : 0.5}
#D_Kunlun_peak4.get_from_wholeblock(peak_idx=4,
#                                   heights_instead=heights_instead)
#
#D_Kunlun_peak5 = nams.diffusion.Diffusivities()
#D_Kunlun_peak5.description = 'Kunlun 3355 cm$^{-1}$'
#D_Kunlun_peak5.wholeblocks = [K3wb_6days, K4wb_154hr, K5wb_75hr]
##D_Kunlun_peak5.wholeblocks = [K4wb_154hr]
#D_Kunlun_peak5.basestyle = {'color' : 'darksalmon', 'marker' : 'p', 
#                            'markersize' :  12, 'markeredgewidth' : 1,
#                            'markeredgecolor' : 'black',
#                            'linestyle' : 'none', 'alpha' : 0.5}
#D_Kunlun_peak5.get_from_wholeblock(peak_idx=5,
#                                   heights_instead=heights_instead)
#
#### Jaipur
D_Jaipur_bulkH = nams.diffusion.Diffusivities()
D_Jaipur_bulkH.description = 'Jaipur bulk H (this work)'
D_Jaipur_bulkH.wholeblocks = [J1wb]
D_Jaipur_bulkH.basestyle = {'color' : 'k', 'marker' : 's', 
                            'markerfacecolor' : 'turquoise',
                            'markersize' :  10, 
                            'linestyle' : 'none', 'alpha' : 1.,
                            'markeredgewidth' : 2}
D_Jaipur_bulkH.get_from_wholeblock(wholeblock=True, 
                                   heights_instead=heights_instead)
D_Jaipur_bulkH.logDx = [-9.9]
D_Jaipur_bulkH.logDy = [-10.9]
D_Jaipur_bulkH.logDz = [-9.9]

D_Jaipur_peak0 = nams.diffusion.Diffusivities()
D_Jaipur_peak0.description = 'Jaipur 3645 cm-1'
D_Jaipur_peak0.wholeblocks = [J1wb]
D_Jaipur_peak0.basestyle = {'color' : 'purple', 'marker' : 'o', 
                            'markersize' :  10, 'linewidth' : 1,
                            'linestyle' : 'none', 'alpha' : 0.5,
                            'markeredgewidth' : 2}
D_Jaipur_peak0.get_from_wholeblock(peak_idx=0, wholeblock=True, 
                                   heights_instead=heights_instead)

D_Jaipur_peak1 = nams.diffusion.Diffusivities()
D_Jaipur_peak1.description = 'Jaipur 3617 cm-1'
D_Jaipur_peak1.wholeblocks = [J1wb]
D_Jaipur_peak1.basestyle = {'color' : 'blue', 'marker' : '^', 
                            'markersize' : 10, 
                            'linestyle' : 'none', 'alpha' : 0.5,
                            'markeredgewidth' : 2}
D_Jaipur_peak1.get_from_wholeblock(peak_idx=1,
                                   heights_instead=heights_instead)

D_Jaipur_peak2 = nams.diffusion.Diffusivities()
D_Jaipur_peak2.description = 'Jaipur 3540 cm-1'
D_Jaipur_peak2.wholeblocks = [J1wb]
D_Jaipur_peak2.basestyle = {'color' : 'indigo', 'marker' : 'd', 
                            'markersize' :  10, 'markeredgewidth' : 2,
                            'linestyle' : 'none', 'alpha' : 0.5,}
D_Jaipur_peak2.get_from_wholeblock(peak_idx=2,
                                   heights_instead=heights_instead)

D_Jaipur_peak3 = nams.diffusion.Diffusivities()
D_Jaipur_peak3.description = 'Jaipur 3460 cm-1'
D_Jaipur_peak3.wholeblocks = [J1wb]
D_Jaipur_peak3.basestyle = {'color' : 'blue', 'marker' : '*', 
                            'markersize' :  12, 'markeredgewidth' : 2,
                            'markeredgecolor' : 'black',
                            'linestyle' : 'none', 'alpha' : 0.5,}
D_Jaipur_peak3.get_from_wholeblock(peak_idx=3,
                                   heights_instead=heights_instead)

D_Jaipur_peak4 = nams.diffusion.Diffusivities()
D_Jaipur_peak4.description = 'Jaipur 3443 cm-1'
D_Jaipur_peak4.wholeblocks = [J1wb]
D_Jaipur_peak4.basestyle = {'color' : 'teal', 'marker' : '*', 
                            'markersize' :  12, 'markeredgewidth' : 1,
                            'markeredgecolor' : 'black',
                            'linestyle' : 'none', 'alpha' : 0.5,}
D_Jaipur_peak4.get_from_wholeblock(peak_idx=4,
                                   heights_instead=heights_instead)

D_Jaipur_peak5 = nams.diffusion.Diffusivities()
D_Jaipur_peak5.description = 'Jaipur 3355 cm-1'
D_Jaipur_peak5.wholeblocks = [J1wb]
D_Jaipur_peak5.basestyle = {'color' : 'green', 'marker' : 'p', 
                            'markersize' :  10, 'markeredgewidth' : 2,
                            'markeredgecolor' : 'black',
                            'linestyle' : 'none', 'alpha' : 0.5,}
D_Jaipur_peak5.get_from_wholeblock(peak_idx=5,
                                   heights_instead=heights_instead)

#H_diopside_difference = nams.diffusion.Diffusivities()
#H_diopside_difference.description = 'bulk Jaipur - bulk Kunlun'
#H_diopside_difference.celsius_x = D_Jaipur_bulkH.celsius_all * 2
#a = 10.**(np.array(D_Jaipur_bulkH.logDx*2))
#b = 10.**(np.array(D_Kunlun_bulkH.logDx[1:3]))
#H_diopside_difference.logDx
