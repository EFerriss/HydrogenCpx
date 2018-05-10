# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 14:20:16 2015

@author: Ferriss

Position data, sample data, FTIR spectra profiles and groups of profiles, 
and estimated diffusivities for hydrogen in clinopyroxene. 
Ferriss et al. (2016) Contributions to Mineralogy and Petrology 
https://link.springer.com/article/10.1007/s00410-016-1262-8

Note this script uses pynams v0.1.0, which has since been updated for clarity.
"""
from pynams import pynams as nams
from pynams import diffusion
import numpy as np
reload(nams)

# Set folder on your computer where to find .CSV FTIR spectra files 
folder = 'C:\\Users\\Ferriss\\Documents\\Code\\Python\\HydrogenCpx\\FTIR_CSV_FILES\\'

## SAMPLES
# IGSN: International Geo Sample Number; http://www.geosamples.org/aboutigsn

# Kunlun diopside
K3 = nams.Sample(IGSN='IEFERKUN3')
K4 = nams.Sample(IGSN='IEFERKUN4')
K5 = nams.Sample(IGSN='IEFERKUN5')
K5slice = nams.Sample(IGSN='IEFERK5SL')
K6 = nams.Sample(IGSN='IEFERKUN6')

# Jaipur diopside
J1 = nams.Sample(IGSN='IEFERJAI1')
J2 = nams.Sample(IGSN='IEFERJAI2')

# Augite PMR-53
PMR = nams.Sample(IGSN='IEFERPMR1')

#### SAMPLE THICKNESSES #####
K3.twoA_list = [1998, 1988, 2002, 1998, 1997]
K3.twoB_list = [1481, 1484, 1479, 1482, 1482]
K3.twoC_list = [1803, 1801, 1800, 1796, 1804]
K3.sample_thick_microns = nams.get_3thick(K3)

K4.twoA_list = 7000
K4.twoB_list = [2185, 2190, 2188, 2185, 2188]
K4.twoC_list = [1546, 1551, 1536, 1548, 1548]
K4.sample_thick_microns = nams.get_3thick(K4)

K5.twoA_list = 3450
K5.twoB_list = [1614, 1609, 1612, 1602, 1607]
K5.twoC_list = [1756, 1759, 1748, 1763, 1759]
K5.sample_thick_microns = nams.get_3thick(K5)

K5slice.twoA_list = [955, 955, 954, 957, 953]
K5slice.twoB_list = K5.twoB_list
K5slice.twoC_list = K5.twoC_list

K6.twoA_list = [6912, 6913, 6917, 6913, 6917]
K6.twoB_list = [2731, 2739, 2741, 2723, 2705]
K6.twoC_list = [1524, 1511, 1500, 1488, 1517]
K6.sample_thick_microns = nams.get_3thick(K6)

# Jaipur diopside J1 was originally mis-oriented, so the notes and filenames 
# all have somewhat misleading filenames.
# a -> c
# b -> a
# c -> b
J1.twoC_list = [2465, 2460, 2467, 2470, 2452] # originally twoA
J1.twoA_list = [4367, 4367, 4366, 4365, 4364] # originally twoB
J1.twoB_list = [3218, 3236, 3190, 3232, 3231] # originally twoC
J1.sample_thick_microns = nams.get_3thick(J1)

J2.twoA_list = [5129., 5125., 5117., 5102., 5136.]
J2.twoB_list = [3398., 3381., 3397., 3364., 3399.]

PMR.thickness_microns = np.mean([885., 867., 873., 878., 879.])

#%% K3 Kunlun diopside initial spectra
### unpolarized down three different ray paths


leng = np.mean(K3.twoA_list)
profile_K3_trueInit_raypathA = nams.Profile(sample=K3, folder=folder,
                                            profile_name = 'K3 initial R || a*',
                                            direction = 'a', raypath = 'b',
                                            fname_list = ['K3_sp13_Bda', 
                                                          'K3_sp14_Cda'],
                                            positions_microns = [leng/2., 
                                                                 leng/2.])

leng = np.mean(K3.twoB_list)
profile_K3_trueInit_raypathB = nams.Profile(sample = K3, folder=folder,
                                            profile_name = 'K3 initial R || b',
                                            direction = 'b', raypath = 'c',
                                            fname_list = ['K3_sp1_Adb', 
                                                          'K3_sp3_Cdb'],
                                            positions_microns = [leng/2., 
                                                                 leng/2.])

leng = np.mean(K3.twoC_list)
profile_K3_trueInit_raypathC = nams.Profile(sample = K3,
                                            profile_name = 'K3 initial R || c',
                                            direction = 'c', raypath = 'a',
                                            base_low_wn = 3500,
                                            fname_list = ['K3_sp11_Adc', 
                                                          'K3_sp12_Bdc'],
                                            positions_microns = [leng/2., 
                                                                 leng/2.])

K3wb_trueInit = nams.WholeBlock(name = ''.join(('Initial measurements for',
                                                ' Kunlun diopside K3 along',
                                                ' 3 ray paths')),
                                profiles = [profile_K3_trueInit_raypathA, 
                                            profile_K3_trueInit_raypathB, 
                                            profile_K3_trueInit_raypathC],
                                time_seconds = 10., 
                                worksheetname = 'K3 initial')

#%% K3 heated at 696 C for 2 hours "pre-anneal", "initial" profiles
leng = np.mean(K3.twoA_list)
profile_K3_init_A = nams.Profile(sample = K3, 
                                 profile_name = 'K3 heated at 696 C for 2hr || a*',
                                 initial_profile = profile_K3_trueInit_raypathB,
                                 direction = 'a', raypath = 'b',
                                 fname_list = ['K3_cdb06'],
                                 positions_microns = [leng/2.])

leng = np.mean(K3.twoB_list)
profile_K3_init_B = nams.Profile(sample = K3,
                                 profile_name = ''.join(('K3 heated at 696 C',
                                                 'for 2hr || b, OFF CENTER')),
                                 initial_profile = profile_K3_trueInit_raypathC,
                                 direction = 'b', raypath = 'c', 
                                 base_low_wn=3500.,
                                 fname_list = ['K3_bdc01', 'K3_bdc02', 'K3_bdc03', 
                                               'K3_bdc04', 'K3_bdc05', 'K3_bdc06', 
                                               'K3_bdc07', 'K3_bdc08', 'K3_bdc09', 
                                               'K3_bdc10'],
                                 positions_microns = [])
for k in range(10):
    profile_K3_init_B.positions_microns.append((k+1)*leng/11)

leng = np.mean(K3.twoC_list)
profile_K3_init_C = nams.Profile(sample = K3,
                                 profile_name = 'K3 heated at 696 C for 2hr || c',
                                 initial_profile = profile_K3_trueInit_raypathB,
                                 direction = 'c', raypath = 'b',
                                 fname_list = ['K3_cdb01', 'K3_cdb02', 'K3_cdb03', 
                                               'K3_cdb04', 'K3_cdb05', 'K3_cdb06', 
                                               'K3_cdb07', 'K3_cdb08', 'K3_cdb09', 
                                               'K3_cdb10', 'K3_cdb11', 'K3_cdb12', 
                                               'K3_cdb13'],
                                 positions_microns = np.linspace(50, leng-50, 13))

K3wb_init = nams.WholeBlock(name = ''.join(('"Initial" (pre-anneal at 696 C',
                                            ' for 2 hrs) profiles for Kunlun',
                                            ' diopside K3')),
                            profiles = [profile_K3_init_A, 
                                        profile_K3_init_B, 
                                        profile_K3_init_C],
                            time_seconds = 2.*3600., 
                            worksheetname = 'K3 696C 2hr')

#%% K3 heated at 696 C for an additional 17hr 15m = 19hr 15m total
leng = np.mean(K3.twoA_list)
profile_K3_696C_19hr_A = nams.Profile(sample = K3,
                                      profile_name = 'K3 696 C for 19hr 15m || a',
                                      initial_profile = profile_K3_init_A,
                                      direction = 'a', raypath = 'b',
                                      fname_list = ['K3h_cdb02'],
                                      positions_microns = [leng/2.])

leng = np.mean(K3.twoB_list)
profile_K3_696C_19hr_B = nams.Profile(sample = K3,
                                      profile_name = 'K3 696 C for 19hr 15m || b OFF-CENTER',
                                      initial_profile = profile_K3_init_B,
                                      direction = 'b', raypath = 'c',
                                      base_low_wn=3500.,
                                      fname_list = ['K3h_bdc01', 'K3h_bdc02', 
                                                    'K3h_bdc03',
                                                    'K3h_bdc04', 'K3h_bdc05'],
                                      positions_microns = [leng/5., 2.*leng/5., 
                                                           leng/2., 3.*leng/5., 
                                                           4.*leng/5.])

leng = np.mean(K3.twoC_list)
profile_K3_696C_19hr_C = nams.Profile(sample = K3, 
                                      profile_name = 'K3 696C for 19hr 15m || c',
                                      initial_profile = profile_K3_init_C,
                                      direction = 'c', raypath = 'b',
                                      fname_list = ['K3h_cdb01', 'K3h_cdb02', 
                                                    'K3h_cdb03'],
                                      positions_microns = [leng/4., leng/2., 
                                                           3.*leng/4.])

K3wb_700C_19hr = nams.WholeBlock(name = 'Kunlun diopside K3 at 696 C for 19hr 15 min',
                                 profiles = [profile_K3_696C_19hr_A, profile_K3_696C_19hr_B,
                                             profile_K3_696C_19hr_C],
                                 time_seconds = (19.*3600.) + (15.*60.),
                                 worksheetname = 'K3 696C 19hr')

#%% K3 heated at 696 C for an additional 16hr = 35hr 15m total
leng = np.mean(K3.twoA_list)
profile_K3_696C_35hr_A = nams.Profile(sample = K3,
                                      profile_name = 'K3 696C for 35hr 15m || a',
                                      initial_profile = profile_K3_init_A,
                                      direction = 'a', raypath = 'b',
                                      fname_list = ['K3h700_cdb01'],
                                      positions_microns = [leng/2.])

leng = np.mean(K3.twoB_list)
profile_K3_696C_35hr_B = nams.Profile(sample = K3,
                                      profile_name = 'K3 696C for 35hr 15m || b, OFF CENTER',
                                      initial_profile = profile_K3_init_B,
                                      direction = 'b', raypath = 'c',
                                      base_low_wn = 3500,
                                      fname_list = ['K3h700_bdc01', 
                                                    'K3h700_bdc02', 
                                                    'K3h700_bdc03'], 
                                     positions_microns = [leng/4., leng/2., 
                                                          3.*leng/4.])

leng = np.mean(K3.twoC_list)
profile_K3_696C_35hr_C = nams.Profile(sample = K3,
                                      profile_name = 'K3 696C for 35hr 15m || c',
                                      initial_profile = profile_K3_init_C,
                                      direction = 'c', raypath = 'b',
                                      fname_list = ['K3h700_cdb01', 
                                                    'K3h700_cdb02',
                                                    'K3h700_cdb03'],
                                      positions_microns = [leng/4., leng/2., 
                                                           3.*leng/4.])

K3wb_700C_35hr = nams.WholeBlock(name = 'Kunlun diopside K3 at 696 C for 35hr 15 min',
                                 profiles = [profile_K3_696C_35hr_A, profile_K3_696C_35hr_B,
                                             profile_K3_696C_35hr_C],
                                 time_seconds = (19.*3600.) + (15.*60.) + (16.*3600.), 
                                 worksheetname = 'K3 696C 35hr')

#%% K3 heated at 796 C for 15 hr 40 m
leng = np.mean(K3.twoA_list)
profile_K3_800C_15hr40m_A = nams.Profile(initial_profile = profile_K3_init_A,
                                         sample = K3,
                                         profile_name=''.join(('Kunlun diopside',
                                             ' K3 at 796C for 15hr 40m || a')),
                                         direction = 'a', raypath = 'b',
                                         fname_list = ['K3h800_cdb02'],
                                         positions_microns = [leng/2.])

leng = np.mean(K3.twoB_list)
profile_K3_800C_15hr40m_B = nams.Profile(initial_profile = profile_K3_init_B,
                                         sample = K3, 
                                         profile_name = ('Kunlun diopside K3 at 796C for 15hr 40m || b'),
                                         direction = 'b', raypath = 'c',
                                         base_low_wn = 3500.,
                                         fname_list = ['K3h800_bdc01', 'K3h800_bdc02', 
                                                       'K3h800_bdcMID'],
                                         positions_microns = [50., 150., 
                                                              leng/2.])

leng = np.mean(K3.twoC_list)
profile_K3_800C_15hr40m_C = nams.Profile(initial_profile = profile_K3_init_C,
                                         sample = K3,
                                         profile_name = 'Kunlun diopside K3 at 796C for 15hr 40m || c',
                                         direction = 'c', raypath = 'b',
                                         fname_list = ['K3h800_cdb01', 'K3h800_cdb02', 
                                                       'K3h800_cdbMID'],
                                         positions_microns = [50., 150., 
                                                              leng/2.])

K3wb_800C_15hr = nams.WholeBlock(name = 'Kunlun diopside K3 at 796 C for 15hr 40min',
                                 profiles = [profile_K3_800C_15hr40m_A, 
                                             profile_K3_800C_15hr40m_B,
                                             profile_K3_800C_15hr40m_C],
                                 time_seconds = (15.*3600) + 40.,
                                 worksheetname = 'K3 800C 15hr')

#%% K3 heated at 817 C for 6 days
leng = np.mean(K3.twoA_list)
profile_K3_817C_6days_a = nams.Profile(sample = K3,
                                       profile_name = 'K3 initial || a*',
                                       short_name = 'K3g_adb',
                                       initial_profile = profile_K3_init_A,
                                       direction = 'a', raypath = 'b',
                                       fname_list = ['K3g_adb01', 'K3g_adb02', 
                                                     'K3g_adb05', 'K3g_adb07', 
                                                     'K3g_adb10', 'K3g_adb15', 
                                                     'K3g_adb19', 'K3g_adb20'],
                                       positions_microns = [50., 150., 450., 
                                                            650., 950., 1450., 
                                                            leng-150., leng-50.])

leng = np.mean(K3.twoB_list)
profile_K3_817C_6days_b = nams.Profile(sample = K3,
                                       profile_name = 'K3 initial || b',
                                       short_name = 'K3g_bdc',
                                       initial_profile = profile_K3_init_B,
                                       direction = 'b' , raypath = 'c',
                                       base_low_wn=3500., 
                                       fname_list = ['K3g_bdc01', 'K3g_bdc02', 
                                                     'K3g_bdc06', 
                                                     'K3g_bdc10', 'K3g_bdc12'],
                                       positions_microns = [1150., 950., 550.,
                                                            150., 50.])

leng = np.mean(K3.twoC_list)
profile_K3_817C_6days_c = nams.Profile(sample = K3,
                                       profile_name = 'K3 initial || c',
                                       short_name = 'K3g_cdb',
                                       initial_profile = profile_K3_init_C,
                                       direction = 'c' , raypath = 'b',
                                       fname_list = ['K3g_cdb01', 'K3g_cdb02', 
                                                     'K3g_cdb03', 
                                                     'K3g_cdb05', 'K3g_cdb06', 
                                                     'K3g_cdb08', 
                                                     'K3g_adb10', 'K3g_cdb12', 
                                                     'K3g_cdb14', 
                                                     'K3g_cdb17', 'K3g_cdb18'],
                                         # CHECK POSITION ORIENTATION HERE
                                         positions_microns = [50., 150., 250., 
                                                              450., 550., 750.,
                                                              950., 1150., 1350., 
                                                              leng-150., leng-50.])

K3wb_6days = nams.WholeBlock(name = 'Kunlun diopside K3 heated at 817 C for 6 days',
                             profiles = [profile_K3_817C_6days_a, 
                                         profile_K3_817C_6days_b, 
                                         profile_K3_817C_6days_c],
                             time_seconds = 6.*24.*3600.,
                             temperature_celsius = 817.,
                             worksheetname = 'K3 817C 6days')

#%% K4 Kunlun diopside true initial profiles K4*I
leng = np.mean(K4.twoA_list)
profile_K4_init_A = nams.Profile(sample = K4,
                                 fname_list = ["K4_adcIL", "K4_dcImid", 
                                               "K4_adcIr4"],
                                 base_low_wn = 3500.,
                                 direction = 'a', raypath = 'c',
                                 positions_microns = np.array([1000, leng/2.0, 
                                                               leng-1000]),
                                 profile_name = 'K4 initial || a')

leng = np.mean(K4.twoB_list)
profile_K4_init_B = nams.Profile(sample = K4, direction = 'b', raypath = 'c',
                                 base_low_wn = 3500.,
                                 profile_name = 'K4 initial || b',
                                 positions_microns = np.array([leng/2.0, 
                                                               leng/6.0, 
                                                               leng*(2./3.), 
                                                               leng*(5./6.)]), 
                                fname_list = ['K4_dcImid', 'K4_dcIleft', 'K4_dcIr2', 
                                              'K4_dcIright'])

leng = np.mean(K4.twoC_list)
profile_K4_init_C = nams.Profile(sample = K4, profile_name = 'K4 initial || c',
                                 direction = 'c', raypath = 'b',
                                 fname_list = ['K4_cdbIleft', 'K4_cdbImid', 
                                               'K4_cdbIright'],
                                 positions_microns = [leng/6.0, leng/2.0, 
                                                      leng*(5.0/6.0)])

K4wb_init = nams.WholeBlock(name = 'Initial profiles for Kunlun diopside K4',
                            profiles = [profile_K4_init_A, 
                                        profile_K4_init_B, 
                                        profile_K4_init_C],
                            time_seconds = 10., 
                            worksheetname = 'Kunlun 904C initial')

#%% K4 0-time experiement (quenched immediately by falling to bottom) ###
### Used as initial for other experiments ### K4*Q
leng = np.mean(K4.twoA_list)
profile_K4_quench_A = nams.Profile(sample = K4, 
                                   profile_name = 'Kunlun 0-time experiment || a',
                                   direction = 'a', raypath = 'c',
                                   initial_profile = profile_K4_init_A,
                                   fname_list = ['K4q_adc05', 'K4q_bdcMID', 
                                                 'K4q_adc65'],
                                   positions_microns = [525., leng/2., 6525.],
                                   base_low_wn=3500.)

leng = np.mean(K4.twoB_list)
profile_K4_quench_B = nams.Profile(sample = K4, 
                                   profile_name = 'Kunlun 0-time experiment || b',
                                   direction = 'b', raypath = 'c', 
                                   initial_profile = profile_K4_init_B,
                                   fname_list = ['K4q_bdc01', 'K4q_bdcMID', 
                                                 'K4q_bdc02'],
                                   positions_microns = [120., leng/2., 
                                                        leng-120.],
                                   base_low_wn=3500.)

leng = np.mean(K4.twoC_list)
profile_K4_quench_C = nams.Profile(sample = K4, 
                                   profile_name = 'Kunlun 0-time experiment || c',
                                   direction = 'c', raypath = 'b',
                                   initial_profile = profile_K4_init_C,
                                   fname_list = ['K4q_cdb01', 'K4q_cdbMID', 
                                                 'K4q_cdb02', 'K4q_cdb03'],
                                   positions_microns = [100., leng/2., 
                                                        leng-110., leng-220.])

K4wb_quench = nams.WholeBlock(name = '0-time (37m at 480C) profiles at 904 C for Kunlun diopside K4',
                              profiles = [profile_K4_quench_A, 
                                          profile_K4_quench_B, 
                                          profile_K4_quench_C],
                              time_seconds = 30.,
                              worksheetname = 'Kunlun 904C 0-time')

%# K4 heated 1 hour at 904C, K4*H
leng = np.mean(K4.twoA_list)
profile_K4_904C_1hr_A = nams.Profile(sample = K4, direction = 'a', raypath = 'c',
                                     profile_name = 'Kunlun heated 1 hr at 904 C || a',
                                     fname_list = ['K4h_adc01i', 'K4h_adc02', 
                                                   'K4h_adc15', 'K4h_adc22', 
                                                   'K4h_adc29', 'K4h_dcMID',
                                                   'K4h_adc42', 'K4h_adc56', 
                                                   'K4h_adc66',],
                                                   #'K4h_adc68e'
                                    base_low_wn=3500,
                                    positions_microns = [100, 200, 1500, 2200, 2900, 
                                                         leng/2.0, 4200, 5600, 6600, 
                                                         #leng-100
                                                         ],
                                    initial_profile = profile_K4_quench_A)

leng = np.mean(K4.twoB_list)
profile_K4_904C_1hr_B = nams.Profile(sample = K4, direction = 'b', raypath = 'c',
                                     profile_name = 'Kunlun heated 1 hr at 904 C || b',
                                     base_low_wn=3500.,
                                     initial_profile = profile_K4_quench_B,
                                     fname_list = ['K4h_bdc01e', 'K4h_bdc01i', 
                                                   'K4h_bdc02', 
                                                   'K4h_bdc04', 'K4h_bdc06', 'K4h_dcMID', 
                                                   'K4h_bdc13', 'K4h_bdc15', 'K4h_bdc17', 
                                                   'K4h_bdc20', 'K4h_bdc21i', 'K4h_bdc21e'],
                                    positions_microns = [50., 100., 200., 
                                                         400., 600., leng/2., 
                                                         1300., 1500., 1700., 
                                                         leng-200, leng-100., leng-50.])

leng = np.mean(K4.twoC_list)
profile_K4_904C_1hr_C = nams.Profile(sample = K4, direction = 'c', raypath = 'b',
                                     profile_name = 'Kunlun heated 1 hr at 904 C || c',
                                     initial_profile = profile_K4_quench_C,
                                     fname_list = ['K4h_cdb01', 'K4h_cdb02', 
                                                   'K4h_cdb03', 
                                                   'K4h_cdb04', 'K4h_cdb08', 'K4h_cdb10', 
                                                   'K4h_cdb12', 'K4h_cdb14', 'K4h_cdb15', 
                                                   'K4h_cdb16'],
                                    positions_microns = [115., 215., 315., 
                                                         415., leng/2., 1015., 
                                                         1215., 1415., 1515., 
                                                         leng-115])

K4wb_1hr = nams.WholeBlock(name = 'Kunlun diopside K4 heated at 904 C for 1 hour',
                           profiles = [profile_K4_904C_1hr_A, 
                                       profile_K4_904C_1hr_B, 
                                       profile_K4_904C_1hr_C], 
                           time_seconds = 3600.,
                           worksheetname = 'Kunlun 904C 1 hour')

#%% K4 heated 91 hours (K4P) ###
profile_K4_904C_91hr_A = nams.Profile(sample = K4, 
                                      profile_name = 'Kunlun K4 heated 904 C for 91 hr || a',
                                      short_name = 'K4adc_P',
                                      direction = 'a', raypath = 'c',
                                      base_low_wn=3500.,
                                      initial_profile = profile_K4_quench_A,
                                      fname_list = ['K4p_adc01', 'K4p_adc02', 'K4p_adc08', 
                                                    'K4p_adc14', 'K4p_adc22', 'K4p_adc29', 
                                                    'K4p_adc35', 'K4p_adc42', 'K4p_adc49', 
                                                    'K4p_adc56', 'K4p_adc61', 'K4p_adc67'],
                                      positions_microns = [100., 200., 800., 
                                                           1400., 2200., 2900., 
                                                           3500., 4200., 4900., 
                                                           5600., 6100., 6700.])

profile_K4_904C_91hr_B = nams.Profile(sample = K4,
                                      profile_name = 'Kunlun K4 heated 904 C for 91 hr || b',
                                      short_name = 'K4bdc_P', direction = 'b',
                                      raypath = 'c',
                                      base_low_wn=3500,
                                      initial_profile = profile_K4_quench_B,
                                      fname_list = ['K4p_bdc01', 'K4p_bdc02', 'K4p_bdc04', 
                                                    'K4p_bdc06', 'K4p_bdc09', 'K4p_adc35', 
                                                    'K4p_bdc13', 'K4p_bdc15', 'K4p_bdc17', 
                                                    'K4p_bdc19', 'K4p_bdc21'],
                                      positions_microns = [100., 200., 400., 
                                                           600., 900., 1200., 
                                                           1400., 1500., 1700., 
                                                           1900., 2100.])

leng = np.mean(K4.twoC_list)
profile_K4_904C_91hr_C = nams.Profile(sample = K4,
                                      profile_name = 'Kunlun K4 heated 904 C for 91 hr || c',
                                      short_name = 'K4cdb_P', direction = 'c',
                                      raypath = 'b',
                                      initial_profile = profile_K4_quench_C,
                                      fname_list = ['K4p_cdb01', 'K4p_cdb02', 'K4p_cdb04', 
                                                    'K4p_cdb06', 'K4p_cdb09', 'K4p_cdb11', 
                                                    'K4p_cdb13', 'K4p_cdb15', 'K4p_cdb16'],
                                      positions_microns = [50., 150., 400., 
                                                           600., 900., 1100., 
                                                           1300., leng-150., leng-50.])

K4wb_91hr = nams.WholeBlock(name = 'Kunlun diopside K4 heated at 904 C for 91 hours',
                            worksheetname = 'Kunlun 904C 91 hours',
                            temperature_celsius = 904.,
                            profiles = [profile_K4_904C_91hr_A, 
                                        profile_K4_904C_91hr_B, 
                                        profile_K4_904C_91hr_C],
                            time_seconds = 91.*3600.)

#%% K4 heated 154 hours ###
profile_K4_904C_154hr_A = nams.Profile(sample = K4, direction = 'a',
                                       raypath = 'c', base_low_wn=3500.,
                                       initial_profile = profile_K4_quench_A,
                                       profile_name = 'Kunlun heated 154 hours at 904C || a*',
                                       short_name = 'K4adcF',
                                       fname_list = ['K4f_adc01', 'K4f_adc02', 'K4f_adc04',
                                                     'K4f_adc06', 'K4f_adc08', 'K4f_adc11',
                                                     'K4f_adc14', 'K4f_adc18', 'K4f_adc22',
                                                     'K4f_adc25', 'K4f_adc29', 'K4f_adc35',
                                                     'K4f_adc42', 'K4f_adc49', 'K4f_adc56',
                                                     'K4f_adc61', 'K4f_adc67', 'K4f_adc68'],
                                       positions_microns = [100, 200, 400, 
                                                            600, 800, 1100,
                                                            1400, 1800, 2200, 
                                                            2500, 2900, 3500, 
                                                            4200, 4900, 5600, 
                                                            6100, 6700, 6800])

leng = np.mean(K4.twoB_list)
profile_K4_904C_154hr_B = nams.Profile(sample = K4, direction = 'b',
                                       raypath = 'c', base_low_wn=3500.,
                                       initial_profile = profile_K4_quench_B,
                                       profile_name = 'Kunlun heated 154 hours at 904C || b',
                                       short_name = 'K4bdcF',
                                       fname_list = ['K4f_bdc01', 'K4f_bdc02', 'K4f_bdc03', 
                                                     'K4f_bdc04', 'K4f_bdc05', 'K4f_bdc06', 
                                                     'K4f_bdc07', 'K4f_bdc08', 'K4f_bdc09', 
                                                     'K4f_bdc10', 'K4f_adc35', 'K4f_bdc15', 
                                                     'K4f_bdc17', 'K4f_bdc19', 'K4f_bdc21', 
                                                     'K4f_bdc22'],
                                        positions_microns = [100., 200., 300., 
                                                             400., 300., 600., 
                                                             700., 800., 900., 
                                                             1000., 1200., 1500., 
                                                             1700., 1900., leng-150., 
                                                             leng-50.])

leng = np.mean(K4.twoC_list)                                             
profile_K4_904C_154hr_C = nams.Profile(sample = K4, direction = 'c',
                                       raypath = 'b',
                                       initial_profile = profile_K4_quench_C,
                                       profile_name = 'Kunlun heated 154 hours at 904C || c',
                                       short_name = 'K4cdbF',
                                       fname_list = ['K4f_cdb01', 'K4f_cdb02', 'K4f_cdb03',
                                                     'K4f_cdb04', 'K4f_cdb05', 'K4f_cdb06',
                                                     'K4f_cdb07', 'K4f_cdb09', 'K4f_cdb11',
                                                     'K4f_cdb12', 'K4f_cdb13', 'K4f_cdb14',
                                                     'K4f_cdb16'],
                                       positions_microns = [50, 120, 300, 400, 500, 600, 700,
                                                            900, 1100, 1200, 1300, 1400, 
                                                            leng-50])

K4wb_154hr = nams.WholeBlock(name = 'Kunlun diopside K4 heated at 904 C for 154 hours',
                             worksheetname = 'Kunlun 904C 154 hours',
                             temperature_celsius = 904.,
                             profiles = [profile_K4_904C_154hr_A, 
                                         profile_K4_904C_154hr_B, 
                                         profile_K4_904C_154hr_C],
                             time_seconds = 154.*3600.)

#%% K5 Kunlun diopside initial profiles
leng = np.mean(K5.twoA_list)
profile_K5_initial_A = nams.Profile(sample = K5, direction = 'a', raypath = 'b',
                                    profile_name = 'Kunlun K5 initial || a',
                                    fname_list = ['K5i_adb01', 'K5i_dbMID', 'K5i_adb32'],
                                    positions_microns = [100., leng/2., leng-200])

leng = np.mean(K5.twoB_list)
profile_K5_initial_B = nams.Profile(sample = K5, direction = 'b', raypath = 'c',
                                    fname_list = ['K5i_bdc01', 'K5i_bdcMID', 
                                                  'K5i_bdc16'],
                                    profile_name = 'Kunlun K5 initial || b',
                                    positions_microns = [50., leng/2., leng-50],
                                    base_low_wn=3500.)

leng = np.mean(K5.twoC_list)
profile_K5_initial_C = nams.Profile(sample = K5, direction = 'c', raypath = 'b',
                                    profile_name = 'Kunlun K5 initial || c',
                                    fname_list = ['K5i_cdb01', 'K5i_dbMID', 'K5i_cdb16'],
                                    positions_microns = [100., leng/2., 1600.])

K5wb_init = nams.WholeBlock(name = 'Kunlun diopside K5 initial',
                            worksheetname = 'Kunlun 1000C initial',
                            profiles = [profile_K5_initial_A,
                                        profile_K5_initial_B, 
                                        profile_K5_initial_C],
                            time_seconds = 10.)

#%% K5 heated at 1000C for 75 hours;  whole-block
## Diffusivities from Ferriss et al. 2015 American Mineralogist
profile_K5_1000C_75hr_A = nams.Profile(sample = K5, direction = 'a',
                                       raypath = 'b', time_seconds = 75.*3600.,
                                       profile_name = 'Kunlun heated at 1000C for 75hr || a*',
                                       short_name = 'K5adb',
                                       initial_profile = profile_K5_initial_A,
                                       diffusivity_log10m2s = -13.04,
                                       diff_error = 0.12,
                                       peak_diffusivities = [0., -13.27, -13.08, 0., 0., 0.],
                                       peak_diff_error = [0., 0.29, 0.25, 0., 0., 0.],
                                       fname_list = ['K5_adb01', 'K5_adb02', 'K5_adb03', 
                                                     'K5_adb04', 'K5_adb05', 'K5_adb06',
                                                     'K5_adb07', 'K5_adb09', 'K5_adb11',
                                                     'K5_adb14', 'K5_cdb09', 'K5_adb21',
                                                     'K5_adb24', 'K5_adb25', 'K5_adb26',
                                                     'K5_adb27', 'K5_adb28', 'K5_adb29',
                                                     'K5_adb30', 'K5_adb31'],
                                       positions_microns = [50., 150., 300., 
                                                            400., 500., 600.,
                                                            700., 900., 1100., 
                                                            1400., 1800., 2100., 
                                                            2400., 2500., 2600., 
                                                            2700., 2800., 2900.,
                                                            3000., 3100.])

leng = np.mean(K5.twoB_list)
profile_K5_1000C_75hr_B = nams.Profile(sample = K5, direction = 'b', 
                                       raypath = 'c', time_seconds = 75.*3600.,
                                       profile_name = 'Kunlun heated at 1000C for 75hr || b',
                                       short_name = 'K5bdc', base_low_wn=3500.,
                                       initial_profile = profile_K5_initial_B,
                                       diffusivity_log10m2s = -13.41,
                                       diff_error = 0.08,
                                       peak_diffusivities = [0., -12.38, -12.07, 0., 0., 0.],
                                       peak_diff_error = [0., 0.03, 0.02, 0., 0., 0.],
                                       fname_list = ['K5_bdc01', 'K5_bdc02', 'K5_bdc03',
                                                     'K5_bdc04', 'K5_bdc05', 'K5_bdc06',
                                                     'K5_bdc07', 'K5_bdc08', 'K5_bdc09',
                                                     'K5_bdc10', 'K5_bdc11', 'K5_bdc12',
                                                     'K5_bdc13', 'K5_bdc14', 'K5_bdc15',
                                                     'K5_bdc16'],
                                       positions_microns = [50., 150., 300., 
                                                            400., 500., 600.,
                                                            700., 800., 900., 
                                                            1000., 1100., 1200., 
                                                            1300., 1400., leng-150., 
                                                            leng-50.])

profile_K5_1000C_75hr_C = nams.Profile(sample = K5, direction = 'c', 
                                       raypath = 'b', time_seconds = 75.*3600.,
                                       profile_name = 'Kunlun heated at 1000C for 75hr || c',
                                       short_name = 'K5cdb',
                                       diffusivity_log10m2s = -13.58,
                                       diff_error = 0.12,
                                       peak_diffusivities = [0., -13.36, -13.08, 0., 0., 0.],
                                       peak_diff_error = [0., 0.18, 0.17, 0., 0., 0.],
                                       initial_profile = profile_K5_initial_C,
                                       fname_list = ['K5_cdb01', 'K5_cdb02', 'K5_cdb03', 
                                                     'K5_cdb04', 'K5_cdb05', 'K5_cdb06', 
                                                     'K5_cdb07', 'K5_cdb08', 'K5_cdb09', 
                                                     'K5_cdb10', 'K5_cdb11', 'K5_cdb12', 
                                                     'K5_cdb13', 'K5_cdb14', 'K5_cdb15', 
                                                     'K5_cdb16'],
                                       positions_microns = [50., 150., 300., 
                                                            400., 500., 600., 
                                                            700., 800., 900., 
                                                            1000., 1100., 1200., 
                                                            1300., 1400., 1500., 
                                                            1600])

K5wb = nams.WholeBlock(name = 'Kunlun diopside K5 heated at 1000 C for 75 hours',
                       worksheetname = 'Kunlun 1000C 75hr',
                       temperature_celsius = 1000.,
                       profiles = [profile_K5_1000C_75hr_A,
                                   profile_K5_1000C_75hr_B, 
                                   profile_K5_1000C_75hr_C],
                    time_seconds = 75.*3600.)

#%% profile K5 SLICE K5sbda
#initial state used K3_sp13_Bda.CSV with thickness transforms
profile_K5_slice_1000C_75hr_B = nams.Profile(sample = K5slice, 
                                             direction = 'b',
                                             raypath = 'a',
                                             profile_name = 'Kunlun slice data || b',
                                             fname_list = ['K5c_bda01', 'K5c_bda02', 'K5c_bda03', 'K5c_bda04', 
                                                           'K5c_bda06', 'K5c_cda09', 'K5c_bda10',              
                                                           'K5c_bda13', 'K5c_bda14', 'K5c_bda15', 'K5c_bda16'],
                                             positions_microns = [ 50., 140., 290., 390., 
                                                                  590.,       790., 
                                                                  990.,
                                                                  1290., 1390., 1490., 1538.])

leng = np.mean(K5slice.twoC_list)
profile_K5_slice_1000C_75hr_C = nams.Profile(sample = K5slice, direction = 'c',
                                             raypath = 'a',
                                             profile_name = 'Kunlun slice data || c',
                                             fname_list = ['K5c_cda01', 'K5c_cda02', 'K5c_cda03', 'K5c_cda04', 
                                                           'K5c_cda05', 'K5c_cda06', 'K5c_cda07', 'K5c_cda08', 
                                                           'K5c_cda09', 'K5c_cda10', 'K5c_cda11', 'K5c_cda12', 
                                                           'K5c_cda13', 'K5c_cda14', 'K5c_cda15', 'K5c_cda16', 
                                                           'K5c_cda17'],
                                             positions_microns = [49.1, 149.1, 249.1, 349.1, 
                                                                  449.1, 549.1, 649.1, 749.1, 
                                                                  849.1, 949.1, 1049.1, 1149.1, 
                                                                  1249.1, 1349.1, 1446.7, 1546.7, 
                                                                  leng-50])


#%% K6: polarized IR spectra on untreated Kunlun diopside

# E || a mean absorbance of 4 spectra
K6_db_Ea = nams.Spectrum(sample=K6, thick_microns=K6.sample_thick_microns[1],
                         fname = 'K6_db_Ea', polar = 'E || a')
K6_db_Ea_2 = nams.Spectrum(sample=K6, thick_microns=K6.sample_thick_microns[1],
                           fname = 'K6_db_Ea_2', polar = 'E || a')
K6_dc_Ea = nams.Spectrum(sample=K6, thick_microns=K6.sample_thick_microns[2],
                         fname = 'K6_dc_Ea', polar = 'E || a')
K6_dc_Ea_2 = nams.Spectrum(sample=K6, thick_microns=K6.sample_thick_microns[2],
                           fname = 'K6_dc_Ea_2', polar = 'E || a')
ave_K6_Ea = nams.Spectrum(sample=K6, fname = 'ave_K6_Ea',
                          other_name = 'Average of 4 initial K6 spectra, E || a',
                          polar = 'E || a', folder=folder)
ave_K6_Ea.make_average_spectra([K6_db_Ea, K6_db_Ea_2, K6_dc_Ea, K6_dc_Ea_2],)

# E || b = mean absorbance of 2 spectra
K6_dc_Eb = nams.Spectrum(sample=K6, thick_microns=K6.sample_thick_microns[2],
                         fname = 'K6_dc_Eb', polar = 'E || b')
K6_dc_Eb_2 = nams.Spectrum(sample=K6, thick_microns=K6.sample_thick_microns[2],
                           fname = 'K6_dc_Eb_2', polar = 'E || b')
ave_K6_Eb = nams.Spectrum(sample=K6, fname = 'ave_K6_Eb',
                          other_name = 'Average of 2 initial K6 spectra, E || b',
                          polar = 'E || b', folder=folder)
ave_K6_Eb.make_average_spectra([K6_dc_Eb, K6_dc_Eb_2])

# E || c only one spectrum
K6_db_Ec = nams.Spectrum(sample=K6, thick_microns=K6.sample_thick_microns[1],
                         fname = 'K6_db_Ec_2', polar = 'E || c')
ave_K6_Ec = nams.Spectrum(sample=K6, fname = 'ave_K6_Ec',
                          other_name = 'Average of initial K6 spectra, E || c',
                          polar = 'E || c', folder=folder)
ave_K6_Ec.make_average_spectra([K6_db_Ec])

#%% diopside from north of Jaipur, India; provided by Stephen Mackwell

# Polarized IR spectra digitized from Woods et al. 2000
J_Ea = nams.Spectrum(thick_microns=0.8e6, fname = 'J_Ea', polar = 'a') 
J_Eb = nams.Spectrum(thick_microns=0.8e6, fname = 'J_Eb', polar = 'b')
J_Ec = nams.Spectrum(thick_microns=0.8e6, fname = 'J_Ec', polar = 'c') # c*

# Polarized IR spectra measured on untreated, unoriented sample J2
J2_Ea = nams.Spectrum(fname='J2_dshort_Emid2', thick_microns=np.mean(J2.twoB_list))
J2_Eb = nams.Spectrum(fname='J2_dmid_Eshort', thick_microns=np.mean(J2.twoA_list))
J2_Ec = nams.Spectrum(fname='J2_dshort_Elong2', thick_microns=np.mean(J2.twoB_list))

#%% J1 initial profiles
leng = np.mean(J1.twoA_list)
profile_J1_init_A = nams.Profile(fname_list = ['J1_dcI_t', 'J1_sp2', 'J1_dcI_bot'],
                                 direction = 'a', raypath = 'b', sample=J1,
                                 positions_microns = [leng*0.25, leng/2.0, leng*0.75],
                                 profile_name = 'Jaipur diopside J1 initial || a*')

leng = np.mean(J1.twoB_list)
profile_J1_init_B = nams.Profile(fname_list = ['J1_sp3'], direction = 'b',
                                 raypath = 'a', positions_microns = [leng/2.0],
                                 base_high_wn = 3600., sample=J1,
                                 profile_name = 'Jaipur diopside J1 initial || b')

leng = np.mean(J1.twoC_list)                                 
profile_J1_init_C = nams.Profile(sample=J1, direction = 'c', raypath = 'b',
                                 fname_list = ['J1_dcI_L', 'J1_sp2', 'J1_dcI_R'],
                                 positions_microns = [leng*0.25, leng/2.0, leng*0.75],
                                 profile_name = 'Jaipur diopside J1 initial || c')

J1wb_initial = nams.WholeBlock(name = 'Jaipur diopside initial',
                               time_seconds = 10., 
                               worksheetname = 'Jaipur 904C initial',
                               profiles = [profile_J1_init_A, 
                                           profile_J1_init_B, 
                                           profile_J1_init_C])

#%% J1 after heating 34 minutes at 904 C       
leng = np.mean(J1.twoA_list)
profile_J1_904C_30m_A = nams.Profile(sample=J1,
                                     profile_name = 'Jaipur diopside at 904 C for 30 m || a*',
                                     short_name = 'J1a', direction = 'a',
                                     raypath = 'b', time_seconds = 30.*60.,
                                     initial_profile = profile_J1_init_A,
                                     fname_list = ['J1_bdc01', 'J1_bdc02', 'J1_bdc06', 
                                                   'J1_bdc10', 'J1_bdc14', 'J1_bdc18', 
                                                   'J1_dcMID', 'J1_bdc25', 'J1_bdc29', 
                                                   'J1_bdc34', 'J1_bdc39', 'J1_bdc40', 
                                                   'J1_bdc41', 'J1h_bdc05', 'J1h_bdc06', 
                                                   'J1h_bdc07', 'J1h_bdc16', 'J1h_bdc17', 
                                                   'J1h_bdc18', 'J1h_bdc19', 'J1h_bdc21', 
                                                   'J1h_bdc22', 'J1h_bdc24', 'J1h_bdc26', 
                                                   'J1h_bdc28',               'J1h_bdc33', 
                                                   'J1h_bdc35', 'J1h_bdc37', 'J1h_bdc40', 
                                                   'J1h_bdc41', 'J1h_bdc42'],
                                     positions_microns = [100., 200., 600., 
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
                                                          4100., 4200.])

leng = np.mean(J1.twoB_list)
profile_J1_904C_30m_B = nams.Profile(sample=J1,
                                     profile_name = 'Jaipur diopside at 904 C for 30 m || b',
                                     short_name = 'J1b', direction = 'b',
                                     raypath = 'a', base_high_wn = 3600.,
                                     fname_list = ['J1_cdb01', 'J1_cdb02', 'J1_cdb05', 
                                                   'J1_cdb08', 'J1_cdb12', 'J1_cdb15', 
                                                   'J1_cdb18', 'J1_cdb22', 'J1_cdb25', 
                                                   'J1_cdb28'],
                                     positions_microns = [100., 200., 500., 
                                                          800., 1200., leng/2., 
                                                          1800., 2200., leng-600., 
                                                          2800.],
                                     initial_profile = profile_J1_init_B)

leng = np.mean(J1.twoC_list)
profile_J1_904C_30m_C = nams.Profile(profile_name = 'Jaipur diopside at 904 C for 30 m || c',
                                     short_name = 'J1c', direction = 'c',
                                     raypath = 'b', sample=J1,
                                     fname_list = ['J1_adc01', 'J1_adc02', 'J1_adc05', 
                                                   'J1_adc10', 'J1_dcMID', 'J1_adc15', 
                                                   'J1_adc20', 'J1_adc24', 'J1_adc25', 
                                                   'J1h_adc09', 'J1h_adc10', 'J1h_adc11', 
                                                   'J1h_adc12', 'J1h_adc13', 'J1h_adc15', 
                                                   'J1h_adc16', 'J1h_adc17', 'J1h_adc18',
                                                   'J1h_adc19', 'J1h_adc20', 'J1h_adc21', 
                                                   'J1h_adc22', 'J1h_adc23', 'J1h_adc24', 
                                                   'J1h_adc25'],
                                     positions_microns = [75., 75.+100., 75+400., 
                                                          75.+900., leng/2., (leng/2.) + 200., 
                                                          (leng/2)+700., leng-175., leng-75., 
                                                            900., 1000., 1100., 
                                                            1200., 1300., 1500., 
                                                            1600., 1700., 1800., 
                                                            1900., 2000., 2100., 
                                                            2200., 2300., 2400., 
                                                            leng-50],
                                     initial_profile = profile_J1_init_C)

J1wb = nams.WholeBlock(name = 'Jaipur diopside heated at 904 C for 30m',
                       time_seconds = 30.*60., worksheetname = 'Jaipur 904C 30m',
                       temperature_celsius = 904.,
                       profiles = [profile_J1_904C_30m_A, 
                                   profile_J1_904C_30m_B, 
                                   profile_J1_904C_30m_C])

#%% augite PMR-53 untreated polarized spectra

# polarized IR data digitized from Bell et al. 2004 
PMR_Ea = nams.Spectrum(fname = 'PMR_Ea', polar = 'E || alpha', thick_microns = 1e4)
PMR_Eb = nams.Spectrum(fname = 'PMR_Eb', polar = 'E || beta', thick_microns = 1e4)
PMR_Ec = nams.Spectrum(fname = 'PMR_Ec', polar = 'E || gamma', thick_microns = 1e4)

## The data has since also become available on the PULI database
# spectra from PULI
PULI_PMR53_alpha = nams.Spectrum(fname = "PULI_PMR53_alpha", thick_microns = 1e3) 
PULI_PMR53_beta = nams.Spectrum(fname = "PULI_PMR53_beta", thick_microns = 1e3) 
PULI_PMR53_gamma = nams.Spectrum(fname = "PULI_PMR53_gamma", thick_microns = 1e3) 
# baselines from PULI
PULI_PMR53_alpha_baseline = nams.Spectrum(fname = 'PULI_PMR53_alpha_baseline', thick_microns = 1e3)
PULI_PMR53_beta_baseline = nams.Spectrum(fname = 'PULI_PMR53_beta_baseline', thick_microns = 1e3)
PULI_PMR53_gamma_baseline = nams.Spectrum(fname = 'PULI_PMR53_gamma_baseline', thick_microns = 1e3)

#%% Thin slab of PMR-53 augite dehydrated by time series, unpolarized spectra
PMR_unpol = nams.TimeSeries(fname_list=['P_0_unpol', 'P_1_unpol', 'P_2_unpol',
                                        'P_3_unpol', 'P_4_unpol', 'P_5_unpol', 
                                        'P_6_unpol',],
                            time_hours=[0., 0.25, 0.5, 0.75, 1., 2., 3.],
                            thick_microns=PMR.thickness_microns, 
                            sample=PMR)


#%%################# DIFFUSIVITIES ######################################
### Saved preferred initial values (i), diffusivities (D), and errors (e) 
### for Kunlun diopside K3, K4 at 91 hours of heating (K91) and 154 hours
### of heating (K154), K5, and Jaipur diopside J1 (J)
### In peak wavenumber order: 3645, 3617, 3540, 3443, 3350-->3355, BULK H #####

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

K_dlist = np.array([D_K3, D_K91, D_K154, D_K5])
K_elist = np.array([e_K3, e_K91, e_K154,e_K5])
K_celsius = [816., 904., 904., 1000.]

K_bulk = diffusion.Diffusivities(description='bulk H', 
                                 logD_all=K_dlist[:, -1], 
                                 logD_all_error=K_elist[:, -1], 
                                 celsius_all=K_celsius,
                                 color='k', marker='o', mew=1.5, markersize=10)
K_3645 = diffusion.Diffusivities(description='3645 cm$^{-1}$', 
                                 logD_all=K_dlist[:, 0], 
                                 logD_all_error=K_elist[:, 0],
                                 color='red', marker='x', mew=1.5, markersize=10, 
                                 celsius_all=K_celsius)
K_3617 = diffusion.Diffusivities(description='3617 cm$^{-1}$', 
                                 logD_all=K_dlist[:, 1], 
                                 logD_all_error=K_elist[:, 1],  
                                 color='orange', marker='+', mew=2, markersize=10, 
                                 celsius_all=K_celsius)
K_3540 = diffusion.Diffusivities(description='3540 cm$^{-1}$',
                                 logD_all=K_dlist[:, 2], 
                                 logD_all_error=K_elist[:, 2],
                                 color='teal', marker='s', markersize=8, 
                                 celsius_all=K_celsius)
K_3443 = diffusion.Diffusivities(description='3443 cm$^{-1}$', 
                                 logD_all=K_dlist[:, 3], 
                                 logD_all_error=K_elist[:, 3],  
                                 color='blue', marker='3', mew=2, markersize=10, 
                                 celsius_all=K_celsius)
K_3350 = diffusion.Diffusivities(description='3355 cm$^{-1}$', 
                                 logD_all=K_dlist[:, 4], 
                                 logD_all_error=K_elist[:, 4], alpha=0.5, 
                                 color='indigo', marker='d', markersize=8, 
                                 celsius_all=K_celsius)

J_bulk = nams.diffusion.Diffusivities(description = 'Jaipur bulk H (this work)',
#                                      basestyle = {'color' : 'k', 'marker' : 's', 
#                                                    'markerfacecolor' : 'turquoise',
#                                                    'markersize' :  10, 
#                                                    'linestyle' : 'none', 'alpha' : 1.,
#                                                    'markeredgewidth' : 2},
                                      logDx=[D_J[-1]], logDy=[D_J[-1]-1.], 
                                      logDz=[D_J[-1]], celsius_all=[904.]) 

# augite PMR diffusivities see Fig9_augitePMR.py
D_PMR = nams.diffusion.Diffusivities(celsius_all = [800.],
                                     logD_unoriented = [-11.])
