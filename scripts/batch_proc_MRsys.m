%% Batch MRF image processing and parameter estimation for MR system phantom

%% TBW4/C1

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_MRsys.mat';
fn_MRF_raw = 'MRF105_TBW4_FA60_crush1.mat';
fn_B1 = 'B1_MRsys.mat';
fn_MRF_water_dict_epg = 'MRF_dict_105_FA60_epg.mat';
fn_MRF_water_dict_pepg = 'MRF_dict_pepg_105_FA60_TBW4.mat';
fn_MRF_water_dict_ssepg = 'MRF_dict_105_FA60_TBW4_crush1.mat';
fn_MRF_seq_params = 'MRF105.csv';
fn_MRF_img_recon = 'MRF105_TBW4_FA60_crush1_img_recon.mat';
fn_MRF_proc_no_fat_sep_epg = 'MRF105_TBW4_FA60_crush1_proc_epg.mat';
fn_MRF_proc_no_fat_sep_pepg = 'MRF105_TBW4_FA60_crush1_proc_pepg.mat';
fn_MRF_proc_no_fat_sep_ssepg = 'MRF105_TBW4_FA60_crush1_proc_ssepg.mat';

skiprecon = 0;
doMagdict_epg = 1; % if 1, create dictionary
doMagdict_pepg = 0; %see batch_proc_dict.m
doMagdict_ssepg = 0; %see batch_proc_dict.m
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [1.0:0.025:1.35]; %
input_dict.sfrac = 1 - 1e-4; % 

MRF_master_proc


%% TBW4/C4

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_MRsys.mat';
fn_MRF_raw = 'MRF105_TBW4_FA60_crush4.mat';
fn_B1 = 'B1_MRsys.mat';
fn_MRF_water_dict_epg = 'MRF_dict_105_FA60_epg.mat';
fn_MRF_water_dict_pepg = 'MRF_dict_pepg_105_FA60_TBW4.mat';
fn_MRF_water_dict_ssepg = 'MRF_dict_105_FA60_TBW4_crush4.mat';
fn_MRF_seq_params = 'MRF105.csv';
fn_MRF_img_recon = 'MRF105_TBW4_FA60_crush4_img_recon.mat';
fn_MRF_proc_no_fat_sep_epg = 'MRF105_TBW4_FA60_crush4_proc_epg.mat';
fn_MRF_proc_no_fat_sep_pepg = 'MRF105_TBW4_FA60_crush4_proc_pepg.mat';
fn_MRF_proc_no_fat_sep_ssepg = 'MRF105_TBW4_FA60_crush4_proc_ssepg.mat';

skiprecon = 0;
doMagdict_epg = 0; % if 1, create dictionary
doMagdict_pepg = 0; %see batch_proc_dict.m
doMagdict_ssepg = 0; %see batch_proc_dict.m
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [1.0:0.025:1.35]; %
input_dict.sfrac = 1 - 1e-4; %

MRF_master_proc

%% TBW8/C1

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_MRsys.mat';
fn_MRF_raw = 'MRF105_TBW8_FA60_crush1.mat';
fn_B1 = 'B1_MRsys.mat';
fn_MRF_water_dict_epg = 'MRF_dict_105_FA60_epg.mat';
fn_MRF_water_dict_pepg = 'MRF_dict_pepg_105_FA60_TBW8.mat';
fn_MRF_water_dict_ssepg = 'MRF_dict_105_FA60_TBW8_crush1.mat';
fn_MRF_seq_params = 'MRF105.csv';
fn_MRF_img_recon = 'MRF105_TBW8_FA60_crush1_img_recon.mat';
fn_MRF_proc_no_fat_sep_epg = 'MRF105_TBW8_FA60_crush1_proc_epg.mat';
fn_MRF_proc_no_fat_sep_pepg = 'MRF105_TBW8_FA60_crush1_proc_pepg.mat';
fn_MRF_proc_no_fat_sep_ssepg = 'MRF105_TBW8_FA60_crush1_proc_ssepg.mat';

skiprecon = 0;
doMagdict_epg = 0; % if 1, create dictionary
doMagdict_pepg = 0; %see batch_proc_dict.m
doMagdict_ssepg = 0; %see batch_proc_dict.m
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [1.0:0.025:1.35]; %
input_dict.sfrac = 1 - 1e-4; %

MRF_master_proc

%% TBW8/C4

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_MRsys.mat';
fn_MRF_raw = 'MRF105_TBW8_FA60_crush4.mat';
fn_B1 = 'B1_MRsys.mat';
fn_MRF_water_dict_epg = 'MRF_dict_105_FA60_epg.mat';
fn_MRF_water_dict_pepg = 'MRF_dict_pepg_105_FA60_TBW8.mat';
fn_MRF_water_dict_ssepg = 'MRF_dict_105_FA60_TBW8_crush4.mat';
fn_MRF_seq_params = 'MRF105.csv';
fn_MRF_img_recon = 'MRF105_TBW8_FA60_crush4_img_recon.mat';
fn_MRF_proc_no_fat_sep_epg = 'MRF105_TBW8_FA60_crush4_proc_epg.mat';
fn_MRF_proc_no_fat_sep_pepg = 'MRF105_TBW8_FA60_crush4_proc_pepg.mat';
fn_MRF_proc_no_fat_sep_ssepg = 'MRF105_TBW8_FA60_crush4_proc_ssepg.mat';

skiprecon = 0;
doMagdict_epg = 0; % if 1, create dictionary
doMagdict_pepg = 0; %see batch_proc_dict.m
doMagdict_ssepg = 0; %see batch_proc_dict.m
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [1.0:0.025:1.35]; %
input_dict.sfrac = 1 - 1e-4; %

MRF_master_proc






