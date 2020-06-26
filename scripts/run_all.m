%% runs all scripts used in MRM ssEPG manuscript

% add dependencies

addpath('../src')
addpath('../contrib/altmany-export_fig-bb6c842/')
bart_path = ''; % specify path to BART
addpath([bart_path 'matlab']);
setenv('TOOLBOX_PATH', bart_path);
addpath('../contrib/sdc3_nrz_11aug/') % sample density correction


%% scripts for main text 
% (figure save quality is display dependent)
% most scripts produce intermediate figures that are later assembled into
% final figures by figure_gen.m

% ssEPG v. Bloch (Fig. 1; Fig. S1) (Windows)
clear; close all; clc;
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize', 20);
validation_ss_epg_ussfp

% ssEPG, pEPG, and EPG against measured slice profiles (Figs. 2 and 3; Fig.
% S2) (Windows)
clear, close all, clc;
slice_profile_measures

% MR system phantom (Fig. 4) (Windows, Linux, WSL)
clear, close all, clc;
ssepg_optimization_check % (Windows) delta Q / delta N check, results hardcoded in following scripts
batch_proc_dict_MRsys % (Linux) dictionary generation
batch_proc_MRsys % (WSL) MRF image recon and parameter estimation
clear, close all, clc;
analyze_MRsys_data % (Windows)

% B0 effect modeling (Figs. 5 and 6, Table 1) (Windows and Linux)
clear, close all, clc;
batch_proc_dict_b0_examples % (Linux) dictionary generation
clear, close all, clc;
set(0, 'DefaultLineLineWidth', 3);
set(0, 'DefaultAxesFontSize', 28);
b0_effects % (Windows)

% MRF in the calf (Figs. 7 and 8) (Windows, WSL, Linux)
clear, close all, clc;
batch_proc_dict_leg % (Linux) dictionary generation
clear, close all, clc;
batch_proc_leg % (WSL)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize', 20);
analyze_leg % (Windows)

%% Supporting information

% Additional steady-state slice profile and signal modeling (Fig. S3; Table
% S1) (Windows)
clear, close all, clc;
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize', 20);
slice_profile_add_modeling

% Compare pEPG and ssEPG under most optimal TBW and crusher strength in MR
% system phantom (Fig. S4) (Windows)
clear, close all, clc;
compare_pepg_ssepg_MRF_modeling_TBW4C4

% FBIRN phantom (Figs. S5 and S6) (WSL, Windows)
batch_proc_fbirn % (WSL)
clear, close all, clc;
analyze_fbirn % (Windows)

% composite figure generation (Windows)
clear, close all, clc;
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize', 20);
figure_gen




