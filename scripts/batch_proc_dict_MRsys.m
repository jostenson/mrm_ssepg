%% Dictionary generation for MR system phantom measurements

%% dummy section to incure 1 time CUBLAS error (?)

d_temp = gpuArray( rand(3) );
d_temp = pagefun( @mtimes, d_temp, d_temp );

%% ssEPG
% TBW4/C1

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF105.csv';
nomFlip_deg = 60;

shape_in.fov_factor = 4;
shape_in.sl_thick_m = 8/1000;

TRbase_s = 16/1000;
TE_s = 3/1000;
nreps = 1250;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [1.0:0.025:1.35];
input_dict.sfrac = 1 - 1e-4; %

% setup dictionary parameters
data_struct = importdata([dir_in fn_MRF_seq_params]);
input_dict.FA_v = data_struct.data(:,1)*nomFlip_deg; % deg
input_dict.phi_v = data_struct.data(:,2); % deg
input_dict.TR_v = TRbase_s*1000 + data_struct.data(:,3); % ms
input_dict.TE_v = TE_s*1000 + data_struct.data(:,4); % ms
input_dict.nreps = nreps;
input_dict.reduce = 1;

% plot sequence parameters
figure(1); clf;
subplot(411)
plot(input_dict.FA_v); ylabel('degrees'); title(['FA, TI is ' num2str(input_dict.TI) ' ms'])
subplot(412)
plot(input_dict.phi_v); ylabel('degrees'); title('phase')
subplot(413)
plot(input_dict.TR_v); ylabel('msec'); title('TR')
subplot(414)
plot(input_dict.TE_v); ylabel('msec'); title('TE')
drawnow

fn_dict = 'MRF_dict_105_FA60_TBW4_crush1.mat'; % end fn w/crush #
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
shape_in.tau_s = 1.155/1000;

shape_in.n_cycles_per_crush = str2num(fn_dict(end-4));
shape_in.N_t = 128;
shape_in.N_res = 48;

% get rf_shape and gen dict
rf_shape = get_rf_shape( shape_in,nomFlip_deg,1 );   
input_dict.delk = rf_shape.N_cycle;
input_dict.szomega =  rf_shape.Q;
input_dict.rf_shape = rf_shape;
output_dict = MRF_dict_generator_ssEPG_par( input_dict );
save( [dir_out fn_dict], 'output_dict' );

%% ssEPG
% TBW4/C4

fn_dict = 'MRF_dict_105_FA60_TBW4_crush4.mat'; % end fn w/crush #
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
shape_in.tau_s = 1.155/1000;

shape_in.n_cycles_per_crush = str2num(fn_dict(end-4));
shape_in.N_t = 128;
shape_in.N_res = 192;

% get rf_shape and gen dict
rf_shape = get_rf_shape( shape_in,nomFlip_deg,1 );   
input_dict.delk = rf_shape.N_cycle;
input_dict.szomega =  rf_shape.Q;
input_dict.rf_shape = rf_shape;
output_dict = MRF_dict_generator_ssEPG_par( input_dict );
save( [dir_out fn_dict], 'output_dict' );

%% ssEPG
% TBW8/C1

fn_dict = 'MRF_dict_105_FA60_TBW8_crush1.mat'; % end fn w/crush #
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
shape_in.tau_s = 2.325/1000;

shape_in.n_cycles_per_crush = str2num(fn_dict(end-4));
shape_in.N_t = 256;
shape_in.N_res = 64;

% get rf_shape and gen dict
rf_shape = get_rf_shape( shape_in,nomFlip_deg,1 );   
input_dict.delk = rf_shape.N_cycle;
input_dict.szomega =  rf_shape.Q;
input_dict.rf_shape = rf_shape;
output_dict = MRF_dict_generator_ssEPG_par( input_dict );
save( [dir_out fn_dict], 'output_dict' );
 
%% ssEPG
% TBW8/C4

fn_dict = 'MRF_dict_105_FA60_TBW8_crush4.mat'; % end fn w/crush #
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
shape_in.tau_s = 2.325/1000;

shape_in.n_cycles_per_crush = str2num(fn_dict(end-4));
shape_in.N_t = 256;
shape_in.N_res = 256;

% get rf_shape and gen dict
rf_shape = get_rf_shape( shape_in,nomFlip_deg,1 );   
input_dict.delk = rf_shape.N_cycle;
input_dict.szomega =  rf_shape.Q;
input_dict.rf_shape = rf_shape;
output_dict = MRF_dict_generator_ssEPG_par( input_dict );
save( [dir_out fn_dict], 'output_dict' );

%% pEPG
% TBW4

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF105.csv';
nomFlip_deg = 60;

shape_in.fov_factor = 4;
shape_in.sl_thick_m = 8/1000;
shape_in.N_res = round( shape_in.fov_factor * shape_in.sl_thick_m / 50e-6 );

TRbase_s = 16/1000;
TE_s = 3/1000;
nreps = 1250;

% dictionary params
input_dict.np = 50;
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [1.0:0.025:1.35];
input_dict.sfrac = 1 - 1e-4; %

% setup dictionary parameters

data_struct = importdata([dir_in fn_MRF_seq_params]);
input_dict.FA_v = data_struct.data(:,1)*nomFlip_deg; % deg
input_dict.phi_v = data_struct.data(:,2); % deg
input_dict.TR_v = TRbase_s*1000 + data_struct.data(:,3); % ms
input_dict.TE_v = TE_s*1000 + data_struct.data(:,4); % ms
input_dict.nreps = nreps;
input_dict.reduce = 1;

% plot sequence parameters
figure(1); clf;
subplot(411)
plot(input_dict.FA_v); ylabel('degrees'); title(['FA, TI is ' num2str(input_dict.TI) ' ms'])
subplot(412)
plot(input_dict.phi_v); ylabel('degrees'); title('phase')
subplot(413)
plot(input_dict.TR_v); ylabel('msec'); title('TR')
subplot(414)
plot(input_dict.TE_v); ylabel('msec'); title('TE')
drawnow

fn_dict = 'MRF_dict_pepg_105_FA60_TBW4.mat'; % end fn w/crush #
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
shape_in.tau_s = 1.155/1000;
shape_in.N_t = 128;
shape_in.n_cycles_per_crush = 1;

% get rf_shape and gen dict
rf_shape = get_rf_shape( shape_in,nomFlip_deg,1 );   
rf_shape.shape_in = shape_in;
input_dict.delk = 1; 
input_dict.szomega =  101;
input_dict.rf_shape = rf_shape;
output_dict = MRF_dict_generator_pEPG( input_dict );
save( [dir_out fn_dict], 'output_dict' );


%% pEPG
% TBW8

fn_dict = 'MRF_dict_pepg_105_FA60_TBW8.mat'; % end fn w/crush #
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
shape_in.tau_s = 2.325/1000;
shape_in.N_t = 256;
shape_in.n_cycles_per_crush = 1;

% get rf_shape and gen dict
rf_shape = get_rf_shape( shape_in,nomFlip_deg,1 );   
rf_shape.shape_in = shape_in;
input_dict.delk = 1; 
input_dict.szomega =  101;
input_dict.rf_shape = rf_shape;
output_dict = MRF_dict_generator_pEPG( input_dict );
save( [dir_out fn_dict], 'output_dict' );



