%% check pEPG dictionary against pEPG function

clear, close all, clc;

addpath('../src')
addpath('../contrib/altmany-export_fig-bb6c842/')
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize', 20);
my_dpi = 250;

dir_in = '../data_in/';
dir_out = '../data_out/';
fn_pulse = 'hanning_sinc_ex_tbw4.txt';

%% params

% RF params

nomFA_deg = 60;
TBW = 4;

shape_in.fn = [dir_in fn_pulse];
shape_in.N_t = 256;
shape_in.sl_thick_m = 8/1000;
shape_in.tau_s = 3.5/1000;
shape_in.n_cycles_per_crush = 8;
shape_in.fov_factor = 4;
shape_in.N_res = 5000;

% sequeence params
fn_MRF_seq_params = 'MRF104.csv';
T1_v = [300 1000 2000]/1000;
T2_v = [30 50 100]/1000;
TR_s = 15/1000;
TE_s = 3/1000;
n_reps = 1000;

N_p = 51;

rf_shape = get_rf_shape( shape_in, nomFA_deg, 1 );

%% derivative params

gammabar_Hz_p_T = 42.57747852e6;
sl_thick_m = shape_in.sl_thick_m;
fov_factor = shape_in.fov_factor;
grad_T_p_m = rf_shape.grad_T_p_m;

% load seq param file
data_struct = importdata([dir_in fn_MRF_seq_params]);

% set dictionary sequence parameters
FA_v = data_struct.data(:,1)*nomFA_deg; % deg
phi_v = data_struct.data(:,2); % deg
TR_v = TR_s + data_struct.data(:,3)/1000; % s
TE_v = TE_s + data_struct.data(:,4)/1000; % s
TI = inf;
delk = 1;
szomega = 101;

n_dict = numel(T1_v)*numel(T2_v);

%% ssEPG

Q = rf_shape.Q;
dtau_s = rf_shape.dtau_s;

omega0 = zeros( 3, Q );
omega0(:,1) = [0 0 1].';

tic; 
omega = ssEPG_RF( dtau_s,rf_shape.B1_T,omega0,0,rf_shape.diso_idx );
fprintf('Time for EPG simulation %.1f ms\n', 1000*toc);

epg_prof_mag = abs(omega(1,1)) * (2*Q - 2);

dkz_per_m = dtau_s * gammabar_Hz_p_T * grad_T_p_m;
z_epg_fov_m = 1/dkz_per_m;
dz_epg_m = z_epg_fov_m / (2*Q-2);
z_epg_m = -z_epg_fov_m/2:dz_epg_m:z_epg_fov_m/2 - dz_epg_m;

Fxy = [fliplr(conj(omega(2,2:end))) omega(1,1:end-1) ];
Fz = [fliplr(conj(omega(3,2:end))) omega(3,1:end-1) ];

Mxy = ifftshift( ifft( fftshift( Fxy ) ) );
Mz = ifftshift( ifft( fftshift( Fz ) ) );

figure(5); clf;
plot(z_epg_m*1000,real(Mxy*numel(Mxy)),'--')
hold on
plot(z_epg_m*1000,imag(Mxy*numel(Mxy)),'--')
plot(z_epg_m*1000,abs(Mxy*numel(Mxy)))
hold off
xlabel('slice position (mm)','FontSize',16); ylabel('M_{x,y}','FontSize',16);
xlim([-fov_factor/2*sl_thick_m fov_factor/2*sl_thick_m]*1000)
legend('M_x','M_y','|M_{xy}|')
title('EPG simulation')

figure(6); clf;
plot(z_epg_m*1000,(Mz*numel(Mz)))
xlabel('slice position (mm)','FontSize',16); ylabel('M_{z}','FontSize',16);
xlim([-fov_factor/2*sl_thick_m fov_factor/2*sl_thick_m]*1000)
title('EPG simulation')
drawnow

%% model time series with pEPG fun

input_fun.TE_v = TE_v;
input_fun.TR_v = TR_v;
input_fun.FA_v = FA_v;
input_fun.delk = delk;
input_fun.n_reps= n_reps;
input_fun.szomega = szomega;
input_fun.phi_v = phi_v;
input_fun.TI = TI;

input_fun.N_p = N_p;

% get profile
Mxy_crop = Mxy( z_epg_m >= -fov_factor/2*sl_thick_m & z_epg_m <= fov_factor/2*sl_thick_m );
Mxy_crop = Mxy_crop*numel( Mxy );

% get FA multiplier for all possible partitions
fa_slice_factors = asind( imag(Mxy_crop) )/nomFA_deg;

% determine partiton spacing
idx_partition = round(1:(numel(fa_slice_factors)/N_p):numel(fa_slice_factors));
input_fun.fa_factors = fa_slice_factors( idx_partition );

pepg_dict_fun = zeros( n_reps, n_dict );
nn = 1;
for ii = 1:numel(T2_v)
    input_fun.T2 = T2_v(ii);
    for jj = 1:numel(T1_v)

        input_fun.T1 = T1_v(jj);
        
        [sig_pepg_v,~] = pEPG_SSFP( input_fun );
        %sig_pepg_v = -sig_pepg_v; % change sign convention
        pepg_dict_fun(:,nn) = sig_pepg_v(:);
        nn = nn + 1;
        
    end    
end

%% model time series with pepg dictionary generator

input_dict = input_fun;
input_dict.T1_v = T1_v;
input_dict.T2_v = T2_v;
input_dict.T1 = [];
input_dict.T2 = [];
input_dict.B1_v = [1];
input_dict.reduce = 0;
input_dict.sfrac = 1;
input_dict.rf_shape = rf_shape;
input_dict.np = N_p;
input_dict.nreps = n_reps;

output_dict = MRF_dict_generator_pEPG( input_dict );
pepg_dict_gen = output_dict.dict_norm;

%% plot results

% example signal magnitudes
figure(10); clf;
plot( abs( pepg_dict_fun(:,1:end) ) );
hold on
plot( abs( pepg_dict_gen(:,1:end) ),'--' );
hold off
xlabel('Excitation')
ylabel('|signal|')
legend('serial function','parallel generator')

figure(11); clf;
plot( log10( abs( pepg_dict_gen - pepg_dict_fun ) ) );
xlabel('Excitation')
ylabel('log10|signal|')
title('difference: generator from serial fun')