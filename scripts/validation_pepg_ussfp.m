%% check pEPG signal against FA <= 30 and crush cycles - 1 > TBW/2

clear, close all, clc;

addpath('../src')
addpath('../contrib/altmany-export_fig-bb6c842/')
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize', 20);
my_dpi = 250;

%% parameters

dir_in = '../data_in/';
dir_out = '../data_out/';
fn_pulse = 'hanning_sinc_ex_tbw4.txt';

% RF params

nomFA_deg = 30;
TBW = 4;

shape_in.fn = [dir_in fn_pulse];
shape_in.N_t = 256;
shape_in.sl_thick_m = 8/1000;
shape_in.tau_s = 3.5/1000;
shape_in.n_cycles_per_crush = 8;
shape_in.fov_factor = 4;
shape_in.N_res = 5000;

% sequeence params
T1_s = 1000/1000;
T2_s = 100/1000;
TR_s = 15/1000;
TE_s = 3/1000;
n_reps = 10;

N_p = 51;

rf_shape = get_rf_shape( shape_in, nomFA_deg, 1 );

%% derivative params

gammabar_Hz_p_T = 42.57747852e6;
tau_s = shape_in.tau_s;
sl_thick_m = shape_in.sl_thick_m;
fov_factor = shape_in.fov_factor;
N_res = shape_in.N_res;
n_cycles_per_crush = shape_in.n_cycles_per_crush;
A_crush_T_s_p_m = n_cycles_per_crush / gammabar_Hz_p_T / sl_thick_m;
grad_T_p_m = rf_shape.grad_T_p_m;
refocus_factor = (rf_shape.diso_idx + rf_shape.N_t/2)/rf_shape.N_t;

%% Bloch simulation

fprintf('Beginning Bloch simulation...\n')

sl_Hz = fov_factor * ( -1/2:1/N_res:1/2 - 1/N_res ) * gammabar_Hz_p_T * grad_T_p_m * sl_thick_m;
t_s = linspace( -tau_s/2,tau_s/2, shape_in.N_t );

% model
M0 = [0 0 1]';

M_solved = zeros( N_res,3 );
tic;
for ii = 1:N_res
    [~,M] = ode45( @(t,m) bloch_eq(t,m,rf_shape.B1_T,sl_Hz(ii),t_s),[min(t_s) max(t_s)],M0);
    M_solved(ii,:) = M(end,:);
end
fprintf('Time for Bloch simulation %.1f ms\n', 1000*toc);

z_m = sl_Hz / (gammabar_Hz_p_T * grad_T_p_m);
Mxy_bloch = M_solved(:,1) + 1i*M_solved(:,2);
Mz_bloch = M_solved(:,3);

Mxy_bloch = Mxy_bloch.* exp( 1i * 2 * pi * gammabar_Hz_p_T * grad_T_p_m * z_m(:) * (tau_s * refocus_factor) ); % refocus

figure(2); clf;
plot( z_m*1e3, real(Mxy_bloch),'--' );
hold on
plot( z_m*1e3, imag(Mxy_bloch),'--' );
plot( z_m*1e3, abs(Mxy_bloch) )
hold off
xlim([-fov_factor/2*sl_thick_m fov_factor/2*sl_thick_m]*1000)
xlabel('slice position (mm)','FontSize',16); ylabel('M_{x,y}','FontSize',16);
legend('M_x','M_y','|M_{xy}|')
title('Bloch simulation')

figure(3); clf;
plot( z_m*1e3, Mz_bloch )
xlim([-fov_factor/2*sl_thick_m fov_factor/2*sl_thick_m]*1000)
xlabel('slice position (mm)','FontSize',16); ylabel('M_{z}','FontSize',16);
title('Bloch simulation')
drawnow

bloch_prof_mag = abs(sum(Mxy_bloch));
fprintf( 'Bloch simulation mag of sum of complex profile %.3f\n', bloch_prof_mag  )

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

fprintf( 'EPG simulation mag of DC component scaled to match Bloch sim resolution %.3f\n', epg_prof_mag  )
fprintf( 'Ratio of EPG to Bloch profile magnitudes %.3f\n', epg_prof_mag/bloch_prof_mag );

%% model time series with pEPG


TE_v = TE_s*ones(n_reps,1);
TR_v = TR_s*ones(n_reps,1);
FA_v = nomFA_deg*ones(n_reps,1);
phi_v = zeros(n_reps,1);
TI = inf;
delk = 1;
szomega = 101;

input_pepg.T1 = T1_s;
input_pepg.T2 = T2_s;
input_pepg.TE_v = TE_v;
input_pepg.TR_v = TR_v;
input_pepg.FA_v = FA_v;
input_pepg.delk = delk;
input_pepg.n_reps= n_reps;
input_pepg.szomega = szomega;
input_pepg.phi_v = phi_v;
input_pepg.TI = TI;

input_pepg.N_p = N_p;

% get profile
Mxy_crop = Mxy( z_epg_m >= -fov_factor/2*sl_thick_m & z_epg_m <= fov_factor/2*sl_thick_m );
Mxy_crop = Mxy_crop*numel( Mxy );

% get FA multiplier for all possible partitions
fa_slice_factors = asind( imag(Mxy_crop) )/nomFA_deg;

% determine partiton spacing
idx_partition = round(1:(numel(fa_slice_factors)/N_p):numel(fa_slice_factors));
input_pepg.fa_factors = fa_slice_factors( idx_partition );

[sig_pepg_v,sig_pepg_full] = pEPG_SSFP( input_pepg );
sig_pepg_v = -sig_pepg_v; % change sign convention


%% model time series with ssEPG

delk = rf_shape.N_cycle;
szomega =  Q;

tic;
[Fxy,sig_epg_v] = EPG_MRF_SSFP_prof( T1_s, T2_s, TE_v, TR_v, FA_v, delk, n_reps,szomega, phi_v, TI,rf_shape );
fprintf('Time for EPG time series %.1f ms\n', 1000*toc);
sig_epg_prof = sig_epg_v * (2*Q - 2);


idx_bw_match = numel( z_epg_m )/2 - shape_in.N_res/2:numel( z_epg_m)/2 + shape_in.N_res/2 - 1;

sig_epg_prof_tot = ifftshift( ifft( fftshift( Fxy, 2 ), [], 2 ), 2 );

sig_epg_prof_tot_crop = zeros( n_reps, N_res );

%% stats

sig_pepg_norm = sig_pepg_v;
sig_ssepg_norm = sig_epg_prof/norm(sig_epg_prof);

resid_pepg_norm = sig_pepg_norm(:) - sig_ssepg_norm(:);

RMSE_pepg_norm = sqrt( resid_pepg_norm'*resid_pepg_norm/ n_reps );

fprintf( 'RMSE of pEPG normed %.3f\n', RMSE_pepg_norm );

figure(23); clf;
set(gcf,'Color','w','WindowState','maximize')
plot(abs(sig_ssepg_norm))
hold on
plot(abs(sig_pepg_norm),'--')
hold off
grid on;
xlabel('excitation'); ylabel('Normed signal (au)');
legend('ssEPG','pEPG')
t1 = text(5,0.45,sprintf( 'RMSE of pEPG = %.3f',RMSE_pepg_norm ),'FontSize',28);
% export_fig('../figures/?.png');


