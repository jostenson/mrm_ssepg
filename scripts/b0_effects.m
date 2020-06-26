%% B0 effect modeling

my_dpi = 250;

%% TBW4

dir_in = '../data_in/';
dir_out = '../data_out/';

T1_ms = 1320;
T2_ms = 30;

fn_MRF_seq_params = 'MRF105.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 1;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

B0_kHz = 0;
[Fxy_0,Ft_0] = EPG_MRF_SSFP_prof_B0( T1_ms,T2_ms,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,B0_kHz,i_sign );

B0_kHz_1 = 1/TRbase_s/1000/4;
[Fxy_B0,Ft_B01] = EPG_MRF_SSFP_prof_B0( T1_ms,T2_ms,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,B0_kHz_1,i_sign );

B0_kHz_2 = 1/TRbase_s/1000/2;
[Fxy_B02,Ft_B02] = EPG_MRF_SSFP_prof_B0( T1_ms,T2_ms,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,B0_kHz_2,i_sign );

B0_kHz_3 = 5/TRbase_s/1000/4;
[Fxy_B03,Ft_B03] = EPG_MRF_SSFP_prof_B0( T1_ms,T2_ms,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,B0_kHz_3,i_sign );

B0_kHz_4 = 3/TRbase_s/1000/2;
[Fxy_B04,Ft_B04] = EPG_MRF_SSFP_prof_B0( T1_ms,T2_ms,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,B0_kHz_4,i_sign );

Ft_B05 = Ft_B01(:).* exp( 1i * 2 * pi * TE_ms(:)/TR_ms(1) );
Ft_B06 = Ft_B02(:).* exp( 1i * 2 * pi * TE_ms(:)/TR_ms(1) );

%% plot results
% {'-',':','-.','--','*','+'}

figure(1); clf;
set(gcf,'Position',[200 100 1100 700],'Color','w')
plot( abs( Ft_0./norm(Ft_0) ) );
hold on
plot( abs( Ft_B01./norm(Ft_B01) ), '-.' );
plot( abs( Ft_B02./norm(Ft_B02) ), '-.' );
plot( abs( Ft_B03./norm(Ft_B03) ), ':' );
plot( abs( Ft_B04./norm(Ft_B04) ), ':' );
hold off
xlabel('Excitation #');
ylabel('Signal mag. (au)');
lgd = legend('0','1/4','1/2','5/4','3/2');
legend boxoff;
set( lgd, 'interpreter','latex');
title(lgd,'$\Delta B_0$ [$\frac{1}{TR}]$','interpreter','latex')
export_fig('../figures/mag_with_B0.png', sprintf('-r%d',my_dpi));

figure(2); clf;
set(gcf,'Position',[200 100 1100 700],'Color','w')
plot( angle( Ft_0./norm(Ft_0) ) );
hold on
plot( angle( Ft_B01./norm(Ft_B01) ), '-.' );
plot( angle( Ft_B02./norm(Ft_B02) ), '-.' );
plot( angle( Ft_B03./norm(Ft_B03) ), ':' );
plot( angle( Ft_B04./norm(Ft_B04) ), ':' );
plot( angle( Ft_B05./norm(Ft_B05) ), '--' );
plot( angle( Ft_B06./norm(Ft_B06) ), '--' );
hold off
xlabel('Excitation #');
ylabel('Signal phase (rad)');
lgd = legend('0','1/4','1/2','5/4','3/2','5/4*','3/2*');
legend boxoff;
set( lgd, 'interpreter','latex');
title(lgd,'$\Delta B_0$ [$\frac{1}{TR}]$','interpreter','latex')
export_fig('../figures/phase_with_B0.png', sprintf('-r%d',my_dpi));

%% examine slice profiles

Mxy_0 = fftshift(ifft(fftshift(Fxy_0,2),[],2),2);
n_full = size( Mxy_0,2 );
Mxy_0 = Mxy_0(:,n_full/2+1-shape_in.N_res/2:n_full/2+shape_in.N_res/2)*n_full;

Mxy_B01 = fftshift(ifft(fftshift(Fxy_B0,2),[],2),2);
Mxy_B01 = Mxy_B01(:,n_full/2+1-shape_in.N_res/2:n_full/2+shape_in.N_res/2)*n_full;

Mxy_B02 = fftshift(ifft(fftshift(Fxy_B02,2),[],2),2);
Mxy_B02 = Mxy_B02(:,n_full/2+1-shape_in.N_res/2:n_full/2+shape_in.N_res/2)*n_full;

Mxy_B03 = fftshift(ifft(fftshift(Fxy_B03,2),[],2),2);
Mxy_B03 = Mxy_B03(:,n_full/2+1-shape_in.N_res/2:n_full/2+shape_in.N_res/2)*n_full;

Mxy_B04 = fftshift(ifft(fftshift(Fxy_B04,2),[],2),2);
Mxy_B04 = Mxy_B04(:,n_full/2+1-shape_in.N_res/2:n_full/2+shape_in.N_res/2)*n_full;

n_slice = shape_in.N_res/shape_in.fov_factor;
x_slice = ( (1:shape_in.N_res) - (shape_in.N_res/2 + 1) )/n_slice;

figure(20); clf;
set(gcf,'Position',[200 100 1100 700],'Color','w')
plot( x_slice, abs(Mxy_0(2,:)) );
hold on
plot( x_slice, abs(Mxy_B01(2,:)), '-.' );
plot( x_slice, abs(Mxy_B02(2,:)), '-.' );
plot( x_slice, abs(Mxy_B03(2,:)), ':' );
plot( x_slice, abs(Mxy_B04(2,:)), ':' );
hold off
xlabel('Nominal slice thickness')
ylabel('Signal mag. (au)')
lgd = legend('0','1/4','1/2','5/4','3/2');
legend boxoff;
set( lgd, 'interpreter','latex');
title(lgd,'$\Delta B_0$ [$\frac{1}{TR}]$','interpreter','latex')
export_fig('../figures/profiles_with_B0.png', sprintf('-r%d',my_dpi));


%% T1 and T2 bias estimation MRF105 TBW4/C1
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

k1 = parallel.gpu.CUDAKernel('ssepg_functions.ptx','ssepg_functions.cu','dephase_gradients_rf_stage1');
k2 = parallel.gpu.CUDAKernel('ssepg_functions.ptx','ssepg_functions.cu','dephase_gradients_rf_stage2');
kerns.k1 = k1; kerns.k2 = k2;
k1.ThreadBlockSize = [128 1 1];
k1.GridSize = [64 1 1];
k2.ThreadBlockSize = [128 1 1];
k2.GridSize = [64 1 1];

fn_MRF_seq_params = 'MRF105.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

T1_list = repmat( [1350, 820, 1350, 820], [1 5] );
T2_list = repmat( [30, 30, 85, 65], [1 5] );
B0_list = repmat( [ 0 0.25 0.5 5/4 3/2]*1/TRbase_s/1000, [4 1] );
B0_list = B0_list(:)';

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 1;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW4C1] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF105 TBW4/C4
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF105.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 4;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW4C4] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF105 TBW4/C8
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF105.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 8;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW4C8] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF105 TBW8/C1
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF105.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 2.325/1000;
shape_in.n_cycles_per_crush = 1;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 256;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW8C1] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF105 TBW8/C4
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF105.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 2.325/1000;
shape_in.n_cycles_per_crush = 4;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 256;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW8C4] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF105 TBW8/C8
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF105.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 2.325/1000;
shape_in.n_cycles_per_crush = 8;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 256;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW8C8] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

Ft_MRF105 = [ Ft_TBW4C1 Ft_TBW4C4 Ft_TBW8C1 Ft_TBW8C4 ];
Ft_MRF105_C8 = [ Ft_TBW4C8 Ft_TBW8C8 ];


%% T1 and T2 bias estimation MRF104 TBW4/C1
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF104.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

T1_list = repmat( [1350, 820, 1350, 820], [1 5] );
T2_list = repmat( [30, 30, 85, 65], [1 5] );
B0_list = repmat( [ 0 0.25 0.5 5/4 3/2]*1/TRbase_s/1000, [4 1] );
B0_list = B0_list(:)';

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 1;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW4C1] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF104 TBW4/C4
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF104.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 4;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW4C4] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF104 TBW4/C8
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF104.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 8;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW4C8] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF104 TBW8/C1
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF104.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 2.325/1000;
shape_in.n_cycles_per_crush = 1;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 256;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW8C1] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF104 TBW8/C4
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF104.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 2.325/1000;
shape_in.n_cycles_per_crush = 4;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 256;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW8C4] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF104 TBW8/C8
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF104.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 2.325/1000;
shape_in.n_cycles_per_crush = 8;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 256;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1250;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW8C8] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

Ft_MRF104 = [ Ft_TBW4C1 Ft_TBW4C4 Ft_TBW8C1 Ft_TBW8C4 ];
Ft_MRF104_C8 = [ Ft_TBW4C8 Ft_TBW8C8 ];


%% T1 and T2 bias estimation MRF001 TBW4/C1
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

fn_MRF_seq_params = 'MRF001.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

T1_list = repmat( [1350, 820, 1350, 820], [1 5] );
T2_list = repmat( [30, 30, 85, 65], [1 5] );
B0_list = repmat( [ 0 0.25 0.5 5/4 3/2]*1/TRbase_s/1000, [4 1] );
B0_list = B0_list(:)';

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 1;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1000;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW4C1] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF001 TBW4/C4
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

fn_MRF_seq_params = 'MRF001.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 4;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1000;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW4C4] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF001 TBW4/C8
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

fn_MRF_seq_params = 'MRF001.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw4.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 1.155/1000;
shape_in.n_cycles_per_crush = 8;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 128;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1000;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW4C8] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF001 TBW8/C1
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF001.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 2.325/1000;
shape_in.n_cycles_per_crush = 1;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 256;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1000;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW8C1] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF001 TBW8/C4
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF001.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 2.325/1000;
shape_in.n_cycles_per_crush = 4;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 256;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1000;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW8C4] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

%% T1 and T2 bias estimation MRF001 TBW8/C8
% TBW 4/8, Crush 1/4, T1/T2 1350/30, 820/30, 1350/85, 820/65

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF_seq_params = 'MRF001.csv';
shape_in.fn = [dir_in 'hanning_sinc_ex_tbw8.txt']; % end fn w/tbw
nomFlip_deg = 60;
TRbase_s = 16/1000;

shape_in.tau_s = 2.325/1000;
shape_in.n_cycles_per_crush = 8;
shape_in.fov_factor = 4;
shape_in.sl_thick_m = 5/1000;
shape_in.N_t = 256;
shape_in.N_res = 1024;

TE_s = 3/1000;
nreps = 1000;

% sequence params
TI_ms = 40.0; % ms
data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_deg = data_struct.data(:,1)*nomFlip_deg; % deg
phi_deg = data_struct.data(:,2); % deg
TR_ms = TRbase_s*1000 + data_struct.data(:,3); % ms
TE_ms = TE_s*1000 + data_struct.data(:,4); % ms

% plot sequence parameters
figure(9); clf;
subplot(411)
plot(FA_deg(1:nreps)); ylabel('degrees'); title(['FA, TI is ' num2str(TI_ms) ' ms'])
subplot(412)
plot(phi_deg(1:nreps)); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_ms(1:nreps)); ylabel('msec'); title('TR')
subplot(414)
plot(TE_ms(1:nreps)); ylabel('msec'); title('TE')
drawnow

% get rf_shape for ssEPG
i_sign = 1;
rf_shape = get_rf_shape(shape_in,nomFlip_deg,1);   
delk = rf_shape.N_cycle;
szomega =  rf_shape.Q;

[~,Ft_TBW8C8] = EPG_MRF_SSFP_B0_prof_par( T1_list, T2_list, B0_list,TE_ms,TR_ms,FA_deg,delk,nreps,szomega,phi_deg,TI_ms,rf_shape,numel(T1_list),kerns,0);

Ft_MRF001 = [ Ft_TBW4C1 Ft_TBW4C4 Ft_TBW8C1 Ft_TBW8C4 ];
Ft_MRF001_C8 = [ Ft_TBW4C8 Ft_TBW8C8 ];

save( [dir_out 'b0_mrf_signals.mat'],'Ft_MRF105','Ft_MRF104','Ft_MRF001','Ft_MRF105_C8','Ft_MRF104_C8','Ft_MRF001_C8' );

%% match MRF signals with prof and B0 effects against dictionaries with prof without B0 effects

n_set = numel( T1_list );

Ft_MRF105 = [Ft_MRF105(:,1:n_set*2) Ft_MRF105_C8(:,1:n_set) Ft_MRF105(:,n_set*2+1:end) Ft_MRF105_C8(:,n_set+1:end)];
Ft_MRF104 = [Ft_MRF104(:,1:n_set*2) Ft_MRF104_C8(:,1:n_set) Ft_MRF104(:,n_set*2+1:end) Ft_MRF104_C8(:,n_set+1:end)];
Ft_MRF001 = [Ft_MRF001(:,1:n_set*2) Ft_MRF001_C8(:,1:n_set) Ft_MRF001(:,n_set*2+1:end) Ft_MRF001_C8(:,n_set+1:end)];

fn_d_105_tbw4c1 = 'MRF_105_dict_example_TBW4_crush1.mat';
fn_d_105_tbw4c4 = 'MRF_105_dict_example_TBW4_crush4.mat';
fn_d_105_tbw4c8 = 'MRF_105_dict_example_TBW4_crush8.mat';
fn_d_105_tbw8c1 = 'MRF_105_dict_example_TBW8_crush1.mat';
fn_d_105_tbw8c4 = 'MRF_105_dict_example_TBW8_crush4.mat';
fn_d_105_tbw8c8 = 'MRF_105_dict_example_TBW8_crush8.mat';

fn_d_104_tbw4c1 = 'MRF_104_dict_example_TBW4_crush1.mat';
fn_d_104_tbw4c4 = 'MRF_104_dict_example_TBW4_crush4.mat';
fn_d_104_tbw4c8 = 'MRF_104_dict_example_TBW4_crush8.mat';
fn_d_104_tbw8c1 = 'MRF_104_dict_example_TBW8_crush1.mat';
fn_d_104_tbw8c4 = 'MRF_104_dict_example_TBW8_crush4.mat';
fn_d_104_tbw8c8 = 'MRF_104_dict_example_TBW8_crush8.mat';

fn_d_001_tbw4c1 = 'MRF_001_dict_example_TBW4_crush1.mat';
fn_d_001_tbw4c4 = 'MRF_001_dict_example_TBW4_crush4.mat';
fn_d_001_tbw4c8 = 'MRF_001_dict_example_TBW4_crush8.mat';
fn_d_001_tbw8c1 = 'MRF_001_dict_example_TBW8_crush1.mat';
fn_d_001_tbw8c4 = 'MRF_001_dict_example_TBW8_crush4.mat';
fn_d_001_tbw8c8 = 'MRF_001_dict_example_TBW8_crush8.mat';

%

load( [dir_out fn_d_105_tbw4c1] );
list_4_1 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_4_1 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_105_tbw4c4] );
list_4_4 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_4_4 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_105_tbw4c8] );
list_4_8 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_4_8 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_105_tbw8c1] );
list_8_1 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_8_1 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_105_tbw8c4] );
list_8_4 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_8_4 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_105_tbw8c8] );
list_8_8 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_8_8 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

%

T1_fits_105 = zeros( n_set, 6 );
T2_fits_105 = T1_fits_105;
for ii = 1:6
  
    my_idx = (ii-1)*n_set + 1:ii*n_set;
    my_sigs = Ft_MRF105(:,my_idx);
    
    switch ii
        case 1
            my_dict = dict_4_1;
            my_list = list_4_1;
        case 2
            my_dict = dict_4_4;
            my_list = list_4_4;
        case 3
            my_dict = dict_4_8;
            my_list = list_4_8;
        case 4
            my_dict = dict_8_1;
            my_list = list_8_1;
        case 5
            my_dict = dict_8_4;
            my_list = list_8_4;
        case 6
            my_dict = dict_8_8;
            my_list = list_8_8;
    end
    
    ip = my_sigs' * my_dict;
    [~,idx_max] = max( ip, [], 2 );
    T1_fits_105( :,ii ) = my_list( idx_max,1 );
    T2_fits_105( :,ii ) = my_list( idx_max,2 );
    
end

T1_bias_105 = T1_fits_105 - repmat( T1_list(:), [1 6] );
T2_bias_105 = T2_fits_105 - repmat( T2_list(:), [1 6] );
T1_rel_bias_105 = T1_bias_105./repmat( T1_list(:), [1 6] );
T2_rel_bias_105 = T2_bias_105./repmat( T2_list(:), [1 6] );

%

load( [dir_out fn_d_104_tbw4c1] );
list_4_1 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_4_1 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_104_tbw4c4] );
list_4_4 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_4_4 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_104_tbw4c8] );
list_4_8 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_4_8 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_104_tbw8c1] );
list_8_1 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_8_1 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_104_tbw8c4] );
list_8_4 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_8_4 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_104_tbw8c8] );
list_8_8 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_8_8 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

T1_fits_104 = zeros( n_set, 6 );
T2_fits_104 = T1_fits_104;
for ii = 1:6
  
    my_idx = (ii-1)*n_set + 1:ii*n_set;
    my_sigs = Ft_MRF104(:,my_idx);
    
    switch ii
        case 1
            my_dict = dict_4_1;
            my_list = list_4_1;
        case 2
            my_dict = dict_4_4;
            my_list = list_4_4;
        case 3
            my_dict = dict_4_8;
            my_list = list_4_8;
        case 4
            my_dict = dict_8_1;
            my_list = list_8_1;
        case 5
            my_dict = dict_8_4;
            my_list = list_8_4;
        case 6
            my_dict = dict_8_8;
            my_list = list_8_8;
    end
    
    ip = my_sigs' * my_dict;
    [~,idx_max] = max( ip, [], 2 );
    T1_fits_104( :,ii ) = my_list( idx_max,1 );
    T2_fits_104( :,ii ) = my_list( idx_max,2 );
    
end

T1_bias_104 = T1_fits_104 - repmat( T1_list(:), [1 6] );
T2_bias_104 = T2_fits_104 - repmat( T2_list(:), [1 6] );
T1_rel_bias_104 = T1_bias_104./repmat( T1_list(:), [1 6] );
T2_rel_bias_104 = T2_bias_104./repmat( T2_list(:), [1 6] );

%

load( [dir_out fn_d_001_tbw4c1] );
list_4_1 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_4_1 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_001_tbw4c4] );
list_4_4 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_4_4 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_001_tbw4c8] );
list_4_8 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_4_8 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_001_tbw8c1] );
list_8_1 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_8_1 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_001_tbw8c4] );
list_8_4 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_8_4 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

load( [dir_out fn_d_001_tbw8c8] );
list_8_8 = output_dict.dict_list( abs( output_dict.dict_list(:,3) - 1 ) < eps, : );
dict_8_8 = output_dict.dict_norm( :, abs( output_dict.dict_list(:,3) - 1 ) < eps );

T1_fits_001 = zeros( n_set, 6 );
T2_fits_001 = T1_fits_001;
for ii = 1:6
  
    my_idx = (ii-1)*n_set + 1:ii*n_set;
    my_sigs = Ft_MRF001(:,my_idx);
    
    switch ii
        case 1
            my_dict = dict_4_1;
            my_list = list_4_1;
        case 2
            my_dict = dict_4_4;
            my_list = list_4_4;
        case 3
            my_dict = dict_4_8;
            my_list = list_4_8;
        case 4
            my_dict = dict_8_1;
            my_list = list_8_1;
        case 5
            my_dict = dict_8_4;
            my_list = list_8_4;
        case 6
            my_dict = dict_8_8;
            my_list = list_8_8;
    end
    
    ip = my_sigs' * my_dict;
    [~,idx_max] = max( ip, [], 2 );
    T1_fits_001( :,ii ) = my_list( idx_max,1 );
    T2_fits_001( :,ii ) = my_list( idx_max,2 );
    
end

T1_bias_001 = T1_fits_001 - repmat( T1_list(:), [1 6] );
T2_bias_001 = T2_fits_001 - repmat( T2_list(:), [1 6] );
T1_rel_bias_001 = T1_bias_001./repmat( T1_list(:), [1 6] );
T2_rel_bias_001 = T2_bias_001./repmat( T2_list(:), [1 6] );

save( [dir_out 'b0_mrf_rel_bias.mat'], 'T1_rel_bias_105','T1_rel_bias_104','T1_rel_bias_001',...
    'T2_rel_bias_105','T2_rel_bias_104','T2_rel_bias_001' );
%% get statistics of relative bias

n_set = 4;

nn = 1;
T1_mean_rel_bias_105 = zeros( size(T1_rel_bias_105,1)/n_set, 6 );
T2_mean_rel_bias_105 = T1_mean_rel_bias_105;
T1_mean_rel_bias_104 = T1_mean_rel_bias_105;
T2_mean_rel_bias_104 = T1_mean_rel_bias_105;
T1_mean_rel_bias_001 = T1_mean_rel_bias_105;
T2_mean_rel_bias_001 = T1_mean_rel_bias_105;
T1_std_rel_bias_105 = T1_mean_rel_bias_105;
T2_std_rel_bias_105 = T1_mean_rel_bias_105;
T1_std_rel_bias_104 = T1_mean_rel_bias_105;
T2_std_rel_bias_104 = T1_mean_rel_bias_105;
T1_std_rel_bias_001 = T1_mean_rel_bias_105;
T2_std_rel_bias_001 = T1_mean_rel_bias_105;
for ii = 1:n_set:size(T1_rel_bias_105,1)
    
    T1_mean_rel_bias_105(nn,:) = mean( T1_rel_bias_105(ii:ii+n_set-1,:), 1 );
    T2_mean_rel_bias_105(nn,:) = mean( T2_rel_bias_105(ii:ii+n_set-1,:), 1 );
    T1_std_rel_bias_105(nn,:) = std( T1_rel_bias_105(ii:ii+n_set-1,:),[], 1 );
    T2_std_rel_bias_105(nn,:) = std( T2_rel_bias_105(ii:ii+n_set-1,:),[], 1 );
    
    T1_mean_rel_bias_104(nn,:) = mean( T1_rel_bias_104(ii:ii+n_set-1,:), 1 );
    T2_mean_rel_bias_104(nn,:) = mean( T2_rel_bias_104(ii:ii+n_set-1,:), 1 );
    T1_std_rel_bias_104(nn,:) = std( T1_rel_bias_104(ii:ii+n_set-1,:),[], 1 );
    T2_std_rel_bias_104(nn,:) = std( T2_rel_bias_104(ii:ii+n_set-1,:),[], 1 );
    
    T1_mean_rel_bias_001(nn,:) = mean( T1_rel_bias_001(ii:ii+n_set-1,:), 1 );
    T2_mean_rel_bias_001(nn,:) = mean( T2_rel_bias_001(ii:ii+n_set-1,:), 1 );
    T1_std_rel_bias_001(nn,:) = std( T1_rel_bias_001(ii:ii+n_set-1,:),[], 1 );
    T2_std_rel_bias_001(nn,:) = std( T2_rel_bias_001(ii:ii+n_set-1,:),[], 1 );
    
    nn = nn + 1;
    
end


T1_TBW4C1_mean = [T1_mean_rel_bias_105(:,1) T1_mean_rel_bias_104(:,1) T1_mean_rel_bias_001(:,1)];
T2_TBW4C1_mean = [T2_mean_rel_bias_105(:,1) T2_mean_rel_bias_104(:,1) T2_mean_rel_bias_001(:,1)];
T1_TBW4C1_std = [T1_std_rel_bias_105(:,1) T1_std_rel_bias_104(:,1) T1_std_rel_bias_001(:,1)];
T2_TBW4C1_std = [T2_std_rel_bias_105(:,1) T2_std_rel_bias_104(:,1) T2_std_rel_bias_001(:,1)];
T1_TBW4C1_lo = T1_TBW4C1_std;
T1_TBW4C1_hi = T1_TBW4C1_std;
T2_TBW4C1_lo = T2_TBW4C1_std;
T2_TBW4C1_hi = T2_TBW4C1_std;

T1_TBW4C4_mean = [T1_mean_rel_bias_105(:,2) T1_mean_rel_bias_104(:,2) T1_mean_rel_bias_001(:,2)];
T2_TBW4C4_mean = [T2_mean_rel_bias_105(:,2) T2_mean_rel_bias_104(:,2) T2_mean_rel_bias_001(:,2)];
T1_TBW4C4_std = [T1_std_rel_bias_105(:,2) T1_std_rel_bias_104(:,2) T1_std_rel_bias_001(:,2)];
T2_TBW4C4_std = [T2_std_rel_bias_105(:,2) T2_std_rel_bias_104(:,2) T2_std_rel_bias_001(:,2)];
T1_TBW4C4_lo = T1_TBW4C4_std;
T1_TBW4C4_hi = T1_TBW4C4_std;
T2_TBW4C4_lo = T2_TBW4C4_std;
T2_TBW4C4_hi = T2_TBW4C4_std;

T1_TBW4C8_mean = [T1_mean_rel_bias_105(:,3) T1_mean_rel_bias_104(:,3) T1_mean_rel_bias_001(:,3)];
T2_TBW4C8_mean = [T2_mean_rel_bias_105(:,3) T2_mean_rel_bias_104(:,3) T2_mean_rel_bias_001(:,3)];
T1_TBW4C8_std = [T1_std_rel_bias_105(:,3) T1_std_rel_bias_104(:,3) T1_std_rel_bias_001(:,3)];
T2_TBW4C8_std = [T2_std_rel_bias_105(:,3) T2_std_rel_bias_104(:,3) T2_std_rel_bias_001(:,3)];
T1_TBW4C8_lo = T1_TBW4C8_std;
T1_TBW4C8_hi = T1_TBW4C8_std;
T2_TBW4C8_lo = T2_TBW4C8_std;
T2_TBW4C8_hi = T2_TBW4C8_std;

T1_TBW8C1_mean = [T1_mean_rel_bias_105(:,4) T1_mean_rel_bias_104(:,4) T1_mean_rel_bias_001(:,4)];
T2_TBW8C1_mean = [T2_mean_rel_bias_105(:,4) T2_mean_rel_bias_104(:,4) T2_mean_rel_bias_001(:,4)];
T1_TBW8C1_std = [T1_std_rel_bias_105(:,4) T1_std_rel_bias_104(:,4) T1_std_rel_bias_001(:,4)];
T2_TBW8C1_std = [T2_std_rel_bias_105(:,4) T2_std_rel_bias_104(:,4) T2_std_rel_bias_001(:,4)];
T1_TBW8C1_lo = T1_TBW8C1_std;
T1_TBW8C1_hi = T1_TBW8C1_std;
T2_TBW8C1_lo = T2_TBW8C1_std;
T2_TBW8C1_hi = T2_TBW8C1_std;

T1_TBW8C4_mean = [T1_mean_rel_bias_105(:,5) T1_mean_rel_bias_104(:,5) T1_mean_rel_bias_001(:,5)];
T2_TBW8C4_mean = [T2_mean_rel_bias_105(:,5) T2_mean_rel_bias_104(:,5) T2_mean_rel_bias_001(:,5)];
T1_TBW8C4_std = [T1_std_rel_bias_105(:,5) T1_std_rel_bias_104(:,5) T1_std_rel_bias_001(:,5)];
T2_TBW8C4_std = [T2_std_rel_bias_105(:,5) T2_std_rel_bias_104(:,5) T2_std_rel_bias_001(:,5)];
T1_TBW8C4_lo = T1_TBW8C4_std;
T1_TBW8C4_hi = T1_TBW8C4_std;
T2_TBW8C4_lo = T2_TBW8C4_std;
T2_TBW8C4_hi = T2_TBW8C4_std;

T1_TBW8C8_mean = [T1_mean_rel_bias_105(:,6) T1_mean_rel_bias_104(:,6) T1_mean_rel_bias_001(:,6)];
T2_TBW8C8_mean = [T2_mean_rel_bias_105(:,6) T2_mean_rel_bias_104(:,6) T2_mean_rel_bias_001(:,6)];
T1_TBW8C8_std = [T1_std_rel_bias_105(:,6) T1_std_rel_bias_104(:,6) T1_std_rel_bias_001(:,6)];
T2_TBW8C8_std = [T2_std_rel_bias_105(:,6) T2_std_rel_bias_104(:,6) T2_std_rel_bias_001(:,6)];
T1_TBW8C8_lo = T1_TBW8C8_std;
T1_TBW8C8_hi = T1_TBW8C8_std;
T2_TBW8C8_lo = T2_TBW8C8_std;
T2_TBW8C8_hi = T2_TBW8C8_std;


%% plot results

c = {'0','1/4','1/2','5/4','3/2'};
co = [0 0 0 ; 0 1 0; 0 0 1];

% TBW4C1
% T1
figure(100); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T1_TBW4C1_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T1_TBW4C1_mean, 1);
nbars = size(T1_TBW4C1_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T1_TBW4C1_mean(:,ii), T1_TBW4C1_lo(:,ii), T1_TBW4C1_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-0.5 1])
ylabel('Relative bias')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
my_leg = legend('variable TE/fixed TR','fixed TE/fixed TR','fixed TE/variable TR','Location','northwest');
my_leg.Position(2) = 0.75;

% TBW4C1
% T2
figure(101); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T2_TBW4C1_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T2_TBW4C1_mean, 1);
nbars = size(T2_TBW4C1_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T2_TBW4C1_mean(:,ii), T2_TBW4C1_lo(:,ii), T2_TBW4C1_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-1 1])
ylabel('Relative bias in T_2')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
legend('variable TE/fixed TR','fixed TE/fixed TR','fixed TE/variable TR','Location','northwest');
legend boxoff;
title('TBW4 C1')
export_fig('../figures/b0_T2_bias_tbw4c1.png',sprintf('-r%d',my_dpi));


% TBW4C4
% T1
figure(110); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T1_TBW4C4_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T1_TBW4C4_mean, 1);
nbars = size(T1_TBW4C4_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T1_TBW4C4_mean(:,ii), T1_TBW4C4_lo(:,ii), T1_TBW4C4_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-0.5 1])
ylabel('Relative bias')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
my_leg = legend('variable TE','fixed TE','variable TR','Location','northwest');
my_leg.Position(2) = 0.75;

% TBW4C4
% T2
figure(111); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T2_TBW4C4_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T2_TBW4C4_mean, 1);
nbars = size(T2_TBW4C4_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T2_TBW4C4_mean(:,ii), T2_TBW4C4_lo(:,ii), T2_TBW4C4_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-1 1])
ylabel('Relative bias in T_2')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
title('TBW4 C4')
export_fig('../figures/b0_T2_bias_tbw4c4.png',sprintf('-r%d',my_dpi));

% TBW4C8
% T1
figure(120); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T1_TBW4C8_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T1_TBW4C8_mean, 1);
nbars = size(T1_TBW4C8_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T1_TBW4C8_mean(:,ii), T1_TBW4C8_lo(:,ii), T1_TBW4C8_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-0.5 1])
ylabel('Relative bias')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
my_leg = legend('variable TE','fixed TE','variable TR','Location','northwest');
my_leg.Position(2) = 0.75;

% TBW4C8
% T2
figure(121); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T2_TBW4C8_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T2_TBW4C8_mean, 1);
nbars = size(T2_TBW4C8_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T2_TBW4C8_mean(:,ii), T2_TBW4C8_lo(:,ii), T2_TBW4C8_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-1 1])
ylabel('Relative bias in T_2')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
title('TBW4 C8')
export_fig('../figures/b0_T2_bias_tbw4c8.png',sprintf('-r%d',my_dpi));

% TBW8C1
% T1
figure(200); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T1_TBW8C1_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T1_TBW8C1_mean, 1);
nbars = size(T1_TBW8C1_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T1_TBW8C1_mean(:,ii), T1_TBW8C1_lo(:,ii), T1_TBW8C1_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-0.5 1])
ylabel('Relative bias')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
my_leg = legend('variable TE','fixed TE','variable TR','Location','northwest');
my_leg.Position(2) = 0.75;

% TBW8C1
% T2
figure(201); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T2_TBW8C1_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T2_TBW8C1_mean, 1);
nbars = size(T2_TBW8C1_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T2_TBW8C1_mean(:,ii), T2_TBW8C1_lo(:,ii), T2_TBW8C1_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-1 1])
ylabel('Relative bias in T_2')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
title('TBW8 C1')
export_fig('../figures/b0_T2_bias_tbw8c1.png',sprintf('-r%d',my_dpi));

% TBW8C4
% T1
figure(210); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T1_TBW8C4_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T1_TBW8C4_mean, 1);
nbars = size(T1_TBW8C4_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T1_TBW8C4_mean(:,ii), T1_TBW8C4_lo(:,ii), T1_TBW8C4_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-0.5 1])
ylabel('Relative bias')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
my_leg = legend('variable TE','fixed TE','variable TR','Location','northwest');
my_leg.Position(2) = 0.75;

% TBW8C4
% T2
figure(211); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T2_TBW8C4_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T2_TBW8C4_mean, 1);
nbars = size(T2_TBW8C4_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T2_TBW8C4_mean(:,ii), T2_TBW8C4_lo(:,ii), T2_TBW8C4_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-1 1])
ylabel('Relative bias in T_2')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
title('TBW8 C4')
export_fig('../figures/b0_T2_bias_tbw8c4.png',sprintf('-r%d',my_dpi));

% TBW8C8
% T1
figure(220); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T1_TBW8C8_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T1_TBW8C8_mean, 1);
nbars = size(T1_TBW8C8_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T1_TBW8C8_mean(:,ii), T1_TBW8C8_lo(:,ii), T1_TBW8C8_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-0.5 1])
ylabel('Relative bias')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
my_leg = legend('variable TE','fixed TE','variable TR','Location','northwest');
my_leg.Position(2) = 0.75;

% TBW8C8
% T2
figure(221); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( T2_TBW8C8_mean, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(T2_TBW8C8_mean, 1);
nbars = size(T2_TBW8C8_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, T2_TBW8C8_mean(:,ii), T2_TBW8C8_lo(:,ii), T2_TBW8C8_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([-1 1])
ylabel('Relative bias in T_2')
set(gca,'XTickLabel',c,'FontSize',18);
xlabel('\DeltaB_0 [(TR)^{-1}]');
title('TBW8 C8')
export_fig('../figures/b0_T2_bias_tbw8c8.png',sprintf('-r%d',my_dpi));

%% get max mag. of relative biases for fixed TE/fixed TR for each TBW/C combination

max_T2_bias_MRF104 = max( abs(T2_rel_bias_104) );
save( [dir_out 'b0_max_T2_bias_fixed_TE_fixed_TR.mat'], 'max_T2_bias_MRF104' );