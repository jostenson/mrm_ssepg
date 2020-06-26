%% validate ssEPG against ODE Bloch simulation

my_dpi = 250;

%% parameters

dir_in = '../data_in/';
dir_out = '../data_out/';
fn_pulse = 'hanning_sinc_ex_tbw4.txt';

% RF params

nomFA_deg = 90;
TBW = 4;

shape_in.fn = [dir_in fn_pulse];
shape_in.N_t = 256;
shape_in.sl_thick_m = 8/1000;
shape_in.tau_s = 3.5/1000;
shape_in.n_cycles_per_crush = 4;
shape_in.fov_factor = 4;
shape_in.N_res = 5000;

% sequeence params
T1_s = 1000/1000;
T2_s = 100/1000;
TR_s = 15/1000;
TE_s = 3/1000;
n_reps = 10;

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

%% ssEPG profile

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


%% repeat flip angle ODE Bloch

fprintf('Doing repeated excitations with crushers using ODE Bloch...\n')

% model
M0 = repmat( [0 0 1]', [1 N_res] );
sig = zeros(n_reps,1);
sig_tot = zeros(n_reps,N_res);
crusher_grad = exp( 1i * 2 * pi * gammabar_Hz_p_T * ( -A_crush_T_s_p_m + grad_T_p_m * tau_s * refocus_factor ) * z_m  ).';
tic;
for nn = 1:n_reps
    
    M_solved = zeros( N_res,3 );
    
    for ii = 1:N_res
        [~,M] = ode45( @(t,m) bloch_eq(t,m,rf_shape.B1_T,sl_Hz(ii),t_s),[min(t_s) max(t_s)],M0(:,ii));
        M_solved(ii,:) = M(end,:);
    end
    
    Mxy = M_solved(:,1) + 1i*M_solved(:,2);
    Mxy = Mxy.* exp( 1i * 2 * pi * gammabar_Hz_p_T * grad_T_p_m * z_m(:) * (tau_s * refocus_factor) );
    
    Mxy = Mxy.* exp(-TE_s/T2_s);

    sig(nn) = sum(Mxy);
    
    figure(10); clf;
    subplot(121)
    plot( z_m*1e3, abs(Mxy) )
    xlabel('slice position (mm)','FontSize',16); ylabel('M_{x,y}','FontSize',16);
    title( sprintf('Trans mag rep %d',nn) );

    subplot(122)
    plot( z_m*1e3, M_solved(:,3) )
    xlabel('slice position (mm)','FontSize',16); ylabel('M_{z}','FontSize',16);
    title( sprintf('Long mag rep %d',nn) );
    
    drawnow();
    
    sig_tot( nn,: ) = Mxy;
    Mxy = Mxy.* crusher_grad;
    M0(1:2,:) = [real(Mxy(:).'); imag(Mxy(:).')] * exp(-(TR_s-TE_s)/T2_s);
    M0(3,:) = M_solved(:,3).';
    M0(3,:) = ones(1,N_res).*( 1 - exp( -TR_s/T1_s ) ) + M0(3,:)*exp( -TR_s/T1_s );
    
end
fprintf('Time for Bloch time series %.1f ms\n', 1000*toc);


%% model time series with EPG

FA_v = nomFA_deg*ones(n_reps,1);
TR_v = TR_s*ones(n_reps,1);
TE_v = TE_s*ones(n_reps,1);
phi_v = zeros(n_reps,1);
TI = inf;
delk = 1;
szomega = 101;

tic;
sig_epg = EPG_MRF_SSFP( T1_s, T2_s, TE_v, TR_v, FA_v, delk, n_reps, szomega, phi_v, TI);
fprintf('Time for EPG time series %.1f ms\n', 1000*toc);

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

%% plot results

co = [255/255 140/255 0; 0 0 0];

figure(15); clf;
set(gcf,'Color','w','WindowState','maximized')
for ii = 1:n_reps
    
    figure(15);
    Mxy = 2 * (Q-1) * sig_epg_prof_tot(ii,:);
    
    if ii == 1
            max_sig = max( [ abs( sig_tot(:) ); abs( 2 * (Q-1) * sig_epg_prof_tot(:) ) ] );
    end
    
    subplot(2,n_reps/2,ii)
    plot( z_m*1e3, abs( sig_tot(ii,:) ), 'Color', co(1,:) )
    hold on
    plot( z_m*1e3,abs( Mxy(idx_bw_match) ), '--', 'Color', co(2,:) )   
    xlabel('slice position (mm)');
    if ii == 1 || ii == 6
        ylabel('|M_{x,y}|');
    end
    title( sprintf('Profile Rep. %d',ii) );
    xlim([-fov_factor/4*sl_thick_m fov_factor/4*sl_thick_m]*1000)
%     ylim([0 max_sig*1.2])
    ylim([0 1])
    hold off
    
    drawnow()
    
    if ii == 3 || ii == 8
       figure(19); clf;
       set( gcf,'Color','w');
       plot( z_m*1e3, abs( sig_tot(ii,:) ), 'Color', co(1,:) )
       hold on
       plot( z_m*1e3,abs( Mxy(idx_bw_match) ), '--', 'Color', co(2,:) )
       xlabel('slice position (mm)');
       ylabel('|M_{x,y}|');
       title( sprintf('Profile Rep. %d',ii) );
       xlim([-fov_factor/4*sl_thick_m fov_factor/4*sl_thick_m]*1000)
       ylim([0 1])
       hold off
       if ii == 8
           my_leg = legend('Bloch','ssEPG','Location','northeast');
           legend boxoff;
       end
       txt = sprintf( '../figures/numerical_sim_detail_%.2d.png',ii);
       export_fig( txt, sprintf('-r%d',my_dpi) );
    end
    
    sig_epg_prof_tot_crop(ii,:) = Mxy( idx_bw_match );
        
end

figure(15);
subplot(2,n_reps/2,5);
my_leg = legend('Bloch','ssEPG','Location','northeast');
legend boxoff;
export_fig('../figures/Figure1.png',sprintf('-r%d',my_dpi));

%% stats and more plots

resid_epg = (-sig_epg(:)*N_res/fov_factor) - sig(:);
resid_epg_prof = sig_epg_prof(:) - sig(:);

sig_norm = sig(:)./norm(sig);
sig_epg_norm = -sig_epg(:)./norm(sig_epg);
sig_epg_prof_norm = sig_epg_prof(:)./norm(sig_epg_prof);

resid_epg_norm = sig_epg_norm - sig_norm;
resid_epg_prof_norm = sig_epg_prof_norm - sig_norm;

RMSE_epg = sqrt( resid_epg'*resid_epg/ n_reps );
RMSE_epg_prof = sqrt( resid_epg_prof'*resid_epg_prof/ n_reps );

RMSE_epg_norm = sqrt( resid_epg_norm'*resid_epg_norm/ n_reps );
RMSE_epg_prof_norm = sqrt( resid_epg_prof_norm'*resid_epg_prof_norm/ n_reps );

fprintf( 'RMSE of EPG w/o profile %.1f\n', RMSE_epg );
fprintf( 'RMSE of EPG w/ profile %.1f\n', RMSE_epg_prof );

fprintf( 'RMSE of EPG normed w/o profile %.3f\n', RMSE_epg_norm );
fprintf( 'RMSE of EPG normed w/ profile %.3f\n', RMSE_epg_prof_norm );

figure(20); clf;
set(gcf,'Position',[200 0 800 600],'Color','w')
plot(abs(sig))
hold on
plot(abs(sig_epg*N_res/2),'--')
plot(abs(sig_epg_prof),'-.')
hold off
title(sprintf('TBW %d : FA %.1f deg : T1 %0.f ms : T2 %0.f ms', TBW, nomFA_deg, T1_s*1000, T2_s*1000));
xlabel('Rep.'); ylabel('signal (au)');
legend('Bloch','EPG','ssEPG')

figure(21); clf;
plot(real(sig))
hold on
plot(real(sig_epg*N_res/2),'--')
plot(real(sig_epg_prof),'-.')
hold off
title('Real signal')
xlabel('rep'); ylabel('signal');
legend('Bloch','EPG','ssEPG')

figure(22); clf;
plot(imag(sig))
hold on
plot(-imag(sig_epg*N_res/2),'--')
plot(imag(sig_epg_prof),'-.')
hold off
title('Imag signal')
xlabel('rep'); ylabel('signal');
legend('Bloch','EPG','ssEPG')

figure(23); clf;
set(gcf,'Color','w','WindowState','maximize')
plot(abs(sig_norm))
hold on
plot(abs(sig_epg_norm),'--')
plot(abs(sig_epg_prof_norm),'-.')
hold off
grid on;
xlabel('Repetition'); ylabel('Normed signal (au)');
legend('Bloch','EPG','ssEPG')
t1 = text(5,0.68,sprintf( 'RMSE of ssEPG = %.3f',RMSE_epg_prof_norm ),'FontSize',28);
t2 = text(5,0.62,sprintf( 'RMSE of EPG =    %.3f',RMSE_epg_norm ),'FontSize',28);
export_fig('../figures/Signals_norm_Bloch_v_ssEPG.png');


