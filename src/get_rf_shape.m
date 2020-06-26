function [rf_shape] = get_rf_shape( shape_in,fa_deg,plot_flag )
% create rf_shape structure for use of ssEPG
% INPUT shape_in struct with fields:
% fn: (string) rf shape input file
% N_t: number of points used to define RF envelope
% sl_thick_m: nominal slice thickness in m
% tau_s: RF pulse duration in s
% n_cycles_p_crush: (int) number of cycles introduced by slice select crusher
% N_res: (int) number of points in fov_factor * sl_thick_m (used for EPG
% state matrix determination)
% fov_factor: see N_res

% fa_deg: nominal flip-angle in degrees
% plot_flag: if 1 then plot RF

% OUTPUT
% rf_shape: struct contains RF amplitude modulation scaled in units T

fn = shape_in.fn;
N_t = shape_in.N_t;
sl_thick_m = shape_in.sl_thick_m;
tau_s = shape_in.tau_s;
n_cycles_p_crush = shape_in.n_cycles_per_crush;
N_res = shape_in.N_res;
fov_factor = shape_in.fov_factor;

gammabar_Hz_p_T = 42.57747852e6;

% load rf shape file
fid = fopen( fn, 'r' );
dat = fscanf( fid, '%f' );
fclose( fid );
fa_rad = fa_deg*pi/180;

dur_s = dat(1)/1000;
ttipdown_s = dat(5)/1000;
TBW = dat(6);
nr_samples = dat(7);
rf_am = dat(8:8+nr_samples);

refocus_factor = ttipdown_s/dur_s;

% interpolate and scale
N_rf_og = numel( rf_am );
rf_am = interp1( (-N_rf_og/2:N_rf_og/2 -1)./N_rf_og, rf_am, (-N_t/2:N_t/2 -1)./N_t, 'pchip', 0 );

rf_am = rf_am./sum(rf_am);

grad_T_p_m = ( TBW/tau_s ) / gammabar_Hz_p_T / sl_thick_m;

dtau_s = tau_s/N_t;

N_cycle = -( N_t * ( n_cycles_p_crush / TBW - refocus_factor ) );

if mod( N_cycle, 1 ) ~= 0
    disp( 'get_rf_shape: Warning! N_cycle is non-integer and is being rounded.' )
    N_cycle = round( N_cycle );
end

t_s = (-N_t/2:N_t/2 - 1)*dtau_s;

b1_T = fa_rad / ( 2 * pi * gammabar_Hz_p_T * dtau_s );
rf_T = b1_T * rf_am; 

diso_s = tau_s * refocus_factor - tau_s / 2; diso_idx = round( diso_s/dtau_s );

Q = round( N_res * N_t / fov_factor / TBW / 2 ) + 1; 

if plot_flag == 1
    
    figure(99); clf;
    plot( t_s*1e3, rf_T*1e6 )
    xlabel('time (ms)','fontsize',16)
    ylabel('B1(t) (\muT)','fontsize',16)
    xlim([min(t_s) max(t_s)]*1000)
    drawnow
    
end

rf_shape.N_t = N_t;
rf_shape.dtau_s = dtau_s;
rf_shape.B1_T = rf_T;
rf_shape.FA_nom_deg = fa_deg;
rf_shape.diso_idx = diso_idx;
rf_shape.Q = Q;
rf_shape.N_cycle = N_cycle;
rf_shape.sl_thick_m = sl_thick_m;
rf_shape.grad_T_p_m = grad_T_p_m;
rf_shape.shape_in = shape_in;
rf_shape.TBW = TBW;

end

