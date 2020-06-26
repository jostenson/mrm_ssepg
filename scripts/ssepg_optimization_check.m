%% Evalulation of ssEPG under reduced resolution
% (delta Q / delta N check)

%% initial parameters

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_pulse{1} = 'hanning_sinc_ex_tbw4.txt';
fn_pulse{2} = 'hanning_sinc_ex_tbw8.txt';
TBW = [4 8];
tau_s = [2.0 3.5]/1000;
sl_thick_m = 5/1000;

fn_MRF_seq_params = 'MRF105.csv';
n_cycles_per_crush = [1 2 4];
nomFlip_deg = 60;
TI_s = 40/1000; % ms
TRbase_s = 16/1000;
TE_s = 3/1000;
nreps = 1250;

N_t0 = [128 256];
fov_factor = 4;
Q0 = 3201;
N_res0 = ceil( fov_factor * 2 / 32 * (Q0-1) ); %  fov_fact * 2 * tbw / n_t * (Q - 1)

n_iter = 100; % rand T1 and T2 check of model accuracy
n_resos = 15; % number of resolutions (Q/deltaNs) to check
T1_max_s = 1500/1000;
T2_max_s = 150/1000;
T1_min_s = 300/1000;
T2_min_s = 5/1000;
delta_check = 0.02; % fractional difference in T1 or T2 allow for rand T1 and T2 on check of optimization

%% derivative params

n_pulses = numel( fn_pulse );
n_crush = numel( n_cycles_per_crush );
s = rng(31465);

%% plot MRF sequence params

data_struct = importdata([dir_in fn_MRF_seq_params]);
FA_v = data_struct.data(:,1)*nomFlip_deg; % deg
phi_v = data_struct.data(:,2); % deg
TR_v = TRbase_s + data_struct.data(:,3)/1000; %
TE_v = TE_s + data_struct.data(:,4)/1000; %

figure(1); clf;
subplot(411)
plot(FA_v); ylabel('degrees'); title(['FA, TI is ' num2str(TI_s*1000) ' ms'])
subplot(412)
plot(phi_v); ylabel('degrees'); title('phase')
subplot(413)
plot(TR_v*1000); ylabel('msec'); title('TR')
subplot(414)
plot(TE_v*1000); ylabel('msec'); title('TE')
drawnow

%% GPU kernels

k2 = parallel.gpu.CUDAKernel('ssepg_functions.ptx','ssepg_functions.cu','dephase_gradients_rf_stage2');
kerns.k2 = k2;

%% check over closely spaced T1s and T2s as function of Q/deltaN

tic;

myT1_s = logspace( log10(T1_min_s), log10(T1_max_s), round(sqrt(n_iter)) );
myT2_s = logspace( log10(T2_min_s), log10(T2_max_s), round(sqrt(n_iter)) );

T1_list_msr = zeros( round(sqrt(n_iter))^2,1 );
T2_list_msr = T1_list_msr;
nn = 1;
for ii = 1:round(sqrt(n_iter))
    for jj = 1:round(sqrt(n_iter))
        T1_list_msr(nn) = myT1_s(ii);
        T2_list_msr(nn) = myT2_s(jj);
        nn = nn + 1;
    end
end

% get T1s and T2s around those
T1_list = [];
T2_list = [];
for kk = 1:3
    T1_factor = 1 - delta_check + (kk-1)*delta_check;
    for mm = 1:3
        T2_factor = 1 - delta_check + (mm-1)*delta_check;
        
        T1_list = [T1_list; T1_list_msr*T1_factor];
        T2_list = [T2_list; T2_list_msr*T2_factor];
        
    end
end

% test optimization
mean_rel_T1_errors = zeros( n_pulses,n_crush,n_resos );
std_rel_T1_errors = mean_rel_T1_errors;
mean_rel_T2_errors = mean_rel_T1_errors;
std_rel_T2_errors = mean_rel_T1_errors;
max_rel_T1_errors = mean_rel_T1_errors;
max_rel_T2_errors = mean_rel_T1_errors;
mean_abs_rel_T1_errors = mean_rel_T1_errors;
mean_abs_rel_T2_errors = mean_rel_T1_errors;
Q_p_dNs = mean_rel_T1_errors;
N_res_optis = zeros( n_pulses,n_crush );
Q_p_dN_optis = N_res_optis;

for ii = 1:n_pulses
    for jj = 1:n_crush
        
        fprintf('Doing T1 & T2 check for %d queries for RF TBW %d and N_crush %d ...\n',...
            n_iter,TBW(ii),n_cycles_per_crush(jj) );
        
        % get measured signals
        shape_in.fn = [dir_in fn_pulse{ii}];
        shape_in.N_t = N_t0(ii);
        shape_in.sl_thick_m = sl_thick_m;
        shape_in.tau_s = tau_s(ii);
        shape_in.n_cycles_per_crush = n_cycles_per_crush(jj);
        shape_in.fov_factor = fov_factor;
        shape_in.N_res = N_res0;
        
        rf_shape = get_rf_shape( shape_in, nomFlip_deg, 1 );
        
        delk = rf_shape.N_cycle;
        szomega = rf_shape.Q;
        [~,sigs_init] = EPG_MRF_SSFP_prof_par( T1_list_msr,T2_list_msr,TE_v,TR_v,FA_v,delk,nreps,szomega,phi_v,TI_s,rf_shape,n_iter,kerns,0);
        sigs_init = sigs_init.* repmat( 1./vecnorm(sigs_init), [nreps 1] );
        
        for kk = 1:n_resos
            
            Q_p_dN_approx = 2*(kk-1) + 4;            
            
            fprintf('Q/dN ~ %d ...\n',Q_p_dN_approx);
            
            Nr = 2 * n_cycles_per_crush(jj) * Q_p_dN_approx;
            my_N_res = fov_factor * Nr;
            
            % get dictionary using optimized values
            shape_in.N_res = my_N_res;
            rf_shape = get_rf_shape( shape_in, nomFlip_deg, 1 );
            
            delk = rf_shape.N_cycle;
            szomega = rf_shape.Q;
            
            [~,Ft] = EPG_MRF_SSFP_prof_par( T1_list,T2_list,TE_v,TR_v,FA_v,delk,nreps,szomega,phi_v,TI_s,rf_shape,9*n_iter,kerns,0);
            D = Ft.* repmat( 1./vecnorm(Ft), [nreps 1] );
            
            [~,idx_max] = max( sigs_init'*D,[],2 );
            
            T1_error = T1_list(idx_max) - T1_list_msr;
            T2_error = T2_list(idx_max) - T2_list_msr;
            
            mean_rel_T1_error = mean( T1_error./T1_list_msr );
            mean_rel_T2_error = mean( T2_error./T2_list_msr );
            mean_abs_rel_T1_error = mean( abs( T1_error./T1_list_msr ) );
            mean_abs_rel_T2_error = mean( abs( T2_error./T2_list_msr ) );
            std_rel_T1_error = std( T1_error./T1_list_msr );
            std_rel_T2_error = std( T1_error./T1_list_msr );
            max_rel_T1_error = max( abs(T1_error./T1_list_msr) );
            max_rel_T2_error = max( abs(T2_error./T2_list_msr) );
            
            fprintf('Random T1 and T2 check complete for RF TBW %d and N_crush %d. \n', ...
                TBW(ii),n_cycles_per_crush(jj) );
            fprintf(' Mean rel. T1 error %.4f :: std rel. T1 error %.4f\n', mean_rel_T1_error, std_rel_T1_error );
            fprintf(' Mean rel. T2 error %.4f :: std rel. T2 error %.4f\n', mean_rel_T2_error, std_rel_T2_error );
            fprintf(' Mean abs rel. T1 error %.4f :: mean abs rel. T2 error %.4f\n', mean_abs_rel_T1_error, mean_abs_rel_T2_error );
            fprintf(' Max. rel. T1 error %.4f :: max rel. T2 error %.4f\n', max_rel_T1_error, max_rel_T2_error );
            
            mean_rel_T1_errors(ii,jj,kk) = mean_rel_T1_error;
            std_rel_T1_errors(ii,jj,kk) = std_rel_T1_error;
            mean_rel_T2_errors(ii,jj,kk) = mean_rel_T2_error;
            std_rel_T2_errors(ii,jj,kk) = std_rel_T2_error;
            max_rel_T1_errors(ii,jj,kk) = max_rel_T1_error;
            max_rel_T2_errors(ii,jj,kk) = max_rel_T2_error;
            mean_abs_rel_T1_errors(ii,jj,kk) = mean_abs_rel_T1_error;
            mean_abs_rel_T2_errors(ii,jj,kk) = mean_abs_rel_T2_error;
            Q_p_dNs(ii,jj,kk) = 1./n_cycles_per_crush(jj) * ( my_N_res/2/fov_factor + TBW(ii)/N_t0(ii) );
            
            if max_rel_T1_error == 0 && max_rel_T2_error == 0
                N_res_optis(ii,jj) = my_N_res;
                Q_p_dN_optis(ii,jj) = Q_p_dNs(ii,jj,kk);
                break;
            end
            
        end
    end
end
t2 = toc;
fprintf('Random check run time %.1f s\n',t2);

%% plot results

figure(10); clf;
figure(11); clf;
nn = 1;
for ii = 1:n_pulses
    switch ii
        case 1
            my_color = 'b';
        case 2
            my_color = 'r';
        case 3
            my_color = 'k';
    end
    for jj = 1:n_crush
        
        switch jj
            case 1
                my_mark = 'x';
            case 2
                my_mark = 'v';
        end
        
        my_mean_T1_error = squeeze( mean_abs_rel_T1_errors(ii,jj,:) );
        my_mean_T2_error = squeeze( mean_abs_rel_T2_errors(ii,jj,:) );
        my_QdN = squeeze( Q_p_dNs(ii,jj,:) );
   
        f10 = figure(10);
        plot( my_QdN( my_QdN > 0 ), my_mean_T1_error( my_QdN > 0 ), [my_color my_mark] );
        xlabel( '$\frac{Q}{\Delta N}$','Interpreter','latex');
        ylabel( 'Mean abs. rel. error' )
        hold on
        
        f11 = figure(11);
        plot( my_QdN( my_QdN > 0 ), my_mean_T2_error( my_QdN > 0 ), [my_color my_mark] );
        xlabel( '$\frac{Q}{\Delta N}$','Interpreter','latex');
        ylabel( 'Mean abs. rel. error' )
        hold on
       
                
        
    end
end

figure(11);
legend({'TBW4/C1','TBW4/C2','TBW4/C4','TBW8/C1','TBW8/C2','TBW8/C8'})



%%
fn_save = sprintf('mrf_ssepg_opti_check_%d.mat',nomFlip_deg);

save( [dir_out fn_save], 'nomFlip_deg','fn_pulse','TBW',...
    'tau_s','sl_thick_m','fn_MRF_seq_params','n_cycles_per_crush','nomFlip_deg',...
    'TI_s','TRbase_s','TE_s','nreps','N_t0','fov_factor','Q0',...
    'N_res0','n_iter','T1_max_s','T2_max_s','T1_min_s','T2_min_s','T1_list_msr','T2_list_msr',...
    'delta_check', 'mean_rel_T1_errors','mean_rel_T2_errors','std_rel_T1_errors','std_rel_T2_errors', ...
    'max_rel_T1_errors','max_rel_T2_errors','mean_abs_rel_T1_errors','mean_abs_rel_T2_errors', ...
    'Q_p_dNs','N_res_optis','Q_p_dN_optis' );