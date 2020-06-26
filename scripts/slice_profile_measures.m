%% ssEPG, pEPG, and EPG against measured slice profiles
% for serveral flip angles, TBWs, and crusher strenghts

my_dpi = 250;

%% parameters

N_p = 256; % number of partitions for pEPG method

dir_in = '../data_in/';

idx_profile_example = 6;

for nn = 1:2 % loop over crusher strengths
    
    if nn == 1
        fn{1} = 'SLPROF_TBW2C1_FA30.mat';
        fn{2} = 'SLPROF_TBW2C1_FA60.mat';
        fn{3} = 'SLPROF_TBW2C1_FA90.mat';
        fn{4} = 'SLPROF_TBW4C1_FA30.mat';
        fn{5} = 'SLPROF_TBW4C1_FA60.mat';
        fn{6} = 'SLPROF_TBW4C1_FA90.mat';
        fn{7} = 'SLPROF_TBW8C1_FA30.mat';
        fn{8} = 'SLPROF_TBW8C1_FA60.mat';
        fn{9} = 'SLPROF_TBW8C1_FA90.mat';
    else
        fn{1} = 'SLPROF_TBW2_FA30.mat';
        fn{2} = 'SLPROF_TBW2_FA60.mat';
        fn{3} = 'SLPROF_TBW2_FA90.mat';
        fn{4} = 'SLPROF_TBW4_FA30.mat';
        fn{5} = 'SLPROF_TBW4_FA60.mat';
        fn{6} = 'SLPROF_TBW4_FA90.mat';
        fn{7} = 'SLPROF_TBW8_FA30.mat';
        fn{8} = 'SLPROF_TBW8_FA60.mat';
        fn{9} = 'SLPROF_TBW8_FA90.mat';
    end
    
    % phantom properties
    T1_s = 419/1000;
    T2_s = 54/1000;
    fov_factor = 32/6;
    
    if nn == 1
        B1_corr = 0.960; %C1
    else
        B1_corr = 0.975; %C4
    end
    
    % sequence properties
    seq_params(1).pulse_type = 'hanning_sinc_ex_tbw2.txt';
    seq_params(1).tau_s = 0.40e-3;
    seq_params(1).sl_thick_m = 6.0e-3;
    seq_params(1).fa_deg = 30;
    seq_params(1).N_bw = 0;
    seq_params(1).n_reps = 10;
    
    seq_params(2).pulse_type = 'hanning_sinc_ex_tbw2.txt';
    seq_params(2).tau_s = 0.768e-3;
    seq_params(2).sl_thick_m = 6e-3;
    seq_params(2).fa_deg = 60;
    seq_params(2).N_bw = 0;
    seq_params(2).n_reps = 10;
    
    seq_params(3).pulse_type = 'hanning_sinc_ex_tbw2.txt';
    seq_params(3).tau_s = 1.152e-3;
    seq_params(3).sl_thick_m = 6e-3;
    seq_params(3).fa_deg = 90;
    seq_params(3).N_bw = 0;
    seq_params(3).n_reps = 10;
    
    seq_params(4).pulse_type = 'hanning_sinc_ex_tbw4.txt';
    seq_params(4).tau_s = 0.8e-3;
    seq_params(4).sl_thick_m = 6e-3;
    seq_params(4).fa_deg = 30;
    seq_params(4).N_bw = 0;
    seq_params(4).n_reps = 10;
    
    seq_params(5).pulse_type = 'hanning_sinc_ex_tbw4.txt';
    seq_params(5).tau_s = 1.246e-3;
    seq_params(5).sl_thick_m = 6e-3;
    seq_params(5).fa_deg = 60;
    seq_params(5).N_bw = 0;
    seq_params(5).n_reps = 10;
    
    seq_params(6).pulse_type = 'hanning_sinc_ex_tbw4.txt';
    seq_params(6).tau_s = 1.870e-3;
    seq_params(6).sl_thick_m = 6e-3;
    seq_params(6).fa_deg = 90;
    seq_params(6).N_bw = 0;
    seq_params(6).n_reps = 10;
    
    seq_params(7).pulse_type = 'hanning_sinc_ex_tbw8.txt';
    seq_params(7).tau_s = 1.60e-3;
    seq_params(7).sl_thick_m = 6e-3;
    seq_params(7).fa_deg = 30;
    seq_params(7).N_bw = 0;
    seq_params(7).n_reps = 10;
    
    seq_params(8).pulse_type = 'hanning_sinc_ex_tbw8.txt';
    seq_params(8).tau_s = 2.512e-3;
    seq_params(8).sl_thick_m = 6e-3;
    seq_params(8).fa_deg = 60;
    seq_params(8).N_bw = 0;
    seq_params(8).n_reps = 10;
    
    seq_params(9).pulse_type = 'hanning_sinc_ex_tbw8.txt';
    seq_params(9).tau_s = 3.766e-3;
    seq_params(9).sl_thick_m = 6e-3;
    seq_params(9).fa_deg = 90;
    seq_params(9).N_bw = 0;
    seq_params(9).n_reps = 10;
    
    if nn == 1
        seq_params(1).n_cycles_per_crush = 1;
        seq_params(2).n_cycles_per_crush = 1;
        seq_params(3).n_cycles_per_crush = 1;
        seq_params(4).n_cycles_per_crush = 1;
        seq_params(5).n_cycles_per_crush = 1;
        seq_params(6).n_cycles_per_crush = 1;
        seq_params(7).n_cycles_per_crush = 1;
        seq_params(8).n_cycles_per_crush = 1;
        seq_params(9).n_cycles_per_crush = 1;
    else
        seq_params(1).n_cycles_per_crush = 4;
        seq_params(2).n_cycles_per_crush = 4;
        seq_params(3).n_cycles_per_crush = 4;
        seq_params(4).n_cycles_per_crush = 4;
        seq_params(5).n_cycles_per_crush = 4;
        seq_params(6).n_cycles_per_crush = 4;
        seq_params(7).n_cycles_per_crush = 4;
        seq_params(8).n_cycles_per_crush = 4;
        seq_params(9).n_cycles_per_crush = 4;
    end
    
    co = [0 0 0; 0.6 0.6 0; 0 0.75 0.75]; % line/bar colors
    
    %% get measured set and examine across channels
    
    RMSE_table = zeros(3,numel(seq_params));
    for idx_select = 1:numel(seq_params)
        
        load([dir_in fn{idx_select}])
        
        data = output_MRF_raw.data;
        
        [n_dyn, n_ch, n_reps,n_rd] = size(data);
        
        data = squeeze( sum( data, 1 ) );
        
        profile_by_ch = zeros( n_ch, n_reps, n_rd-1 );
        for ii = 1:n_ch
            my_ch_data = squeeze( data(ii,:,:) );
            my_Ch_Data = ifftshift( fft( fftshift( my_ch_data(:,2:end), 2 ), [], 2 ), 2 );
            profile_by_ch(ii,:,:) = my_Ch_Data;
        end
        
        data_profs_combined = abs( squeeze( profile_by_ch(1,:,:) - 1i*profile_by_ch(2,:,:) ) ).';
        
        dz_m =fov_factor*seq_params(idx_select).sl_thick_m / (n_rd-1);
        z_measure_m = -fov_factor/2*seq_params(idx_select).sl_thick_m:dz_m:fov_factor/2*seq_params(idx_select).sl_thick_m - dz_m;
        
        figure(2); clf;
        for ii = 1:n_reps
            plot( z_measure_m*1000,data_profs_combined(:,ii));
            hold on;
            xlabel('slice position (mm')
            pause(0);
        end
        hold off
        
        %% plot measurement results
        
        sig_mag = abs( sum( squeeze( profile_by_ch(1,:,:) - 1i * profile_by_ch(2,:,:) ), 2 ) );
        sig_mag = sig_mag./norm(sig_mag);
        
        figure(5); clf;
        plot( sig_mag );
        xlabel('rep')
        
        
        %% EPG w/slice profile (ssEPG)
        
        nomFA_deg = output_MRF_raw.params.nomFlip_deg;
        
        shape_in.fn = [dir_in seq_params(idx_select).pulse_type];
        shape_in.N_t = 256;
        shape_in.sl_thick_m = seq_params(idx_select).sl_thick_m;
        shape_in.tau_s = seq_params(idx_select).sl_thick_m;
        shape_in.n_cycles_per_crush = seq_params(idx_select).n_cycles_per_crush;
        shape_in.fov_factor = fov_factor;
        shape_in.N_res = n_rd - 1;
        
        rf_shape = get_rf_shape( shape_in, nomFA_deg, 1 );
        
        FA_v = nomFA_deg*ones(n_reps,1);
        TR_v = output_MRF_raw.params.TRbase_s*ones(n_reps,1);
        TE_v = output_MRF_raw.params.TE_s*ones(n_reps,1);
        phi_v = zeros(n_reps,1);
        TI = inf;
        delk = rf_shape.N_cycle;
        szomega = rf_shape.Q;
        
        tic;
        [Fxy,sig_epg_slice_v] = EPG_MRF_SSFP_prof( T1_s, T2_s, TE_v, TR_v, B1_corr*FA_v, delk, n_reps, szomega, phi_v, TI, rf_shape );
        toc;
        
        sig_epg_slice_v = sig_epg_slice_v/norm(sig_epg_slice_v);
        Mxy = fftshift( ifft( fftshift( Fxy,2 ), [], 2 ), 2 ) * ( 2 * rf_shape.Q - 2 );
        idx_start = rf_shape.Q - shape_in.N_res/2;
        idx_end = rf_shape.Q - 1 + shape_in.N_res/2;
        Mxy_crop = Mxy(:,idx_start:idx_end);
        
        %% model time series with EPG
        
        delk = 1;
        szomega = 101;
        
        tic;
        sig_epg_v = EPG_MRF_SSFP( T1_s, T2_s, TE_v, TR_v, B1_corr*FA_v, delk, n_reps, szomega, phi_v, TI);
        toc;
        
        sig_epg_v = sig_epg_v/norm(sig_epg_v);
        
        %% model time series with pEPG
        
        % get effective FA for each partition
        first_profile = Mxy_crop(1,:);
        fa_slice_factors = asind( imag(first_profile)*exp(TE_v(1)/T2_s) )/nomFA_deg; % adjust for T2, B1 correction already included here
        
        % determine partiton spacing
        idx_partition = round(1:(numel(fa_slice_factors)/N_p):numel(fa_slice_factors));
        
        % compute pEPG signals
        sig_pepg_v = zeros(n_reps,1);
        sig_pepg_full = zeros(n_reps,N_p);
        tic;
        for ii = 1:N_p
            my_factor = fa_slice_factors( idx_partition(ii) );
            my_pepg_sig_v = EPG_MRF_SSFP( T1_s, T2_s, TE_v, TR_v, my_factor*FA_v, delk, n_reps, szomega, phi_v, TI);
            sig_pepg_v = my_pepg_sig_v(:) + sig_pepg_v;
            sig_pepg_full(:,ii) = my_pepg_sig_v;
        end
        toc;
        
        sig_pepg_v = sig_pepg_v/norm(sig_pepg_v);
        
        
        %% plot results
        
        pulse_type = seq_params(idx_select).pulse_type(1:end-4);
        sl_thick_m = seq_params(idx_select).sl_thick_m;
        
        figure(100); clf;
        plot( abs( sig_mag ) );
        hold on;
        plot( abs( sig_epg_slice_v ),'--','Color',co(1,:));
        plot( abs( sig_pepg_v ),'-.','Color',co(2,:));
        plot( abs( sig_epg_v ),':','Color',co(3,:));
        hold off;
        xlabel('repetition');
        ylabel('signal (au)');
        title( sprintf( 'Pulse %s :: Nominal FA %.1f deg', pulse_type, nomFA_deg),'Interpreter','none');
        legend('measured','ssEPG','pEPG','EPG')
        
        resids_ssEPG = abs(sig_mag(:)) - abs(sig_epg_slice_v(:));
        RMSE_ssEPG = sqrt( sum( resids_ssEPG(:)'*resids_ssEPG(:) / numel(resids_ssEPG) ) );
        fprintf( 'ssEPG RMSE for %s with FA %.1f is %.3f\n', pulse_type, nomFA_deg, RMSE_ssEPG );
        
        resids_pEPG = abs(sig_mag(:)) - abs(sig_pepg_v(:));
        RMSE_pEPG = sqrt( sum( resids_pEPG(:)'*resids_pEPG(:) / numel(resids_pEPG) ) );
        fprintf( 'pEPG RMSE for %s with FA %.1f is %.3f\n', pulse_type, nomFA_deg, RMSE_pEPG );
        
        resids_nosl = abs(sig_mag(:)) - abs(sig_epg_v(:));
        RMSE_nosl = sqrt( sum( resids_nosl(:)'*resids_nosl(:) / numel(resids_nosl) ) );
        fprintf( 'EPG RMSE for %s with FA %.1f is %.3f\n', pulse_type, nomFA_deg, RMSE_nosl );
        
        RMSE_table(:,idx_select) = [RMSE_ssEPG, RMSE_pEPG, RMSE_nosl]';
        
        figure(101); clf;
        set( gcf, 'Color', 'w', 'WindowState', 'maximize' );
        for ii = 1:n_reps
            
            if ii == 1
                norm_model = norm( Mxy_crop(ii,:) );
                norm_pepg = norm( sig_pepg_full(ii,:) );
                norm_measure = norm( data_profs_combined(:,ii) );
                max_sig = max( [ abs(Mxy_crop(:)/norm_model); abs(data_profs_combined(:)/norm_measure); abs(sig_pepg_full(:)/norm_pepg) ] );
            end
            
            figure(101);
            set(gcf,'Color','w');
            subplot(2,n_reps/2,ii)
            plot(z_measure_m*1000,abs(data_profs_combined(:,ii)/norm_measure))
            hold on
            plot(z_measure_m*1000,abs(Mxy_crop(ii,:)/norm_model),'--','Color',co(1,:))
            plot(z_measure_m*1000,abs(sig_pepg_full(ii,:)/norm_pepg),':','Color',co(2,:))
            hold off
            if ii == 1 || ii == 6
                ylabel('|M_{x,y}|');
            end
            xlabel('slice position (mm)');
            xlim([-sl_thick_m*fov_factor/5 sl_thick_m*fov_factor/5]*1000)
            ylim([0 1.1*max_sig])
            title( sprintf('Profile Rep. %d',ii) );
            if ii == 5 && nn == 2
                my_leg = legend('Measured','ssEPG','pEPG','Location','northeast');
                legend boxoff;
            end
            
            if idx_select == idx_profile_example && ( ii == 3 || ii == 8 )
                figure(109); clf;
                set(gcf,'Color','w');
                plot(z_measure_m*1000,abs(data_profs_combined(:,ii)/norm_measure))
                hold on
                plot(z_measure_m*1000,abs(Mxy_crop(ii,:)/norm_model),'--','Color',co(1,:))
                plot(z_measure_m*1000,abs(sig_pepg_full(ii,:)/norm_pepg),':','Color',co(2,:))
                hold off
                ylabel('|M_{x,y}|');
                xlabel('slice position (mm)');
                xlim([-sl_thick_m*fov_factor/5 sl_thick_m*fov_factor/5]*1000)
                ylim([0 1.1*max_sig])
                title( sprintf('Profile Rep. %d',ii) );
                if ii == 8
                    my_leg = legend('Measured','ssEPG','pEPG','Location','northeast');
                    legend boxoff;
                end
                txt = sprintf('../figures/measured_profile_%.2d.png',ii);
                export_fig(txt,sprintf('-r%d',my_dpi));
            end
            
        end
        
        pause(0)
        
        if idx_select == idx_profile_example
            set(gcf,'WindowState','minimize');
            figure(101);
            export_fig('../figures/slice_prof_measurement.png',sprintf('-r%d',my_dpi));
        end
        
    end
    
    %% plot RMSE results

    RMSE_ssepg = reshape( RMSE_table(1,:), [3 3] );
    RMSE_pepg = reshape( RMSE_table(2,:), [3 3] );
    RMSE_epg = reshape( RMSE_table(3,:), [3 3] );
    
    font_sz = 28;
    line_w = 2;
    my_xlabels = repmat( [2,4,8],[1 3] );
    my_ylabels = [30, 60, 90];

    figure(200); clf;
    set(gcf,'Position',[200 100 1400 800],'Color','w');
    imagesc( [RMSE_ssepg inf*ones(3,1) RMSE_pepg inf*ones(3,1) RMSE_epg] ); colormap( hot ); axis image; caxis( [0 0.15] );
    set( gca, 'LineWidth', line_w );
    set(gca,'TickLength',[0 0])
    set( gca,'XTick', [1:3,5:7,9:11], 'XTickLabel', my_xlabels, 'FontSize', font_sz);
    set( gca,'YTick', [1:3], 'YTickLabel', my_ylabels, 'FontSize', font_sz);
    xlabel('TBW','FontSize',font_sz)
    ylabel(['FA (' char(176) ')'],'FontSize',font_sz);
    text(1.5,0.25,'ssEPG','FontSize',font_sz)
    text(5.5,0.25,'pEPG','FontSize',font_sz)
    text(9.65,0.25,'EPG','FontSize',font_sz)
    if nn == 1
        c = colorbar('south');
        c.Position(1) = 0.25;
        c.Position(2) = 0.10;
        c.Position(3) = 0.5;
        c.TickLabels(1) = {'0.0'};
        c.Label.String = 'RMSE';
        c.Label.Position(1) = 0.08;
        c.Label.Position(2) = -2;
        c.LineWidth = line_w;
        export_fig('../figures/Figure3a.png',sprintf('-r%d',my_dpi));
    else
        export_fig('../figures/Figure3b.png',sprintf('-r%d',my_dpi));
    end
    
end