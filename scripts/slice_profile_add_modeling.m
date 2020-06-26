%% Additional profile modeling in the steady state

my_dpi = 250;
font_sz = 28;
line_w = 2;
    
%% parameters

N_p = 256; % number of partitions for pEPG method
n_rd = 257; % match measurement
TRbase_s = 24/1000;
TE_s = 1.1/1000;

dir_in = '../data_in/';

for nn = 1:2 % loop over crusher strengths
    
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
    seq_params(1).n_reps = 100;
    
    seq_params(2).pulse_type = 'hanning_sinc_ex_tbw2.txt';
    seq_params(2).tau_s = 0.768e-3;
    seq_params(2).sl_thick_m = 6e-3;
    seq_params(2).fa_deg = 60;
    seq_params(2).N_bw = 0;
    seq_params(2).n_reps = 100;
    
    seq_params(3).pulse_type = 'hanning_sinc_ex_tbw2.txt';
    seq_params(3).tau_s = 1.152e-3;
    seq_params(3).sl_thick_m = 6e-3;
    seq_params(3).fa_deg = 90;
    seq_params(3).N_bw = 0;
    seq_params(3).n_reps = 100;
    
    seq_params(4).pulse_type = 'hanning_sinc_ex_tbw4.txt';
    seq_params(4).tau_s = 0.8e-3;
    seq_params(4).sl_thick_m = 6e-3;
    seq_params(4).fa_deg = 30;
    seq_params(4).N_bw = 0;
    seq_params(4).n_reps = 100;
    
    seq_params(5).pulse_type = 'hanning_sinc_ex_tbw4.txt';
    seq_params(5).tau_s = 1.246e-3;
    seq_params(5).sl_thick_m = 6e-3;
    seq_params(5).fa_deg = 60;
    seq_params(5).N_bw = 0;
    seq_params(5).n_reps = 100;
    
    seq_params(6).pulse_type = 'hanning_sinc_ex_tbw4.txt';
    seq_params(6).tau_s = 1.870e-3;
    seq_params(6).sl_thick_m = 6e-3;
    seq_params(6).fa_deg = 90;
    seq_params(6).N_bw = 0;
    seq_params(6).n_reps = 100;
    
    seq_params(7).pulse_type = 'hanning_sinc_ex_tbw8.txt';
    seq_params(7).tau_s = 1.60e-3;
    seq_params(7).sl_thick_m = 6e-3;
    seq_params(7).fa_deg = 30;
    seq_params(7).N_bw = 0;
    seq_params(7).n_reps = 100;
    
    seq_params(8).pulse_type = 'hanning_sinc_ex_tbw8.txt';
    seq_params(8).tau_s = 2.512e-3;
    seq_params(8).sl_thick_m = 6e-3;
    seq_params(8).fa_deg = 60;
    seq_params(8).N_bw = 0;
    seq_params(8).n_reps = 100;
    
    seq_params(9).pulse_type = 'hanning_sinc_ex_tbw8.txt';
    seq_params(9).tau_s = 3.766e-3;
    seq_params(9).sl_thick_m = 6e-3;
    seq_params(9).fa_deg = 90;
    seq_params(9).N_bw = 0;
    seq_params(9).n_reps = 100;
    
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
    
    RMSD_table = zeros(2,numel(seq_params));
    ss_table = zeros(2,numel(seq_params));
    for idx_select = 1:numel(seq_params)
        
        
        dz_m =fov_factor*seq_params(idx_select).sl_thick_m / (n_rd-1);
        z_measure_m = -fov_factor/2*seq_params(idx_select).sl_thick_m:dz_m:fov_factor/2*seq_params(idx_select).sl_thick_m - dz_m;
        
        
        %% EPG w/slice profile (ssEPG)
        
        nomFA_deg = seq_params(idx_select).fa_deg;
        n_reps = seq_params(idx_select).n_reps;
        
        shape_in.fn = [dir_in seq_params(idx_select).pulse_type];
        shape_in.N_t = 256;
        shape_in.sl_thick_m = seq_params(idx_select).sl_thick_m;
        shape_in.tau_s = seq_params(idx_select).sl_thick_m;
        shape_in.n_cycles_per_crush = seq_params(idx_select).n_cycles_per_crush;
        shape_in.fov_factor = fov_factor;
        shape_in.N_res = n_rd - 1;
        
        rf_shape = get_rf_shape( shape_in, nomFA_deg, 1 );
        
        FA_v = nomFA_deg*ones(n_reps,1);
        TR_v = TRbase_s*ones(n_reps,1);
        TE_v = TE_s*ones(n_reps,1);
        phi_v = zeros(n_reps,1);
        TI = inf;
        delk = rf_shape.N_cycle;
        szomega = rf_shape.Q;
        
        tic;
        [Fxy,sig_epg_slice_v] = EPG_MRF_SSFP_prof( T1_s, T2_s, TE_v, TR_v, B1_corr*FA_v, delk, n_reps, szomega, phi_v, TI, rf_shape );
        toc;
        
        sig_epg_slice_v = sig_epg_slice_v * size(Fxy,2);
        Mxy = fftshift( ifft( fftshift( Fxy,2 ), [], 2 ), 2 ) * ( 2 * rf_shape.Q - 2 );
        idx_start = rf_shape.Q - shape_in.N_res/2;
        idx_end = rf_shape.Q - 1 + shape_in.N_res/2;
        Mxy_crop = Mxy(:,idx_start:idx_end);
        
        %% model time series with EPG
        
        tic;
        sig_epg_v = EPG_MRF_SSFP( T1_s, T2_s, TE_v, TR_v, B1_corr*FA_v, delk, n_reps, szomega, phi_v, TI);
        toc;
        
        sig_epg_v = sig_epg_v * N_p / fov_factor; % scale to pEPG
        
        %% model time series with pEPG
        
        delk = 1;
        szomega = 101;
        
        % get effective FA for each partition
        first_profile = Mxy_crop(1,:);
        fa_slice_factors = asind( imag(first_profile)*exp(TE_v(1)/T2_s) )/FA_v(1);
        
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
        
        
        %% plot results
        
        pulse_type = seq_params(idx_select).pulse_type(1:end-4);
        sl_thick_m = seq_params(idx_select).sl_thick_m;
        
        figure(100); clf;
        plot( abs( sig_epg_slice_v ),'--','Color',co(1,:));
        hold on;
        plot( abs( sig_pepg_v ),'-.','Color',co(2,:));
        hold off;
        xlabel('repetition');
        ylabel('signal (au)');
        title( sprintf( 'Pulse %s :: Nominal FA %.1f deg', pulse_type, nomFA_deg),'Interpreter','none');
        legend('ssEPG','pEPG')
        
        ref_sig = abs( sig_epg_slice_v(:) );
        
        resids_ssEPG = abs(ref_sig(:)) - abs(sig_epg_slice_v(:));
        RMSD_ssEPG = sqrt( sum( resids_ssEPG(:)'*resids_ssEPG(:) / numel(resids_ssEPG) ) );
        fprintf( 'ssEPG RMSD for %s with FA %.1f is %.3f\n', pulse_type, nomFA_deg, RMSD_ssEPG );
        
        resids_pEPG = abs(ref_sig(:)) - abs(sig_pepg_v(:));
        RMSD_pEPG = sqrt( sum( resids_pEPG(:)'*resids_pEPG(:) / numel(resids_pEPG) ) );
        fprintf( 'pEPG RMSD for %s with FA %.1f is %.3f\n', pulse_type, nomFA_deg, RMSD_pEPG );
        
        RMSD_table(:,idx_select) = [RMSD_ssEPG, RMSD_pEPG ]';
        ss_table(:,idx_select) = [abs(sig_epg_slice_v(end)), abs(sig_pepg_v(end))]';
        
        if idx_select == 1
            max_prof_sig = max( abs(Mxy_crop(end,:)),[],2 );
        end
        
        idx_order = [1,4,7,2,5,8,3,6,9];
        figure(101);
        set( gcf, 'Color', 'w', 'WindowState', 'maximize' );
        subplot(3,3,idx_order(idx_select))
        plot(z_measure_m*1000,abs(Mxy_crop(end,:)),'--','Color',co(1,:))
        hold on
        plot(z_measure_m*1000,abs(sig_pepg_full(end,:)),':','Color',co(2,:))
        hold off
        if idx_select == 1 || idx_select == 2 || idx_select == 3
            ylabel('|M_{x,y}|');
        end
        ylim([0 max_prof_sig*1.5])
        xlabel('slice position (mm)');
        xlim([-sl_thick_m*fov_factor/5 sl_thick_m*fov_factor/5]*1000)
        
        switch idx_select
            case 1
                text(-11,0.065,['FA 30' char(176)],'FontSize',font_sz,'Rotation',90);
                text(-1.5,0.3,'TBW 2','FontSize',font_sz)
                
            case 2
                text(-11,0.065,['FA 60' char(176)],'FontSize',font_sz,'Rotation',90);
                
            case 3
                text(-11,0.065,['FA 90' char(176)],'FontSize',font_sz,'Rotation',90);
            case 4
                text(-1.5,0.3,'TBW 4','FontSize',font_sz)
            case 7
                text(-1.5,0.3,'TBW 8','FontSize',font_sz)
                my_leg = legend('ssEPG','pEPG','Location','northeast');
                legend boxoff;
                my_leg.Position(1) = 0.8411;
                my_leg.Position(2) = 0.8597;
        end
        
        
    end
    if nn == 1
    	export_fig('../figures/ss_profs_c1.png',sprintf('-r%d',my_dpi));
    else
        export_fig('../figures/ss_profs_c4.png',sprintf('-r%d',my_dpi));
    end
    
    
    %% plot difference in steady-state signal
    
    my_xlabels = repmat( [2,4,8],[1 3] );
    my_ylabels = [30, 60, 90];
    
    ss_ssepg = reshape( ss_table(1,:), [3 3] );
    ss_pepg = reshape( ss_table(2,:), [3 3] );
    ss_ratio =  ss_pepg./ss_ssepg;
    
    figure(200); clf;
    set(gcf,'Position',[200 100 1400 800],'Color','w');
    imagesc( ss_ratio );  axis image;
    set( gca, 'LineWidth', line_w );
    set(gca,'TickLength',[0 0])
    set( gca,'XTick', [1:3,5:7,9:11], 'XTickLabel', my_xlabels, 'FontSize', font_sz);
    set( gca,'YTick', [1:3], 'YTickLabel', my_ylabels, 'FontSize', font_sz);
    xlabel('TBW','FontSize',font_sz)
    ylabel(['FA (' char(176) ')'],'FontSize',font_sz);
    c = colorbar();
    c.Label.String = 'Signal mag. ratio';
    c.LineWidth = line_w;
    if nn == 1
        save('../data_out/rel_ss_signal_pepg_to_ssepg_C1.mat','ss_ratio');
        export_fig('../figures/ss_sig_mag_c1.png',sprintf('-r%d',my_dpi));
    else
        save('../data_out/rel_ss_signal_pepg_to_ssepg_C4.mat','ss_ratio');
        export_fig('../figures/ss_sig_mac_c4.png',sprintf('-r%d',my_dpi));
    end
        
end

%% merge profile figures

im1 = imread('../figures/ss_profs_c1.png');
im4 = imread('../figures/ss_profs_c4.png');

im1 = imresize( im1,[size(im4,1) size(im4,2)] );
buff_im = 255*ones( 100, size(im1,2),3 );
im_tot = cat( 1,im1,buff_im,im4 );

figure(300); clf;
imshow(im_tot,'Border','tight');
text(80,80,'a','FontSize',28)
text(80,2625,'b','FontSize',28)
export_fig('../figures/ss_profs_c1and4.png',sprintf('-r%d',2*my_dpi));
