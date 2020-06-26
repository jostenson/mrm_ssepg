%% analyze MR system phantom estimates

my_dpi = 250;


%%
dir_in = '../data_out/';
dir_out = dir_in;
dir_conv = '../data_in/';

fn_epg{1} = 'MRF105_TBW4_FA60_crush1_proc_epg.mat';
fn_pepg{1} = 'MRF105_TBW4_FA60_crush1_proc_pepg.mat';
fn_ssepg{1} = 'MRF105_TBW4_FA60_crush1_proc_ssepg.mat';
fn_epg{2} = 'MRF105_TBW4_FA60_crush4_proc_epg.mat';
fn_pepg{2} = 'MRF105_TBW4_FA60_crush4_proc_pepg.mat';
fn_ssepg{2} = 'MRF105_TBW4_FA60_crush4_proc_ssepg.mat';
fn_epg{3} = 'MRF105_TBW8_FA60_crush1_proc_epg.mat';
fn_pepg{3} = 'MRF105_TBW8_FA60_crush1_proc_pepg.mat';
fn_ssepg{3} = 'MRF105_TBW8_FA60_crush1_proc_ssepg.mat';
fn_epg{4} = 'MRF105_TBW8_FA60_crush4_proc_epg.mat';
fn_pepg{4} = 'MRF105_TBW8_FA60_crush4_proc_pepg.mat';
fn_ssepg{4} = 'MRF105_TBW8_FA60_crush4_proc_ssepg.mat';

legend_txt = {'TBW 4 C 1','TBW 4 C 4','TBW 8 C 1','TBW 8 C 4'};

idx_T2_phys = 7:14; % indices measurements/specs at physiological range of T2

make_conv = 1;
fn_conv = 'Ostenson_MRsys_conventional_images.mat';
fn_conv_proc = 'Ostenson_MRsys_conv_proc.mat';
IR_times_ms_v = [50,100,200,300,500,750,1000,1500,2000,3000,6000];
TE_times_ms_v = 12*(1:16);
TD_ms = 2500;
n_T2_mod_params = 2;

temp_K = 273.15 + 21.5;
temp_spec_K = 273.15 + 20;
idx_T2_temp_corr = 8:12; % ignore long T2/diffusion and short T2/low SNR
idx_T1_temp_corr = 1:14;

make_rois = 0;
fn_roi = 'Ostenson_MRsys_rois.mat';

T1_spec_ms_v = [2480, 2173, 1907, 1604, 1332, 1044, 801.7, 608.6, 458.4, 336.5, 244.2, 176.6, 126.9, 90.9];
T2_spec_ms_v = [581.3, 403.5, 278.1, 190.9, 133.3, 96.9, 64.1, 46.4, 32.0, 22.6, 15.8, 11.2 7.9, 5.6];

n_ROI = numel(T1_spec_ms_v);

%% get data

n_fn = numel(fn_epg);

for ii = 1:n_fn
    
    load( [dir_in fn_epg{ii}] );
    T1_maps_epg(:,:,ii) = output_MRF_match_epg.T1_map;
    T2_maps_epg(:,:,ii) = output_MRF_match_epg.T2_map;
    
    M0_maps(:,:,ii) = abs( output_MRF_match_epg.M0_map );
    
    clear output_MRF_match_epg;
    
    load( [dir_in fn_pepg{ii}] );
    T1_maps_pepg(:,:,ii) = output_MRF_match_pepg.T1_map;
    T2_maps_pepg(:,:,ii) = output_MRF_match_pepg.T2_map;
    
    M0_maps(:,:,ii) = abs( output_MRF_match_pepg.M0_map );
    
    clear output_MRF_match_pepg;
    
    load( [dir_in fn_ssepg{ii}] );
    T1_maps_ssepg(:,:,ii) = output_MRF_match_ssepg.T1_map;
    T2_maps_ssepg(:,:,ii) = output_MRF_match_ssepg.T2_map;
    
    M0_maps(:,:,ii) = M0_maps(:,:,ii) + abs( output_MRF_match_ssepg.M0_map );
    
end

%% get rois

M0 = sum( M0_maps, 3 );
figure(1); clf;
imagesc( M0 ); axis image; colormap(gray);
if make_rois == 1
    for ii = 1:n_ROI
        figure(1);
        title( sprintf('Draw ROI %d',ii) );
        rois(:,:,ii) = roipoly();
    end
    save([dir_conv fn_roi],'rois');
else
    load([dir_conv fn_roi])
end

%% T1 and T2 mapping

if make_conv == 1
    
    load([dir_conv fn_conv])
    T1_data = T1_data(:,:,:,1);
    T2_data = T2_data(:,:,:,1);
    
    img_mask = sum(rois,3);
    
    % IR
    
    lb = [0 -1 1];
    ub = [1.2*max(abs(T1_data(:))) -0.2 5000];
    
    [T1_map_ms, RMSE_T1] = T1_IR_TD( T1_data, IR_times_ms_v, TD_ms, img_mask, lb, ub, 1 );
    
    figure(7); clf;
    imagesc( T1_map_ms ); axis image; colorbar();
    title('Conventional T1 map')
    drawnow;
    
    conventional_maps.SIR.T1_map_ms = T1_map_ms;
    conventional_maps.SIR.RMSE_T1 = RMSE_T1;
    
    
    % T2 mapping
    
    % MSE
    
    if n_T2_mod_params == 3
        lb = [0 0 0];
        ub = [5*max(T2_data(:)) 3000 0.5*max(T2_data(:))];
    elseif n_T2_mod_params == 2
        lb = [0 0];
        ub = [5*max(T2_data(:)) 3000];
    end
    [T2_map_ms, RMSE_T2] = T2_MSE( T2_data, TE_times_ms_v, img_mask, lb, ub, n_T2_mod_params, 1 );
    
    figure(8); clf;
    imagesc( T2_map_ms ); axis image; colorbar();
    title('Conventional T2 map')
    drawnow;
    
    conventional_maps.SE.T2_map_ms = T2_map_ms;
    conventional_maps.SE.RMSE_T2 = RMSE_T2;
    
    save([dir_out fn_conv_proc],'conventional_maps');
    
else
    
    load([dir_in fn_conv_proc]);
    
end

T1_ref = conventional_maps.SIR.T1_map_ms;
T2_ref = conventional_maps.SE.T2_map_ms;

%% get stats

load('../data_in/B0_MRsys.mat');

B0_mean = [];
B0_std = [];
for ii = 1:n_ROI
   
    my_roi = logical( rois(:,:,ii) );
    
    B0_mean(ii) = mean( B0_map(my_roi) );
    B0_std(ii) = std( B0_map(my_roi) );
    
    for jj = 1:n_fn
    
        my_T1 = T1_maps_epg(:,:,jj);
        T1_means_epg(ii,jj) = mean( my_T1( my_roi ) );
        my_T2 = T2_maps_epg(:,:,jj);
        T2_means_epg(ii,jj) = mean( my_T2( my_roi ) );
        
        my_T1 = T1_maps_pepg(:,:,jj);
        T1_means_pepg(ii,jj) = mean( my_T1( my_roi ) );
        my_T2 = T2_maps_pepg(:,:,jj);
        T2_means_pepg(ii,jj) = mean( my_T2( my_roi ) );
        
        my_T1 = T1_maps_ssepg(:,:,jj);
        T1_means_ssepg(ii,jj) = mean( my_T1( my_roi ) );
        my_T2 = T2_maps_ssepg(:,:,jj);
        T2_means_ssepg(ii,jj) = mean( my_T2( my_roi ) );
    
    end
    
    my_T1 = T1_ref;
    T1_median_conv(ii) = median( my_T1( my_roi ) );
    my_T2 = T2_ref;
    T2_median_conv(ii) = median( my_T2( my_roi ) );

end

%% get temp corr for spec based on conventional measurements

corr_temp_T1 = temp_K/temp_spec_K * pinv( T1_spec_ms_v(idx_T1_temp_corr).' ) * T1_median_conv(idx_T1_temp_corr).';
corr_temp_T2 = temp_K/temp_spec_K * pinv( T2_spec_ms_v(idx_T2_temp_corr).' ) * T2_median_conv(idx_T2_temp_corr).';

T1_spec_corr_v = corr_temp_T1 * T1_spec_ms_v;
T2_spec_corr_v = corr_temp_T2 * T2_spec_ms_v;

%% save stats

save( '../data_out/MRF_and_conv_MRsys_estimates.mat','T1_means_epg','T1_means_pepg','T1_means_ssepg','T1_spec_corr_v',...
    'T2_means_epg','T2_means_pepg','T2_means_ssepg','T2_spec_corr_v' );

%% CCC

for ii = 1:n_fn
    
    [CCC_T1_epg(ii), CCC_T1_ci_epg(ii,:)] = ccc( T1_means_epg(:,ii), T1_spec_corr_v );
    [CCC_T1_pepg(ii), CCC_T1_ci_pepg(ii,:)] = ccc( T1_means_pepg(:,ii), T1_spec_corr_v );
    [CCC_T1_ssepg(ii), CCC_T1_ci_ssepg(ii,:)] = ccc( T1_means_ssepg(:,ii), T1_spec_corr_v );
    
    [CCC_T2_epg(ii), CCC_T2_ci_epg(ii,:)] = ccc( T2_means_epg(:,ii), T2_spec_corr_v );
    [CCC_T2_pepg(ii), CCC_T2_ci_pepg(ii,:)] = ccc( T2_means_pepg(:,ii), T2_spec_corr_v );
    [CCC_T2_ssepg(ii), CCC_T2_ci_ssepg(ii,:)] = ccc( T2_means_ssepg(:,ii), T2_spec_corr_v );
    
end

fprintf('Mean CCC T1 -> EPG %.3f :: pEPG %.3f :: ssEPG %.3f\n', mean(CCC_T1_epg), mean(CCC_T1_pepg), mean(CCC_T1_ssepg) );
fprintf('Mean CCC T2 -> EPG %.3f :: pEPG %.3f :: ssEPG %.3f\n', mean(CCC_T2_epg), mean(CCC_T2_pepg), mean(CCC_T2_ssepg) );

%% plot results

figure(10); clf;
plot(T1_spec_corr_v,T1_means_epg,'x')
line([0 max(T1_spec_corr_v)],[0 max(T1_spec_corr_v)],'Color','r');
xlabel('Spec T1 (ms)'); ylabel('T1 (ms)')
title('T1 EPG')
legend( legend_txt, 'Location', 'northwest');

figure(11); clf;
plot(T1_spec_corr_v,T1_means_pepg,'x')
line([0 max(T1_spec_corr_v)],[0 max(T1_spec_corr_v)],'Color','r');
xlabel('Spec T1 (ms)'); ylabel('T1 (ms)')
title('T1 pEPG')
legend( legend_txt, 'Location', 'northwest');

figure(12); clf;
plot(T1_spec_corr_v,T1_means_ssepg,'x')
line([0 max(T1_spec_corr_v)],[0 max(T1_spec_corr_v)],'Color','r');
xlabel('Spec T1 (ms)'); ylabel('T1 (ms)')
title('T1 ssEPG')
legend( legend_txt, 'Location', 'northwest');

figure(13); clf;
plot(T1_spec_corr_v,T1_median_conv,'x')
line([0 max(T1_spec_corr_v)],[0 max(T1_spec_corr_v)],'Color','r');
xlabel('Spec T1 (ms)'); ylabel('T1 (ms)')
title('T1 conv')

figure(20); clf;
set(gcf,'Position',[200 -200 600 600],'Color','w')
plot(T2_spec_corr_v,T2_means_epg(:,1:2:end),'x')
hold on
plot(T2_spec_corr_v,T2_means_epg(:,2:2:end),'^')
hold off
grid on
xlim([0 800]); ylim([0 1100]);
line([0 max(T2_spec_corr_v)],[0 max(T2_spec_corr_v)],'Color','r');
xlabel('Spec T2 (ms)'); ylabel('T2 (ms)')
title('T2 EPG')
legend( {legend_txt{1:2:end}, legend_txt{2:2:end}}, 'Location', 'northwest');

figure(21); clf;
set(gcf,'Position',[200 -200 600 600],'Color','w')
plot(T2_spec_corr_v,T2_means_pepg(:,1:2:end),'x')
hold on
plot(T2_spec_corr_v,T2_means_pepg(:,2:2:end),'^')
hold off
grid on
xlim([0 800]); ylim([0 1100]);
line([0 max(T2_spec_corr_v)],[0 max(T2_spec_corr_v)],'Color','r');
xlabel('Spec T2 (ms)'); ylabel('T2 (ms)')
title('T2 pEPG')
legend( {legend_txt{1:2:end}, legend_txt{2:2:end}}, 'Location', 'northwest');

figure(22); clf;
set(gcf,'Position',[200 -200 600 600],'Color','w')
plot(T2_spec_corr_v,T2_means_ssepg(:,1:2:end),'x')
hold on
plot(T2_spec_corr_v,T2_means_ssepg(:,2:2:end),'^')
hold off
grid on
xlim([0 800]); ylim([0 1100]);
line([0 max(T2_spec_corr_v)],[0 max(T2_spec_corr_v)],'Color','r');
xlabel('Spec T2 (ms)'); ylabel('T2 (ms)')
title('T2 ssEPG')
legend( {legend_txt{1:2:end}, legend_txt{2:2:end}}, 'Location', 'northwest');

figure(23); clf;
plot(T2_spec_corr_v,T2_median_conv,'x')
line([0 max(T2_spec_corr_v)],[0 max(T2_spec_corr_v)],'Color','r');
xlabel('Spec T2 (ms)'); ylabel('T2 (ms)')
title('T2 conv')

%% plot CCC

c = categorical({'TBW4/Cr1','TBW4/Cr4','TBW8/Cr1','TBW8/Cr4'});

figure(30); clf;
set(gcf,'Position',[200 -200 600 600],'Color','w')
bar( c,[CCC_T1_epg(:) CCC_T1_pepg(:) CCC_T1_ssepg(:)] )
ylim([0 1.3])
ylabel('CCC')
title('T1 CCC')
legend('EPG','pEPG','ssEPG')

figure(31); clf;
set(gcf,'Position',[200 -200 600 600],'Color','w')
bar( c,[CCC_T2_epg(:) CCC_T2_pepg(:) CCC_T2_ssepg(:)] )
ylim([0 1.3])
ylabel('CCC')
title('T2 CCC')
legend('EPG','pEPG','ssEPG')

%% plot T2 over physilogical range

co = [0 0 1; 255/255 140/255 0];

figure(120); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
plot(T2_spec_corr_v(idx_T2_phys),T2_means_epg(idx_T2_phys,1:2:end),'x')
hold on
plot(T2_spec_corr_v(idx_T2_phys),T2_means_epg(idx_T2_phys,2:2:end),'^')
hold off
set(gca,'LineWidth',2);
grid on
xlim([0 100]); ylim([0 200]);
line([0 max(T2_spec_corr_v)],[0 max(T2_spec_corr_v)],'Color','r','LineStyle','--');
xlabel('Spec T2 (ms)'); ylabel('T2 (ms)')
title('EPG')
legend( {legend_txt{1:2:end}, legend_txt{2:2:end}}, 'Location', 'northwest');
legend boxoff;
export_fig('../figures/MRsys_T2_epg.png',sprintf('-r%d',4*my_dpi));

figure(121); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
plot(T2_spec_corr_v(idx_T2_phys),T2_means_pepg(idx_T2_phys,1:2:end),'x')
hold on
plot(T2_spec_corr_v(idx_T2_phys),T2_means_pepg(idx_T2_phys,2:2:end),'^')
hold off
set(gca,'LineWidth',2);
grid on
xlim([0 100]); ylim([0 200]);
line([0 max(T2_spec_corr_v)],[0 max(T2_spec_corr_v)],'Color','r','LineStyle','--');
xlabel('Spec T2 (ms)'); ylabel('T2 (ms)')
title('pEPG')
legend( {legend_txt{1:2:end}, legend_txt{2:2:end}}, 'Location', 'northwest');
legend boxoff;
export_fig('../figures/MRsys_T2_pepg.png',sprintf('-r%d',4*my_dpi));

figure(122); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
plot(T2_spec_corr_v(idx_T2_phys),T2_means_ssepg(idx_T2_phys,1:2:end),'x')
hold on
plot(T2_spec_corr_v(idx_T2_phys),T2_means_ssepg(idx_T2_phys,2:2:end),'^')
hold off
set(gca,'LineWidth',2);
grid on
xlim([0 100]); ylim([0 200]);
line([0 max(T2_spec_corr_v)],[0 max(T2_spec_corr_v)],'Color','r','LineStyle','--');
xlabel('Spec T2 (ms)'); ylabel('T2 (ms)')
title('ssEPG')
legend( {legend_txt{1:2:end}, legend_txt{2:2:end}}, 'Location', 'northwest');
legend boxoff;
export_fig('../figures/MRsys_T2_ssepg.png',sprintf('-r%d',4*my_dpi));

%% plot CCC over physiological range of T2

c = {'TBW4/C1','TBW4/C4','TBW8/C1','TBW8/C4'};

co = [0 0 0; 0.6 0.6 0; 0 0.75 0.75];

for ii = 1:n_fn
    
    [CCC_physT2_epg(ii), CCC_physT2_ci_epg(ii,:)] = ccc( T2_means_epg(idx_T2_phys,ii), T2_spec_corr_v(idx_T2_phys) );
    [CCC_physT2_pepg(ii), CCC_physT2_ci_pepg(ii,:)] = ccc( T2_means_pepg(idx_T2_phys,ii), T2_spec_corr_v(idx_T2_phys) );
    [CCC_physT2_ssepg(ii), CCC_physT2_ci_ssepg(ii,:)] = ccc( T2_means_ssepg(idx_T2_phys,ii), T2_spec_corr_v(idx_T2_phys) );
    
end

ccc_cat = [CCC_physT2_ssepg(:) CCC_physT2_pepg(:) CCC_physT2_epg(:)];
ccc_ci_cat_lo = ccc_cat - [CCC_physT2_ci_ssepg(:,1) CCC_physT2_ci_pepg(:,1) CCC_physT2_ci_epg(:,1)];
ccc_ci_cat_hi = [CCC_physT2_ci_ssepg(:,2) CCC_physT2_ci_pepg(:,2) CCC_physT2_ci_epg(:,2)] - ccc_cat;
figure(131); clf;
set(gcf,'Position',[200 100 600 500],'Color','w')
set(gcf,'defaultAxesColorOrder',co)
b = bar( ccc_cat, 'FaceColor','w','LineWidth',2 );
b(1).EdgeColor = co(1,:);
b(2).EdgeColor = co(2,:);
b(3).EdgeColor = co(3,:);
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(ccc_cat, 1);
nbars = size(ccc_cat, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    errorbar(x, ccc_cat(:,ii), ccc_ci_cat_lo(:,ii), ccc_ci_cat_hi(:,ii), 'k', 'linestyle', 'none','LineWidth',2);
end
ylim([0 1.3])
ylabel('CCC')
set(gca,'XTickLabel',c,'FontSize',18);
% title('T2 CCC')
my_leg = legend('ssEPG','pEPG','EPG','Location','northwest');
my_leg.Position(2) = 0.75;
legend boxoff;
export_fig('../figures/MRsys_T2_ccc.png',sprintf('-r%d',4*my_dpi));



