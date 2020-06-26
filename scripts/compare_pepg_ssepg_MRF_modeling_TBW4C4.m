%% Compare pEPG against ssEPG for T2 estimation using TBW4C4, MR system phantom

%% params

fn_D_pepg = '../data_out/MRF_dict_pepg_105_FA60_TBW4.mat';
fn_D_ssepg = '../data_out/MRF_dict_105_FA60_TBW4_crush4.mat';

B1 = 1.25;

%% get pEPG and ssEPG dictionaries

load( fn_D_pepg );
dict_pepg = output_dict;
clear output_dict;

load( fn_D_ssepg );
dict_ssepg = output_dict;
clear output_dict;

%% get sub-dictionaries based on MR system ~B1

idx_B1 = abs( dict_pepg.dict_list(:,3) - B1 ) <= eps; % get B1

if any( idx_B1 - ( abs( dict_ssepg.dict_list(:,3) - B1 ) <= eps ) )
    error('Dictionary index mismatch.')
end

D_pepg = conj( dict_pepg.dict_norm(:,idx_B1) ); %reconcile sign convention between pepg and ssepg dictionaries
D_ssepg = dict_ssepg.dict_norm(:,idx_B1);

T1T2_list = dict_pepg.dict_list(idx_B1,1:2);

%% Fit T1 and T2 using pEPG with ssEPG serving as measured signal

Dcc = D_ssepg' * D_pepg;
[Max_ips,Idx_ips] = max( Dcc,[],2 );
T1T2_fits = T1T2_list( Idx_ips,: );


%% get stats and interpolate

relative_diff = (T1T2_fits - T1T2_list)./T1T2_list;

figure(1); clf;
plot( relative_diff );
xlabel('Dictionary index')
ylabel('relative diff')
legend('T1','T2');

S_T1 = scatteredInterpolant( T1T2_list(:,1),T1T2_list(:,2),relative_diff(:,1) );
S_T2 = scatteredInterpolant( T1T2_list(:,1),T1T2_list(:,2),relative_diff(:,2) );

t1_q = (100:10:2900)';
t2_q = (2:1:650)';
[T1_q,T2_q] = meshgrid( t1_q, t2_q );
Relative_diff_T1_grid = S_T1( T1_q,T2_q );
Relative_diff_T2_grid = S_T2( T1_q,T2_q );

idx_zero = T1_q < T2_q;
Relative_diff_T1_grid( idx_zero ) = 0;
Relative_diff_T2_grid( idx_zero ) = 0;

T1_ticks = [500:500:2500];
T2_ticks = [5 50:100:650];

[~,idx_T1_ticks] = min(abs(t1_q(:) - T1_ticks),[],1);
[~,idx_T2_ticks] = min(abs(t2_q(:) - T2_ticks),[],1);

T1_tick_labels = num2cell( t1_q( idx_T1_ticks ) );
T2_tick_labels = num2cell( t2_q( idx_T2_ticks ) );

figure(2); clf;
set(gcf,'Position',[200 0 800 600],'Color','w');
imagesc( Relative_diff_T1_grid ); colormap(flipud(hot)); colorbar();
set( gca,'XTick',idx_T1_ticks,'XTickLabel',T1_tick_labels );
set( gca,'YTick',idx_T2_ticks,'YTickLabel',T2_tick_labels );
xlabel('T_1 (ms)')
ylabel('T_2 (ms)')

figure(3); clf;
set(gcf,'Position',[200 200 900 800],'Color','w');
imagesc( Relative_diff_T2_grid ); colormap(flipud(hot));
c = colorbar();
decrease_by = 0.05;
axpos = get(gca,'position');
axpos(3) = axpos(3) - decrease_by;
set(gca,'position',axpos);
c.Label.String = '% difference in T_2';
c.Position(1) = 0.82;
set( gca,'XTick',idx_T1_ticks,'XTickLabel',T1_tick_labels );
set( gca,'YTick',idx_T2_ticks,'YTickLabel',T2_tick_labels );
xlabel('T_1 (ms)')
ylabel('T_2 (ms)')

%% plot results

load('../data_out/MRF_and_conv_MRsys_estimates.mat');
markers = {'+','o','*','.','x','s','d','^','v','>','<','o','h'};

T2_pepg_mrsys = T2_means_pepg(:,2);
T1_ssepg_mrsys = T1_means_ssepg(:,2);
T2_ssepg_mrsys = T2_means_ssepg(:,2);

Esimated_pepg_ssepg_relive_diff = (T2_pepg_mrsys(:) - T2_ssepg_mrsys(:))./T2_ssepg_mrsys(:);
Predicted_pepg_T2_diff = S_T2( T1_ssepg_mrsys,T2_ssepg_mrsys );

figure(10); clf;
set(gcf,'Position',[200 200 900 800],'Color','w');
for ii = 1:numel(Predicted_pepg_T2_diff)
    plot(Predicted_pepg_T2_diff(ii)*100,Esimated_pepg_ssepg_relive_diff(ii)*100,markers{mod(ii-1,numel(markers))+1},'MarkerSize',10,'Color',[0,0,255]/255)
    hold on;
end
xlim([-20 1]); ylim([-20 1]); grid on;
xlabel('Modeled relative difference in T_2 (%)')
ylabel('Estimated relative difference in T_2 (%)')
export_fig( '../figures/diff_T2_pepg_v_ssepg_mrsys.png','-r400' );


my_color = [176,224,230]/255;
imh  = figure(3);
hold on;
for ii = 1:numel(Predicted_pepg_T2_diff)
    my_t1 = T1_ssepg_mrsys(ii);
    my_t2 = T2_ssepg_mrsys(ii);
    my_idx_t1 = round((my_t1 - t1_q(1))/(t1_q(end)-t1_q(1))*numel(t1_q));
    if my_idx_t1 < 1
        my_idx_t1 = 1;
    end
    my_idx_t2 = round((my_t2 - t2_q(1))/(t2_q(end)-t2_q(1))*numel(t2_q));
    my_pos = [my_idx_t1 my_idx_t2];
    plot(my_idx_t1,my_idx_t2,'marker',markers{ mod(ii-1,numel(markers))+1},'MarkerSize',10,'Color',my_color)
end
hold off;
export_fig( '../figures/diff_T2_pepg_v_ssepg_model_with_pnts.png','-r400' );

%% generate unified figure

im1 = imread( '../figures/diff_T2_pepg_v_ssepg_model_with_pnts.png' );
im2 = imread( '../figures/diff_T2_pepg_v_ssepg_mrsys.png' );

im2 = imresize( im2, [size(im1,1) size(im1,2)] );

buff_im = 255 * ones( size(im1,1), 150, 3 );

imtot = cat( 2,im1,buff_im, im2 );

figure(10); clf;
imshow( imtot, 'Border', 'tight' );
text(70,170,'a','FontSize',34);
text(3500,170,'b','FontSize',34);
export_fig( '../figures/diff_T2_ssepg_on_pepg.png','-r400');