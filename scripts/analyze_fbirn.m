%% FBIRN B0 fitting and analysis

font_sz = 20;
my_dpi = 250;

%% params

dir_in = '../data_out/';

row_start = 1;
row_end = 240;
col_start = 1;
col_end = 240;

fn_epg{1} = 'MRF_fbirn_TBW4_C1_proc_epg.mat';
fn_epg{2} = 'MRF_fbirn_TBW4_C2_proc_epg.mat';
fn_epg{3} = 'MRF_fbirn_TBW4_C4_proc_epg.mat';

fn_pepg{1} = 'MRF_fbirn_TBW4_C1_proc_pepg.mat';
fn_pepg{2} = 'MRF_fbirn_TBW4_C2_proc_pepg.mat';
fn_pepg{3} = 'MRF_fbirn_TBW4_C4_proc_pepg.mat';

fn_ssepg{1} = 'MRF_fbirn_TBW4_C1_proc_ssepg.mat';
fn_ssepg{2} = 'MRF_fbirn_TBW4_C2_proc_ssepg.mat';
fn_ssepg{3} = 'MRF_fbirn_TBW4_C4_proc_ssepg.mat';

n_fn = numel( fn_epg );

fn_img{1} = 'MRF_fbirn_TBW4_C1_img_recon.mat';
fn_img{2} = 'MRF_fbirn_TBW4_C2_img_recon.mat';
fn_img{3} = 'MRF_fbirn_TBW4_C4_img_recon.mat';

fn_dict{1} = 'MRF_dict_leg_B0_TBW4_crush1.mat';
fn_dict{2} = 'MRF_dict_leg_B0_TBW4_crush2.mat';
fn_dict{3} = 'MRF_dict_leg_B0_TBW4_crush4.mat';

%% get B0 corrected maps

for ii = 1:n_fn
    
    load([dir_in fn_img{ii}]);
    D = load([dir_in fn_dict{ii}]);
    D.dict_compress = [D.dict_compress01 D.dict_compress02];
    
    match = MRF_dict_match_B0( output_img_recon, D );
    
    eval( sprintf( 'match%.2d = match;', ii ) );
    clear match D;
    
end

%%

n_row = row_end - row_start + 1;
n_col = col_end - col_start + 1;

for ii = 1:n_fn
    
    load( [dir_in fn_epg{ii}] );
    if ii == 1
        T1_epg = zeros( n_row, n_col, n_fn );
        T2_epg = T1_epg;
        T1_pepg = T1_epg;
        T2_pepg = T1_epg;
        T1_ssepg = T1_epg;
        T2_ssepg = T1_epg;
        T1_ssepgB0 = T1_epg;
        T2_ssepgB0 = T1_epg;
        mask = zeros( n_row, n_col );
        my_M0 = output_MRF_match_epg.M0_map(row_start:row_end,col_start:col_end);
        mask( abs( my_M0 ) > multithresh( abs( my_M0 ), 1 ) ) = 1;
        mask( output_MRF_match_epg.T2_map(row_start:row_end,col_start:col_end) > 150 ) = 0;
        mask( output_MRF_match_epg.T1_map(row_start:row_end,col_start:col_end) < 200 ) = 0;
        mask = bwareafilt( logical(mask),1 );
        
    end
    T1_epg(:,:,ii) = mask.*output_MRF_match_epg.T1_map(row_start:row_end,col_start:col_end);
    T2_epg(:,:,ii) = mask.*output_MRF_match_epg.T2_map(row_start:row_end,col_start:col_end);
    clear output_MRF_match_epg;
    
    load( [dir_in fn_pepg{ii}] );
    T1_pepg(:,:,ii) = mask.*output_MRF_match_pepg.T1_map(row_start:row_end,col_start:col_end);
    T2_pepg(:,:,ii) = mask.*output_MRF_match_pepg.T2_map(row_start:row_end,col_start:col_end);
    clear output_MRF_match_pepg;
   
    load( [dir_in fn_ssepg{ii}] );
    T1_ssepg(:,:,ii) = mask.*output_MRF_match_ssepg.T1_map(row_start:row_end,col_start:col_end);
    T2_ssepg(:,:,ii) = mask.*output_MRF_match_ssepg.T2_map(row_start:row_end,col_start:col_end);
    clear output_MRF_match_ssepg;
    
    eval( sprintf( 'T1_ssepgB0(:,:,ii) = mask.*match%.2d.T1_map(row_start:row_end,col_start:col_end);',ii ) );
    eval( sprintf( 'T2_ssepgB0(:,:,ii) = mask.*match%.2d.T2_map(row_start:row_end,col_start:col_end);',ii ) );
    eval( sprintf( 'B0_ssepgB0(:,:,ii) = mask.*match%.2d.B0_map(row_start:row_end,col_start:col_end);',ii ) );
    
end

%%

std_T2_epg = std( T2_epg,[],3 );
std_T2_pepg = std( T2_pepg,[],3 );
std_T2_ssepg = std( T2_ssepg,[],3 );
std_T2_ssepgB0 = std( T2_ssepgB0,[],3 );

mean_T2_epg = mean( T2_epg,3 );
mean_T2_pepg = mean( T2_pepg,3 );
mean_T2_ssepg = mean( T2_ssepg,3 );
mean_T2_ssepgB0 = mean( T2_ssepgB0,3 );

mean_B0 = mean( B0_ssepgB0,3 );
mean_B0_zoned = mean_B0;
mean_B0_zoned( abs(mean_B0) < 0.70 * 1/16/2 | abs(mean_B0) > 1.30 * 1/16/2 ) = 0;

%%

n_vox = sum( mask(:) );
T1s_epg = zeros( n_vox, 3 );
T1s_pepg = T1s_epg;
T1s_ssepg = T1s_epg;
T1s_ssepgB0 = T1s_epg;
T2s_epg = T1s_epg;
T2s_pepg = T1s_epg;
T2s_ssepg = T1s_epg;
T2s_ssepgB0 = T1s_epg;
for ii = 1:n_fn
    
    my_T1 = T1_epg(:,:,ii);
    my_T2 = T2_epg(:,:,ii);
    T1s_epg(:,ii) = my_T1( mask );
    T2s_epg(:,ii) = my_T2( mask );
    
    my_T1 = T1_pepg(:,:,ii);
    my_T2 = T2_pepg(:,:,ii);
    T1s_pepg(:,ii) = my_T1( mask );
    T2s_pepg(:,ii) = my_T2( mask );
    
    my_T1 = T1_ssepg(:,:,ii);
    my_T2 = T2_ssepg(:,:,ii);
    T1s_ssepg(:,ii) = my_T1( mask );
    T2s_ssepg(:,ii) = my_T2( mask );
    
    my_T1 = T1_ssepgB0(:,:,ii);
    my_T2 = T2_ssepgB0(:,:,ii);
    T1s_ssepgB0(:,ii) = my_T1( mask );
    T2s_ssepgB0(:,ii) = my_T2( mask );

end

means_T1_epg = mean( T1s_epg,1 );
means_T2_epg = mean( T2s_epg,1 );
stds_T1_epg = std( T1s_epg,[],1 );
stds_T2_epg = std( T2s_epg,[],1 );
ses_T1_epg = stds_T1_epg/sqrt(n_vox);
ses_T2_epg = stds_T2_epg/sqrt(n_vox);

means_T1_pepg = mean( T1s_pepg,1 );
means_T2_pepg = mean( T2s_pepg,1 );
stds_T1_pepg = std( T1s_pepg,[],1 );
stds_T2_pepg = std( T2s_pepg,[],1 );
ses_T1_pepg = stds_T1_pepg/sqrt(n_vox);
ses_T2_pepg = stds_T2_pepg/sqrt(n_vox);

means_T1_ssepg = mean( T1s_ssepg,1 );
means_T2_ssepg = mean( T2s_ssepg,1 );
stds_T1_ssepg = std( T1s_ssepg,[],1 );
stds_T2_ssepg = std( T2s_ssepg,[],1 );
ses_T1_ssepg = stds_T1_ssepg/sqrt(n_vox);
ses_T2_ssepg = stds_T2_ssepg/sqrt(n_vox);

means_T1_ssepgB0 = mean( T1s_ssepgB0,1 );
means_T2_ssepgB0 = mean( T2s_ssepgB0,1 );
stds_T1_ssepgB0 = std( T1s_ssepgB0,[],1 );
stds_T2_ssepgB0 = std( T2s_ssepgB0,[],1 );
ses_T1_ssepgB0 = stds_T1_ssepgB0/sqrt(n_vox);
ses_T2_ssepgB0 = stds_T2_ssepgB0/sqrt(n_vox);

means_T1 = [means_T1_epg; means_T1_pepg; means_T1_ssepg; means_T1_ssepgB0];
means_T2 = [means_T2_epg; means_T2_pepg; means_T2_ssepg; means_T2_ssepgB0];

ses_T1 = [ses_T1_epg; ses_T1_pepg; ses_T1_ssepg; ses_T1_ssepgB0];
ses_T2 = [ses_T2_epg; ses_T2_pepg; ses_T2_ssepg; ses_T2_ssepgB0];

fprintf('Mean T1 TBW4 C4 (ms): EPG %0.f :: pEPG %0.f :: ssEPG %0.f :: ssEPGB0 %0.f\n', means_T1(:,3) );
fprintf('Mean T2 TBW4 C4 (ms): EPG %0.f :: pEPG %0.f :: ssEPG %0.f :: ssEPGB0 %0.f\n', means_T2(:,3) );

%%

[n_r, n_c, ~] = size( T1_epg );

T1_mosaic = [reshape( T1_epg, [n_r n_c*n_fn] ); reshape( T1_pepg, [n_r n_c*n_fn] ); ...
    reshape( T1_ssepg, [n_r n_c*n_fn] ); reshape( T1_ssepgB0, [n_r n_c*n_fn] ) ];
T1_mosaic = [T1_mosaic, zeros(size(T1_mosaic,1),100)];

figure(1); clf;
set(gcf,'Position',[200 100 900 900],'Color','w')
imagesc( T1_mosaic );
axis image; axis off; 
c = colorbar('East','Color','w','FontSize',font_sz);
c.Label.String = 'T_1 (ms)';
c.Position(2) = 0.25;
c.Position(4) = 0.5;
c.Label.Position(1) = -2;
c.Label.Position(2) = 500;
c.LineWidth = 2;
caxis([0 1000]); colormap(hot);
export_fig('../figures/fbirn_T1.png',sprintf('-r%d',my_dpi));

T2_mosaic = [reshape( T2_epg, [n_r n_c*n_fn] ); reshape( T2_pepg, [n_r n_c*n_fn] ); ...
    reshape( T2_ssepg, [n_r n_c*n_fn] ); reshape( T2_ssepgB0, [n_r n_c*n_fn] ) ];
T2_mosaic = [T2_mosaic, zeros(size(T2_mosaic,1),100)];

figure(2); clf;
set(gcf,'Position',[200 100 900 900],'Color','w')
imagesc( T2_mosaic );
axis image; axis off; 
c = colorbar('East','Color','w','FontSize',font_sz);
c.Label.String = 'T_2 (ms)';
c.Position(2) = 0.25;
c.Position(4) = 0.5;
c.Label.Position(1) = -2;
c.Label.Position(2) = 50;
c.LineWidth = 2;
caxis([0 100]); colormap(hot);
export_fig('../figures/fbirn_T2.png',sprintf('-r%d',my_dpi));


figure(10); clf;
set(gcf,'Position',[200 100 550 500],'Color','w')
imagesc( mean_B0*1000 );
axis image; axis off;
c = colorbar('North','Color','w','FontSize',font_sz);
c.Label.String = '\DeltaB_0 (Hz)';
caxis([-45 45]); colormap(gray);
c.Label.Position(2) = -1.3;
export_fig('../figures/fbirn_B0.png',sprintf('-r%d',my_dpi));



figure(11); clf;
set(gcf,'Position',[200 100 550 500],'Color','w')
imagesc( mean_B0*1000 ); axis image; axis off; colormap(gray); caxis([-45 45]);
red = cat(3, ones( size(mask) ), zeros( size(mask) ), zeros( size(mask) ) );
hold on
h = imshow( red );
hold off
set( h,'AlphaData',abs(mean_B0_zoned*1000) );
export_fig('../figures/fbirn_B0_zoned.png',sprintf('-r%d',my_dpi));


COV_mosaic = 100*[ std_T2_epg./mean_T2_epg; std_T2_pepg./mean_T2_pepg; std_T2_ssepg./mean_T2_ssepg; std_T2_ssepgB0./mean_T2_ssepgB0 ];
COV_mosaic = [ COV_mosaic, zeros(size(T2_mosaic,1),75) ];

figure(12); clf;
set(gcf,'Position',[200 100 900 900],'Color','w')
imagesc( COV_mosaic  );
axis image; axis off;
c = colorbar('East','Color','w','FontSize',font_sz);
c.Label.String = 'COV (%)';
c.Position(2) = 0.25;
c.Position(4) = 0.5;
c.LineWidth = 2;
caxis([0 50]); 
c.Label.Position(1) = -1.7;
export_fig('../figures/fbirn_COV.png',sprintf('-r%d',my_dpi));

