%% Figure unification/composite generation

dir_out = '../figures/';
my_dpi = 250;

%% Figures 1-2 already generated

%% Figure 3

font_sz = 28;
buff_sz = 250;

fn_3a = 'Figure3a.png';
fn_3b = 'Figure3b.png';

f3a = imread( [dir_out fn_3a] );
f3b = imread( [dir_out fn_3b] );

my_buff = 255*ones(buff_sz,size(f3a,2),3);

fig3 = [my_buff; f3a; my_buff; f3b];

figure(1); clf;
imshow( fig3,'Border','tight' );
txta = text(750,150,'One cycle per nominal slice thickness','FontSize',font_sz);
txtb = text(750,1850,'Four cycles per nominal slice thickness','FontSize',font_sz);
export_fig([dir_out 'Figure3.png'],sprintf('-r%d',my_dpi));

%% Figure 4

font_sz = 28;
buff_sz = 100;

fn_4a = 'MRsys_T2_epg.png';
fn_4b = 'MRsys_T2_pepg.png';
fn_4c = 'MRsys_T2_ssepg.png';
fn_4d = 'MRsys_T2_ccc.png';

f4a = imread( [dir_out fn_4a] );
f4b = imread( [dir_out fn_4b] );
f4c = imread( [dir_out fn_4c] );
f4d = imread( [dir_out fn_4d] );

my_buff1 = 255 * ones( size(f4c,1) - size(f4d,1), size(f4d,2), 3 );
my_buff2 = 255 * ones( size(f4c,1), size(f4b,2) - size(f4d,2), 3 );
f4dmod = [my_buff2 cat(1,f4d,my_buff1) ];

top_row = cat(2,f4a,f4b);
bottom_row = cat(2,f4c,f4dmod);
my_buff = 255 * ones( buff_sz, size(top_row,2), 3 );
fig4 = [ top_row; my_buff; bottom_row ];

figure(4);
imshow(fig4,'Border','tight');
txta = text(50,175,'a','FontSize',font_sz);
txtb = text(50+size(f4a,2),175,'b','FontSize',font_sz);
txtc = text(50,175+size(f4a,1)+buff_sz,'c','FontSize',font_sz);
txtd = text(50+size(f4a,2),175+size(f4a,1)+buff_sz,'d','FontSize',font_sz);
export_fig([dir_out 'Figure4.png'],sprintf('-r%d',my_dpi));

%% Figure 5

font_sz = 28;
buff_sz = 100;

fn_5a = 'profiles_with_B0.png';
fn_5b = 'mag_with_B0.png';
fn_5c = 'phase_with_B0.png';

f5a = imread( [dir_out fn_5a] );
f5b = imread( [dir_out fn_5b] );
f5c = imread( [dir_out fn_5c] );

buffa = size( f5b,2 ) - size( f5a,2 );
buffc = size( f5b,2 ) - size( f5c,2 );
my_buff = 255*ones( buff_sz, size(f5b,2), 3 );

fig5 = [ [f5a 255*ones(size(f5a,1),buffa,3)]; my_buff; f5b; my_buff; [f5c 255*ones(size(f5c,1),buffc,3)] ];

figure(5); clf;
imshow(fig5,'Border','tight');
txta = text(400,75,'a','FontSize',font_sz);
txtb = text(400,buff_sz+100+size(f5a,1),'b','FontSize',font_sz);
txtc = text(400,2*buff_sz+75+size(f5a,1)+size(f5b,1),'c','FontSize',font_sz);
export_fig([dir_out 'Figure5.png'],sprintf('-r%d',4*my_dpi));

%% Figure 6

font_sz = 28;
buff_sz = 100;

fn_6a = 'b0_T2_bias_tbw4c1.png';
fn_6b = 'b0_T2_bias_tbw4c4.png';
fn_6c = 'b0_T2_bias_tbw4c8.png';
fn_6d = 'b0_T2_bias_tbw8c1.png';
fn_6e = 'b0_T2_bias_tbw8c4.png';
fn_6f = 'b0_T2_bias_tbw8c8.png';

f6a = imread( [dir_out fn_6a] );
f6b = imread( [dir_out fn_6b] );
f6c = imread( [dir_out fn_6c] );
f6d = imread( [dir_out fn_6d] );
f6e = imread( [dir_out fn_6e] );
f6f = imread( [dir_out fn_6f] );

top_row = [ f6a 255*ones(size(f6a,1),buff_sz,3) f6b 255*ones(size(f6a,1),buff_sz,3) f6c];
buff_row = 255*ones( buff_sz, size(top_row,2), 3 );
bottom_row = [ f6d 255*ones(size(f6a,1),buff_sz,3) f6e 255*ones(size(f6a,1),buff_sz,3) f6f ];

fig6 = [top_row; buff_row; bottom_row];
figure(6); clf;
imshow( fig6,'Border','tight')
export_fig([dir_out 'Figure6.png'],sprintf('-r%d',2*my_dpi));

%% Figure 7

font_sz = 28;
buff_sz = 150;

fn_7a = 'leg_T1.png';
fn_7b = 'leg_T2.png';

f7a = imread( [dir_out fn_7a] );
f7b = imread( [dir_out fn_7b] );

f7amod = [ 255*ones(buff_sz,buff_sz+size(f7a,2),3); [255*ones( size(f7a,1), buff_sz,3 ) f7a] ];
figure(99); clf;
imshow( f7amod,'Border','tight' );
txt1 = text(75,400,'EPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt2 = text(75,490+400,'pEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt3 = text(75,2*490+400,'ssEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt4 = text(75,90+3*490+400,'ssEPG w/B0','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt5 = text(250,75,'Crush. 1','FontSize',font_sz,'FontWeight','bold');
txt6 = text(250+510,75,'Crush. 2','FontSize',font_sz,'FontWeight','bold');
txt6 = text(250+2*510,75,'Crush. 4','FontSize',font_sz,'FontWeight','bold');
export_fig([dir_out 'tmp_7a.png'],sprintf('-r%d',my_dpi));
f7a = imread([dir_out 'tmp_7a.png']);

figure(99); clf;
f7bmod = [ 255*ones(buff_sz,buff_sz+size(f7b,2),3); [255*ones( size(f7b,1), buff_sz,3 ) f7b] ];
figure(99); clf;
imshow( f7bmod,'Border','tight' );
txt1 = text(75,400,'EPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt2 = text(75,490+400,'pEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt3 = text(75,2*490+400,'ssEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt4 = text(75,90+3*490+400,'ssEPG w/B0','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt5 = text(250,75,'Crush. 1','FontSize',font_sz,'FontWeight','bold');
txt6 = text(250+510,75,'Crush. 2','FontSize',font_sz,'FontWeight','bold');
txt6 = text(250+2*510,75,'Crush. 4','FontSize',font_sz,'FontWeight','bold');
export_fig([dir_out 'tmp_7b.png'],sprintf('-r%d',my_dpi));
f7b = imread([dir_out 'tmp_7b.png']);

fig7 = [f7a f7b];
figure(7); clf;
imshow( fig7,'Border','tight' );
txta = text( 225, 270, 'a', 'FontSize', 1.5*font_sz, 'Color','w' );
txtb = text( 225 + size(f7a,2), 270, 'b', 'FontSize', 1.5*font_sz, 'Color','w' );
export_fig([dir_out 'Figure7.png'],sprintf('-r%d',2*my_dpi));

%% Figure 8

font_sz = 28;
buff_sz = 150;

fn_8a = 'leg_COV.png';
fn_8b = 'leg_B0.png';
fn_8c = 'leg_B0_zoned.png';

f8a = imread( [dir_out fn_8a] );
f8b = imread( [dir_out fn_8b] );
f8b = repmat( f8b, [1 1 3] );
f8b(1:5,:,:) = 129*ones(5,size(f8b,2),3);
f8c = imread( [dir_out fn_8c] );
f8c(end-4:end,:,:) = 129*ones(5,size(f8b,2),3);

f8amod = [ 255*ones( size(f8a,1), buff_sz, 3 ) f8a ];

figure(99); clf;
imshow( f8amod,'Border','tight' );
txt1 = text(75,300,'EPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt2 = text(75,490+300,'pEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt3 = text(75,2*490+300,'ssEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt4 = text(75,100+3*490+300,'ssEPG w/B0','FontSize',font_sz,'FontWeight','bold','Rotation',90);
export_fig([dir_out 'tmp_8a.png'],sprintf('-r%d',my_dpi));

f8amod = imread( [dir_out 'tmp_8a.png'] );

buffbc = size( f8amod,1 ) - (size(f8b,1) + size(f8c,1));
buffbc1 = round( buffbc/3 );
buffbc2 = buffbc1;
buffbc3 = buffbc - (buffbc1 + buffbc2);

fig8 = [ f8amod [129*ones( buffbc1, size(f8b,2),3 ); f8b; 255*ones( buffbc2, size(f8b,2),3 ); f8c; 129*ones( buffbc3, size(f8b,2),3 )] ];

figure(8); clf;
imshow(fig8,'Border','tight');
txta = text(225,75,'a','FontSize',1.5*font_sz,'Color','w');
txtb = text(75+size(f8amod,2),75,'b','FontSize',1.5*font_sz,'Color','w');
txtc = text(75+size(f8amod,2),75+size(f8b,1)+buffbc1+buffbc2,'c','FontSize',1.5*font_sz,'Color','w');
export_fig([dir_out 'Figure8.png'],sprintf('-r%d',2*my_dpi));

%% supporting information detail of simulation and measurement (Figs. S1 and S2)

font_sz = 56;
buff_sz = 50;

fn_sim_det_a = 'numerical_sim_detail_03.png';
fn_sim_det_b = 'numerical_sim_detail_08.png';

im_a = imread( [dir_out fn_sim_det_a] );
im_b = imread( [dir_out fn_sim_det_b] );

n_r = size( im_a,1 );

SI_sim_det = cat(2, im_a, 255*ones( n_r,buff_sz,3 ), im_b );

figure(19); clf;
imshow( SI_sim_det, 'Border', 'tight'  );
txta = text( 250, 125, 'a', 'FontSize', font_sz );
txtb = text( 250+1350, 125, 'b', 'FontSize', font_sz );
export_fig([dir_out 'Figure_sim_detail.png'],sprintf('-r%d',2*my_dpi));

fn_msr_det_a = 'measured_profile_03.png';
fn_msr_det_b = 'measured_profile_08.png';

im_a = imread( [dir_out fn_msr_det_a] );
im_b = imread( [dir_out fn_msr_det_b] );

n_r = size( im_a,1 );

SI_msr_det = cat(2, im_a, 255*ones( n_r,buff_sz,3 ), im_b );

figure(20); clf;
imshow( SI_msr_det, 'Border', 'tight'  );
txta = text( 300, 125, 'a', 'FontSize', font_sz );
txtb = text( 300+1350, 125, 'b', 'FontSize', font_sz );
export_fig([dir_out 'Figure_msr_detail.png'],sprintf('-r%d',2*my_dpi));

%% Figures S3 and S4 already generated

%% supporting information FBIRN under het. B0: T1 and T2 (Fig. S5)

font_sz = 28;
buff_sz = 150;

fn_7a = 'fbirn_T1.png';
fn_7b = 'fbirn_T2.png';

f7a = imread( [dir_out fn_7a] );
f7b = imread( [dir_out fn_7b] );

f7amod = [ 255*ones(buff_sz,buff_sz+size(f7a,2),3); [255*ones( size(f7a,1), buff_sz,3 ) f7a] ];
figure(99); clf;
imshow( f7amod,'Border','tight' );
txt1 = text(75,400,'EPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt2 = text(75,490+400,'pEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt3 = text(75,2*490+400,'ssEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt4 = text(75,90+3*490+400,'ssEPG w/B0','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt5 = text(250,75,'Crush. 1','FontSize',font_sz,'FontWeight','bold');
txt6 = text(250+510,75,'Crush. 2','FontSize',font_sz,'FontWeight','bold');
txt6 = text(250+2*510,75,'Crush. 4','FontSize',font_sz,'FontWeight','bold');
export_fig([dir_out 'tmp_fbirn_a.png'],sprintf('-r%d',my_dpi));
f7a = imread([dir_out 'tmp_fbirn_a.png']);

figure(99); clf;
f7bmod = [ 255*ones(buff_sz,buff_sz+size(f7b,2),3); [255*ones( size(f7b,1), buff_sz,3 ) f7b] ];
figure(99); clf;
imshow( f7bmod,'Border','tight' );
txt1 = text(75,400,'EPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt2 = text(75,490+400,'pEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt3 = text(75,2*490+400,'ssEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt4 = text(75,90+3*490+400,'ssEPG w/B0','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt5 = text(250,75,'Crush. 1','FontSize',font_sz,'FontWeight','bold');
txt6 = text(250+510,75,'Crush. 2','FontSize',font_sz,'FontWeight','bold');
txt6 = text(250+2*510,75,'Crush. 4','FontSize',font_sz,'FontWeight','bold');
export_fig([dir_out 'tmp_fbirn_b.png'],sprintf('-r%d',my_dpi));
f7b = imread([dir_out 'tmp_fbirn_b.png']);

fig7 = [f7a f7b];
figure(7); clf;
imshow( fig7,'Border','tight' );
txta = text( 225, 240, 'a', 'FontSize', 1.5*font_sz, 'Color','w' );
txtb = text( 225 + size(f7a,2), 240, 'b', 'FontSize', 1.5*font_sz, 'Color','w' );
export_fig([dir_out 'Figure_birn_T1_T2.png'],sprintf('-r%d',2*my_dpi));


%% supporting information FBIRN under het. B0: T2 COV (Fig. S6)

font_sz = 28;
buff_sz = 150;

fn_8a = 'fbirn_COV.png';
fn_8b = 'fbirn_B0.png';
fn_8c = 'fbirn_B0_zoned.png';

f8a = imread( [dir_out fn_8a] );
f8b = imread( [dir_out fn_8b] );
f8b = repmat( f8b, [1 1 3] );
f8b(1:5,:,:) = 129*ones(5,size(f8b,2),3);
f8c = imread( [dir_out fn_8c] );
f8c(end-4:end,:,:) = 129*ones(5,size(f8b,2),3);

f8amod = [ 255*ones( size(f8a,1), buff_sz, 3 ) f8a ];

figure(99); clf;
imshow( f8amod,'Border','tight' );
txt1 = text(75,300,'EPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt2 = text(75,490+300,'pEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt3 = text(75,2*490+300,'ssEPG','FontSize',font_sz,'FontWeight','bold','Rotation',90);
txt4 = text(75,100+3*490+300,'ssEPG w/B0','FontSize',font_sz,'FontWeight','bold','Rotation',90);
export_fig([dir_out 'tmp_fbirn_a.png'],sprintf('-r%d',my_dpi));

f8amod = imread( [dir_out 'tmp_fbirn_a.png'] );

buffbc = size( f8amod,1 ) - (size(f8b,1) + size(f8c,1));
buffbc1 = round( buffbc/3 );
buffbc2 = buffbc1;
buffbc3 = buffbc - (buffbc1 + buffbc2);

fig8 = [ f8amod [129*ones( buffbc1, size(f8b,2),3 ); f8b; 255*ones( buffbc2, size(f8b,2),3 ); f8c; 129*ones( buffbc3, size(f8b,2),3 )] ];

figure(8); clf;
imshow(fig8,'Border','tight');
txta = text(225,75,'a','FontSize',1.5*font_sz,'Color','w');
txtb = text(75+size(f8amod,2),75,'b','FontSize',1.5*font_sz,'Color','w');
txtc = text(75+size(f8amod,2),75+size(f8b,1)+buffbc1+buffbc2,'c','FontSize',1.5*font_sz,'Color','w');
export_fig([dir_out 'Figure_fbirn_cov.png'],sprintf('-r%d',2*my_dpi));