function output_match = MRF_dict_match_B0( input_recon, input_dict )
% do MRF dictionary construction/match
% used in conjunction with MRF processing

%   INPUT: input_recon.MRF_img_stack_coil_combined = Nr x Nc x nSig matrix of MRF img data
%                 ".B1_map = provided B1 map
%          input_dict.dict_compress = SVD compressed dictionary
%                   ".U_r = significant left singular vectors
%                   ".dict_list = table of dictionary parameters
%                   ".dict_norm = normalized full dictionary
%                   ".B1_compress_dict_list_v = list of B1 for cols of
%                       compressed dictionary
%          reduce_flag = 1 then dictionary is reduced by SVD, else no
%   OUTPUT: output.T1_map = T1 map
%                ".T2_map = T2 map
%                ".B0_map = B0 map
%                ".B1_map = B1 map
%                ".M0_map = magnetization (proton) density
%                ".R_map = complex correlation map
%                ".dict_list = list of dictionary parameters that
%                   correspond to the entries in the dictionary
%                ".dict_fn = input dictionary filename

disp('Doing MRF dictionary match...');
tic;

msr = input_recon.MRF_img_stack_coil_combined;
[ n_r,n_c,n_reps ] = size( msr );
n_msr = n_r * n_c;
msr = reshape( permute( msr, [3 1 2] ), [n_reps n_msr] );

B1_v = input_dict.B1_v;
n_B1 = numel( B1_v );
B1_map = input_recon.B1_map;
dict_compress = input_dict.dict_compress;
dict_list = input_dict.dict_list;
U_r = input_dict.U_r;
B1_compress_dict_list_v = input_dict.B1_compress_dict_list_v;

% get discretized B1
B1_map = B1_map(:);
B1_diff = abs( repmat( B1_map, [1 n_B1] ) - B1_v(:)' );
[~,idx_B1_min] = min( B1_diff,[], 2 );
B1_map_disc = B1_v( idx_B1_min );

% match
T1_map = zeros( n_msr, 1 );
T2_map = T1_map;
B0_map = T1_map;
M0_map = T1_map;
for ii = 1:n_B1
    
    fprintf('MRF_dict_match: doing B1 %d of %d ...\n',ii,n_B1 )
    
    my_B1 = B1_v(ii);
    my_idx_msr = my_B1 == B1_map_disc;
    my_idx_u = my_B1 == B1_compress_dict_list_v;
    my_idx_dict = my_B1 == dict_list(:,4);
    
    sub_msr = msr( :,my_idx_msr );
    my_U_r = U_r( :,my_idx_u );
    sub_dict = dict_compress( :,my_idx_dict );
    sub_list = dict_list( my_idx_dict,: );

    sub_msr_proj = my_U_r' * sub_msr;
%     IP = sub_msr_proj' * sub_dict;
%     [ my_M0, my_idx_max ] = max( IP,[],2 );
    [ my_M0, my_idx_max ] = ip_match( sub_msr_proj, sub_dict, 10 );
    
    T1_map( my_idx_msr ) = sub_list( my_idx_max,1 );
    T2_map( my_idx_msr ) = sub_list( my_idx_max,2 );
    B0_map( my_idx_msr ) = sub_list( my_idx_max,3 );
    M0_map( my_idx_msr ) = my_M0;
    
end
T1_map = reshape( T1_map, [n_r n_c] );
T2_map = reshape( T2_map, [n_r n_c] );
B0_map = reshape( B0_map, [n_r n_c] );
M0_map = reshape( M0_map, [n_r n_c] );
B1_map_disc = reshape( B1_map_disc, [n_r n_c] );

output_match.T1_map = T1_map;
output_match.T2_map = T2_map;
output_match.B0_map = B0_map;
output_match.B1_map = B1_map_disc;
output_match.M0_map = M0_map;
output_match.dict_list = dict_list;
% output_match.dict_fn = input_dict.fn;

t = toc;
disp(['Doing dictionary match complete. Elapsed time is ' num2str(t) ' s.']);

function [ maxi, idx_max ] = ip_match( A, B, n_blk )

AH = A';
n_row = size( AH,1 );
n_col_B = size( B,2 );
n_sub_B = ceil( n_col_B/n_blk );

% loop over blocks
cur_max = zeros( n_row,2 );
cur_idx_max = cur_max;
for ii = 1:n_blk

    idx_B = ( (ii-1)*n_sub_B + 1 ):( ii*n_sub_B );
    idx_B( idx_B > n_col_B ) = [];
    sub_B = B(:,idx_B);
    
    % get inner product and max
    sub_IP = AH * sub_B;
    [cur_max(:,2),my_idx_max] = max( sub_IP,[],2 );
    cur_idx_max(:,2) = my_idx_max + (ii-1)*n_sub_B;
    [cur_max(:,1),idx_compare] = max( cur_max,[],2 );
   
    cur_idx_max( idx_compare == 1,1 ) = cur_idx_max( idx_compare == 1,1 );
    cur_idx_max( idx_compare == 2,1 ) = cur_idx_max( idx_compare == 2,2 );

end

maxi = cur_max(:,1);
idx_max = cur_idx_max(:,1);