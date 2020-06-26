function output = MRF_dict_generator_ssEPG_B0_par(input)
% make an MRF dictionary using the slice-selective EPG (ssEPG) method using
% parallelization, ASSUME FIXED TR FOR B0 AND SYMMETRIC RF PULSE to improve
% processing speed 


%   INPUT: input.nreps = number of frames (or TRs) in MRF sequence
%              ".TI = inversion time
%              ".delk = pos integer is step between states equal to a full
%                   dephasing imparted by crusher gradient
%              ".szomega = pos integer is number of factors of k to 
%                   include in phase history
%              ".T1_v = vector of T1s in dictionary
%              ".T2_v = vector of T2s in dictionary
%              ".B1_v = vector of B1 multipliers in dictionary
%              ".TR_v = vector of repetition times
%              ".TE_v = vector of echo times
%              ".FA_v = vector of flip angles in degrees
%              ".phi_v = vector of RF phases
%              ".reduce = 1 then do SVD compression of dictionary
%              ".sfrac = fraction of dictionary energy to retain if
%                   compressing
%              ".rf_shape = soft pulse RF information, see also
%                   get_rf_shape.m

%   OUTPUT: output.dict_list = matrix describing T1, T2, B1, etc. for all 
%                       realisitic combinations
%                ".dict_norm = normalized dictionary w/columns
%                       parameterized by dict_list
%                ".dict = dictionary w/columns parameterized by dict_list
%                ".V_red = reduced right singular vectors of dictionary
%                ".dict_red = compressed dictionary
%                ". = all input params copied to output params
%                ".B1_compress_dict_list_v = vector of B1s describing column
%                       B1s of compressed dictionary and associate matrices


disp('Constructing MRF dictionary...')
tic
%% get gpu kernels

k1 = parallel.gpu.CUDAKernel('ssepg_functions.ptx','ssepg_functions.cu','dephase_gradients_rf_stage1');
k2 = parallel.gpu.CUDAKernel('ssepg_functions.ptx','ssepg_functions.cu','dephase_gradients_rf_stage2');
kerns.k1 = k1; kerns.k2 = k2;

%% declare parameters
nreps = input.nreps;
TI = input.TI;
delk = input.delk; 
szomega = input.szomega;
T1_v = input.T1_v;
T2_v = input.T2_v;
B0_v = input.B0_v;
B1_v = input.B1_v;
TR_v = input.TR_v;
TE_v = input.TE_v;
FA_v = input.FA_v;
phi_v = input.phi_v;
rf_shape = input.rf_shape;
i_sign = input.i_sign;

nT2 = numel(T2_v);
nT1 = numel(T1_v);
nB0 = numel(B0_v);
nB1 = numel(B1_v);

%%

if any( TR_v(1) ~= TR_v )
    error('TRs should all be equal to allow speedup of B0 entries.')
end

%% dictionary init

dict_list = zeros(nT2*nT1*nB0*nB1,4);

% determine dictionary list
nn = 1;
for kk = 1:nB1
    
    B1 = B1_v(kk);
    
    for mm = 1:nB0
        
        B0 = B0_v(mm);
        
        for jj = 1:nT2
            
            T2 = T2_v(jj);
            
            for ii = 1:nT1
                
                T1 = T1_v(ii);
                if T1 >= T2
                    
                    dict_list(nn,:) = [T1, T2, B0, B1];
                    nn = nn+1;
                    
                end
                
            end
        end
    end
end

% remove zero rows from dict_list etc. since initial matrices are
% too big
dict_list(~any(dict_list,2),:) = [];

dict_length = size( dict_list,1 );

dict_norm = complex( single( zeros(nreps,dict_length) ) );

%% dictionary generation

dict_list_temp = dict_list; 
dict_list_temp(:,1:3) = 0;

B1_unique_v = unique( dict_list(:,4) );

B0_sub_min = 0;
d_B0 = B0_v(2) - B0_v(1);
B0_sub_max = 1/TR_v(1) - d_B0;
n_B0_rep = max( abs(B0_v) )/(B0_sub_max + d_B0); % number of 1/TR in half B0 range

for ii = 1:nB1
    
    my_idx_B1 = abs( B1_unique_v(ii) - dict_list(:,4)) < eps;
    my_T1_list = dict_list(my_idx_B1,1);
    my_T2_list = dict_list(my_idx_B1,2);
    my_B0_list = dict_list(my_idx_B1,3);
    
    my_idx_subB0 = my_B0_list >= B0_sub_min-eps & my_B0_list <= B0_sub_max+eps;
    
    my_B1 = B1_unique_v(ii);
    disp(['MRF dictionary construction: B1 ', num2str(my_B1)]);
    
    myFA_v = FA_v.*my_B1;
    
    % run sequence using EPG
    n_set = sum( my_idx_subB0 );
    my_T1_sub_list = my_T1_list( my_idx_subB0 );
    my_T2_sub_list = my_T2_list( my_idx_subB0 );
    my_B0_sub_list = my_B0_list( my_idx_subB0 );
    [~,sigs] = EPG_MRF_SSFP_B0_prof_par( my_T1_sub_list, my_T2_sub_list, my_B0_sub_list, TE_v, TR_v, myFA_v, delk, nreps, szomega, phi_v, TI, rf_shape, n_set, kerns, 0);
    
    % determine other B0 entries based on calculated basis
    sigs_full = zeros( nreps, 2*n_B0_rep*n_set );
    sub_list = zeros(2*n_B0_rep*n_set,3);
    for jj = 1:n_B0_rep
        
        my_idx_sigs_full = ( (jj-1)*2*n_set + 1 ):( jj*2*n_set );
        my_B0s_p = (jj-1)/TR_v(1) * ones( size(my_B0_sub_list(:)') );
        my_B0s_n = -jj/TR_v(1) * ones( size(my_B0_sub_list(:)') );
        
        E_modp = exp( i_sign * 1i * 2 * pi * TE_v(:) * my_B0s_p(:)' );
        E_modn = exp( i_sign * 1i * 2 * pi * TE_v(:) * my_B0s_n(:)' );
        sigs_full( :,my_idx_sigs_full ) = [ (sigs.* E_modp) (sigs.* E_modn) ];
        
        % modify dict_list_temp
        sub_list(my_idx_sigs_full,:) = [[my_T1_sub_list(:);my_T1_sub_list(:)] [my_T2_sub_list(:);my_T2_sub_list(:)] [my_B0s_p(:)+my_B0_sub_list(:);my_B0s_n(:)+my_B0_sub_list(:)]];
 
    end
           
    tmp = 1./vecnorm(sigs_full);
    dict_norm(:,my_idx_B1) = sigs_full.*repmat( tmp, [nreps 1] );
    dict_list_temp(my_idx_B1,1:3) = sub_list;

end


% create output structure
output = input;

% output.dict_list = dict_list;
output.dict_list = dict_list_temp;
output.dict_norm = dict_norm;

%% determine reduced dictionary space

% compress dictionary via SVD method by McGivney et al,
% IEEE MI, 2014

if input.reduce == 1
    disp('Compressing MRF dictionary by SVD...')
    
    output.U_r = [];
    output.dict_compress = [];
    output.B1_compress_dict_list_v = [];
    for ii = 1:numel(B1_v); % compress dictionary by B1 discretizations
        
        B1 = B1_v(ii);
        my_dict_norm = dict_norm(:, (dict_list(:,4) == B1) );
        
        [U,S,~] = svd(my_dict_norm,'econ');
        
        if ii == 1
            s_v = diag(S);
            fNRG_v = cumsum(s_v.^2)./sum(s_v.^2);
            nDictSpace = sum(fNRG_v <= input.sfrac);
        end
        
        U_r = U(:,1:nDictSpace);
        my_dict_compress = U_r'*my_dict_norm;
        
        output.U_r = [output.U_r U_r];
        output.dict_compress = [output.dict_compress my_dict_compress];
        output.B1_compress_dict_list_v = [output.B1_compress_dict_list_v B1*ones(1,nDictSpace)];
    
    end
    
    disp('Compression of MRF dictionary by SVD complete.')
end


%%

t = toc;
output.gen_time = t;
disp(['Construction of MRF dictionary complete. Elapsed time is ' num2str(t) ' s.'])

end