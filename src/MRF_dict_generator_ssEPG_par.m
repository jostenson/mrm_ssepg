function output = MRF_dict_generator_ssEPG_par(input)
% make an MRF dictionary using the slice-selective EPG (ssEPG) method using
% parallelization


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
B1_v = input.B1_v;
TR_v = input.TR_v;
TE_v = input.TE_v;
FA_v = input.FA_v;
phi_v = input.phi_v;
rf_shape = input.rf_shape;

nT2 = numel(T2_v);
nT1 = numel(T1_v);
nB1 = numel(B1_v);

%% dictionary init

dict_list = zeros(nT2*nT1*nB1,3);

% determine dictionary list
nn = 1;
for kk = 1:nB1
    
    B1 = B1_v(kk);
    
    for jj = 1:nT2
        
        T2 = T2_v(jj);
        
        for ii = 1:nT1
            
            T1 = T1_v(ii);
            if T1 >= T2
                
                dict_list(nn,:) = [T1, T2, B1];
                nn = nn+1;
                
            end
            
        end
    end
end

% remove zero rows from listT1T2_m, normFt_m, Ft_m since init matrices are
% too big
dict_list(~any(dict_list,2),:) = [];

dict_length = size( dict_list,1 );

dict_norm = zeros(nreps,dict_length);

%% dictionary generation


B1_unique_v = unique( dict_list(:,3) );
nT1T2 = sum( B1_unique_v(1) == dict_list(:,3) );


for ii = 1:nB1
    
    my_idx_B1 = B1_unique_v(ii) == dict_list(:,3);
    my_T1_list = dict_list(my_idx_B1,1);
    my_T2_list = dict_list(my_idx_B1,2);
    my_B1 = B1_unique_v(ii);
    disp(['MRF dictionary construction: B1 ', num2str(my_B1)]);
    
    myFA_v = FA_v.*my_B1;
    
    % run sequence using EPG
    [~,sigs] = EPG_MRF_SSFP_prof_par( my_T1_list, my_T2_list, TE_v, TR_v, myFA_v, delk, nreps, szomega, phi_v, TI, rf_shape, nT1T2, kerns, 0);
    
    dict_norm(:,my_idx_B1) = sigs.*repmat( 1./vecnorm(sigs), [nreps 1] );   
        
end


% create output structure
output = input;

output.dict_list = dict_list;
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
        my_dict_norm = dict_norm(:, (dict_list(:,3) == B1) );
        
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