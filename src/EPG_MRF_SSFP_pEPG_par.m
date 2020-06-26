function Ft = EPG_MRF_SSFP_pEPG_par( T1_list,T2_list,TE_v,TR_v,FA_v,delk,nreps,szomega,phi_v,TI,n_T1T2)
% returns the complex net transverse magnetization using the extended phase graph approach
% for use use in the parcellated/partitioned EPG (pEPG) method, parallelized to compute
% a large subset of atoms simultaneously


%   INPUT: T1_list/T2_list = lists of T1/T2 long./transverse time
%           constants, the lists should repeat so that [T1_list(:)
%           T2_list(:)] list every combination of T1 and T2
%          TE_v = vector of echo times
%          TR_v = vector of repetition times
%          FA_v = vector of flip angles in degrees
%          delk = dephasing step (shoud be 1)
%          nreps = number of excitations
%          szomega = pos integer length of the state vectors used to
%               calculate the magnetization by the EPG formalism (too short
%               a vector will change the answer since some
%               magnetization with be lost by a truncated state vector)
%          phi_v = vecotr of RF phase angles in degress
%          TI = inversion time
%          n_T1T2 = number of elements in the subset (number of elements in
%          T1_list)


%   OUTPUT: Ft_v = vector of length nreps containting the complex net trans
%                   magnetization

d_T1_v = gpuArray( T1_list );
d_T2_v = gpuArray( T2_list );

o1 = 3; o2 = szomega; o3 = n_T1T2;
o123 = o1*o2*o3; o12 = o1*o2; o23 = o2*o3;

% lookups
idx0 = reshape( 1:o123, [o1 o2 o3] );
idx = idx0;
idx(1,:,:) = circshift( idx0(1,:,:), 1, 2 );
idx(2,:,:) = circshift( idx0(2,:,:), -1, 2 );
idx(3,:,:) = idx0(3,:,:);
d_idx = gpuArray( idx(:) );

d_idx_conj_to = gpuArray( 1:o12:o123 );
d_idx_conj_from = gpuArray( 2:o12:o123 );
d_idx_zero = gpuArray( ( o12 - 1 ):o12:o123 );

% convert phi_v to radians
phi_v = phi_v.*pi/180;

% initialize magnetization w/inversion
d_omega = gpuArray( single( complex( zeros( o1,o2,o3 ) ) ) );
d_omega(3,1,:) = -1;

% inversion delay
d_omega = relax( d_omega, TI, d_T1_v, d_T2_v );

d_Ft = gpuArray( zeros(nreps,n_T1T2) );
for ii=1:nreps
    
    alpha = FA_v(ii);
    TR = TR_v(ii);
    phi = phi_v(ii);
    TE = TE_v(ii);
    
    % flip angle transition
    d_T = gpuArray( rfTransEPG(alpha,phi) );
    d_omega = reshape( d_T * reshape( d_omega, [o1 o23] ), [o1 o2 o3] );
    
    % decay for readout at TE
    d_omega = relax( d_omega, TE, d_T1_v, d_T2_v );
    
    % debug
    %     figure(99); clf;
    %     plot(log10(abs(omega(1,:))))
    %     hold on
    %     plot(log10(abs(omega(2,:))))
    %     hold off
    %     drawnow
    
    % store readout transverse magnetization
    d_Ft(ii,:) = d_omega(1,1,:);
    
    % apply gradient dephasing specified by delk
    d_tmp = d_omega( d_idx );
    d_omega = reshape( d_tmp, [o1 o2 o3] );
    d_omega( d_idx_conj_to ) = conj( d_omega( d_idx_conj_from ) );
    d_omega( d_idx_zero ) = 0;
    
    % decay for flip angle at TR
    tau = TR-TE;
    d_omega = relax( d_omega, tau, d_T1_v, d_T2_v );
    
end

Ft = gather( d_Ft );

end

function w = relax( w,tau,T1_v,T2_v )


e2_v = exp( -tau./T2_v(:).' );
e1_v = exp( -tau./T1_v(:).' );

w(1,:,:) = squeeze( w(1,:,:) ).* e2_v;
w(2,:,:) = squeeze( w(2,:,:) ).* e2_v;
w(3,:,:) = squeeze( w(3,:,:) ).* e1_v;

w(3,1,:) = squeeze( w(3,1,:) ) + (1-e1_v(:));

end
