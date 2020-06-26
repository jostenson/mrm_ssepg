function [Fxy,Ft] = EPG_MRF_SSFP_prof_par( T1_list,T2_list,TE_v,TR_v,FA_v,delk,nreps,szomega,phi_v,TI,rf_shape,n_T1T2,kerns,ksp_flag)
% returns the complex net transverse magnetization using the extended phase graph approach
% for use use in the slice-selective EPG (ssEPG) method, parallelized to compute
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
%          phi_v = vector of RF phase angles in degress
%          TI = inversion time
%          n_T1T2 = number of elements in the subset (number of elements in
%          T1_list)
%          rf_shape = structure containing RF_shape, see also
%               get_rf_shape.m
%          kerns = GPU kernel for 
%          ksp_flag = if 1 then compute and return k-space of slice profile

%   OUTPUT: Ft_v = vector of length nreps containting the complex net trans
%                   magnetization
%           Fxy = array nreps x 2*szomega - 2 of the through=slice
%                   magnetization in k-space

k = kerns.k2;
k.ThreadBlockSize = [128 1 1];
k.GridSize = [64 1 1];

diso_idx = rf_shape.diso_idx;
dtau_s = rf_shape.dtau_s;
B1_T = rf_shape.B1_T;
FA_nom_deg = rf_shape.FA_nom_deg;

% convert phi_v to radians
phi_v = phi_v.*pi/180;

% initialize magnetization w/inversion
omega = zeros( 3,szomega,n_T1T2 );
omega(3,1,:) = -1;

% inversion delay
omega = relax( omega, TI, T1_list, T2_list );

if ksp_flag == 1
    Fxy = zeros( nreps,2*(szomega - 1),n_T1T2 );
end
Ft = zeros(nreps,n_T1T2);
for ii=1:nreps
        
    alpha = FA_v(ii);
    TR = TR_v(ii);
    phi = phi_v(ii);
    TE = TE_v(ii);
    
    % flip angle transition
    omega = ssEPG_RF_par( dtau_s,alpha/FA_nom_deg*B1_T,omega,phi,diso_idx,k );
    
    % decay for readout at TE
    omega = relax( omega, TE, T1_list, T2_list );
    
    % debug
%     figure(99); clf;
%     plot(log10(abs(omega(1,:))))
%     hold on
%     plot(log10(abs(omega(2,:))))
%     hold off
%     drawnow
        
    % store readout transverse magnetization
    Ft(ii,:) = omega(1,1,:);
    if ksp_flag == 1
        Fxy(ii,:,:) = [ fliplr( conj(omega(2,2:end,:)) ), omega(1,1:end-1,:) ];
    end    
    
    % apply gradient dephasing specified by delk
    if delk >= 0
        omega(1,:,:) = circshift(omega(1,:,:),delk,2);
        omega(1,1:delk,:) = fliplr(conj(omega(2,2:delk+1,:)));
        omega(2,:,:) = circshift(omega(2,1:end,:),-delk,2);
        omega(2,(end-delk+1):end,:) = 0;
    else
        omega(2,:,:) = circshift(omega(2,:,:),-delk,2);
        omega(2,1:-delk,:) = fliplr(conj(omega(1,2:-delk+1,:)));
        omega(1,:,:) = circshift(omega(1,1:end,:),delk,2);
        omega(1,(end+delk+1):end,:) = 0;
    end
    
    % decay for flip angle at TR
    tau = TR-TE;
    omega = relax( omega, tau, T1_list, T2_list );
    
    
end

if ksp_flag ~= 1
    Fxy = [];
end

end

function w = relax( w,tau,T1_v,T2_v )


e2_v = exp( -tau./T2_v(:).' );
e1_v = exp( -tau./T1_v(:).' );

w(1,:,:) = squeeze( w(1,:,:) ).* e2_v;
w(2,:,:) = squeeze( w(2,:,:) ).* e2_v;
w(3,:,:) = squeeze( w(3,:,:) ).* e1_v;

w(3,1,:) = squeeze( w(3,1,:) ) + (1-e1_v(:));

end