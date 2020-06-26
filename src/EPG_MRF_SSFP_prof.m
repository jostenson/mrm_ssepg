function [Fxy,Ft_v] = EPG_MRF_SSFP_prof( T1,T2,TE_v,TR_v,FA_v,delk,nreps,szomega,phi_v,TI,rf_shape)
% returns the complex net transverse magnetization and through-slice k-space representation
% of the magnetization using the extended phase graph approach based on a pre-calculated soft RF pulse

%   INPUT: T1/T2 = T1/T2 long./transverse time constants (T2' is neglected)
%          TE_v = vector of echo times
%          TR_v = vector of repetition times
%          FA_v = vector of flip angles in degrees
%          delk = specified in rf_shape is the number of dk (states)
%               shifted for one crusher
%          nreps = number of excitations
%          szomega = pos integer length of the state vectors used to
%               calculate the magnetization by the EPG formalism (too short
%               a vector will change the answer since some
%               magnetization with be lost by a truncated state vector)
%          phi_v = vector of RF phase angles in degress
%          TI = inversion time
%          rf_shape = structure containing RF_shape, see also
%               get_rf_shape.m

%   OUTPUT: Ft_v = vector of length nreps containting the complex net trans
%                   magnetization
%           Fxy = array nreps x 2*szomega - 2 of the through=slice
%                   magnetization in k-space

diso_idx = rf_shape.diso_idx;
dtau_s = rf_shape.dtau_s;
B1_T = rf_shape.B1_T;
FA_nom_deg = rf_shape.FA_nom_deg;

% convert phi_v to radians
phi_v = phi_v.*pi/180;

% initialize magnetization w/inversion
omega = [[0,0,-1]' zeros(3,szomega-1)];

% inversion delay
e1 = exp(-TI/T1);
E = diag([0 0 e1]);
omega(:,1) = E*omega(:,1) + [0 0 (1-e1)]';

Fxy = zeros( nreps, 2*(szomega - 1) );
Ft_v = zeros( 1,nreps );
for ii=1:nreps
    
    alpha = FA_v(ii);
    TR = TR_v(ii);
    phi = phi_v(ii);
    TE = TE_v(ii);
    
    % flip angle transition
    omega = ssEPG_RF( dtau_s,alpha/FA_nom_deg*B1_T,omega,phi,diso_idx );
    
    % decay for readout at TE
    e2 = exp(-TE/T2);
    e1 = exp(-TE/T1);
    E = diag([e2 e2 e1]);
    omega = E*omega;
    omega(:,1) = omega(:,1) + [0 0 (1-e1)]';
    
    % debug
%     figure(99); clf;
%     plot(log10(abs(omega(1,:))))
%     hold on
%     plot(log10(abs(omega(2,:))))
%     hold off
%     drawnow
        
    % store readout transverse magnetization
    Ft_v(ii) = omega(1,1);
    Fxy(ii,:) = [fliplr(conj(omega(2,2:end))) omega(1,1:end-1) ];
    
    % apply gradient dephasing specified by delk
    if delk >= 0
        omega(1,:) = circshift(omega(1,:),delk,2);
        omega(1,1:delk) = fliplr(conj(omega(2,2:delk+1)));
        omega(2,:) = circshift(omega(2,:),-delk,2);
        omega(2,(end-delk+1):end) = 0;
    else
        omega(2,:) = circshift(omega(2,:),-delk,2);
        omega(2,1:-delk) = fliplr(conj(omega(1,2:-delk+1)));
        omega(1,:) = circshift(omega(1,:),delk,2);
        omega(1,(end+delk+1):end) = 0;
    end
    
    % decay for flip angle at TR
    e2 = exp(-(TR-TE)/T2);
    e1 = exp(-(TR-TE)/T1);
    E = diag([e2 e2 e1]);
    omega = E*omega;
    omega(:,1) = omega(:,1) + [0 0 (1-e1)]';
    
    
end

end

