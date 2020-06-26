function Ft_v = EPG_MRF_SSFP( T1,T2,TE_v,TR_v,FA_v,delk,nreps,szomega,phi_v,TI)
% returns the complex net transverse magnetization using the extended phase graph approach 

%   INPUT: T1/T2 = T1/T2 long./transverse time constants (T2' is neglected)
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


%   OUTPUT: Ft_v = vector of length nreps containting the complex net trans
%                   magnetization


% convert phi_v to radians
phi_v = phi_v.*pi/180;

% initialize magnetization w/inversion
omega = [[0,0,-1]' zeros(3,szomega-1)];

% inversion delay
e1 = exp(-TI/T1);
E = diag([0 0 e1]);
omega(:,1) = E*omega(:,1) + [0 0 (1-e1)]';

% iterate for defined reps
Ft_v = zeros(1,nreps);
for j=1:nreps
    
    alpha = FA_v(j);
    TR = TR_v(j);
    phi = phi_v(j);
    TE = TE_v(j);
    
    % flip angle transition
    T_m = rfTransEPG(alpha,phi);
    omega = T_m*omega;
    
    % decay for readout at TE
    e2 = exp(-TE/T2);
    e1 = exp(-TE/T1);
    E = diag([e2 e2 e1]);
    omega = E*omega;
    omega(3,1) = omega(3,1) + (1-e1);
    
    % store readout transverse magnetization
    Ft_v(j) = omega(1,1);
    
    % apply gradient dephasing specified by delk
    omega(1,:) = circshift(omega(1,:),delk,2);
    omega(2,:) = circshift(omega(2,:),-delk,2);
    omega(2,(end-delk+1):end) = 0;
    omega(1,1) = conj(omega(2,1));
    
    % decay for flip angle at TR
    e2 = exp(-(TR-TE)/T2);
    e1 = exp(-(TR-TE)/T1);
    E = diag([e2 e2 e1]);
    omega = E*omega;
    omega(3,1) = omega(3,1) + (1-e1);
    
end

end

