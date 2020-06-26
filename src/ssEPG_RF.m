function omega = ssEPG_RF( dtau_s,B1_T,omega0,phi,iso_shift )

% returns state matrix omega after rf pulse approximated by N hard pulse steps

% omega = 3 x Q state matrix that describres transverse and longitudinal
% magnetization in the Fourier transform domain of the through slice
% direction; steps in omega are dkz_per_length
% dtau_s = tau_s/(N-1) where tau_s is the RF pulse direction
% B1_T = the RF pulse in T
% omega0 = the starting state
% phi = the nominal RF phase in radians
% iso_shift = integer correction for change in isodelay from (N/2 + 1), see
% note below

% NOTE: the isodelay (with iso_shift = 0) should be at point N/2 + 1, N
% should be even

% NOTE: if omega is later used in EPG simulation, the distance in spatial
% frequency between any two columns in omega is coupled to dtau_s, dkz =
% dtau_s * gammabar_Hz_per_T * grad_T_per_m and so coupled to shifting due
% to crushers

gamma_rads_per_T = 2 * pi * 42.577479e6;
N = numel(B1_T);
omega = omega0;

for ii = 1:N   
    
    w1_deg = dtau_s * gamma_rads_per_T * B1_T(ii) * 180 / pi;
    
    T = rfTransEPG(-w1_deg*0.5,phi);

    omega = T*omega;
    
    omega(1:2,:) = [circshift( omega(1,:),-1,2 ); circshift( omega(2,:),1,2 )];
    omega(1,end) = 0;
    omega(2,1) = conj(omega(1,1));
    
    omega = T*omega;
    
end

N_shift = round(N/2 - iso_shift);
omega(1,:) = circshift( omega(1,:), N_shift, 2 );
omega(1,1:N_shift) = fliplr(conj(omega(2,2:N_shift+1)));
omega(2,:) = circshift( omega(2,:), -N_shift, 2 );
omega(2,(end-N_shift+1):end) = 0;

end

