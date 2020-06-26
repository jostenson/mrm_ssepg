function [sig_norm,sig_full] = pEPG_SSFP(input)

T1 = input.T1;
T2 = input.T2;
TE_v = input.TE_v;
TR_v = input.TR_v;
FA_v = input.FA_v;
delk = input.delk;
n_reps = input.n_reps;
szomega = input.szomega;
phi_v = input.phi_v;
TI = input.TI;
N_p = input.N_p;
fa_factors = input.fa_factors;


sig = zeros( n_reps, 1 );
sig_full = zeros( n_reps, N_p );
for ii = 1:N_p
    my_factor = fa_factors(ii);
    my_sig = EPG_MRF_SSFP( T1, T2, TE_v, TR_v, my_factor*FA_v, delk, n_reps, szomega, phi_v, TI);
    sig = my_sig(:) + sig;
    sig_full(:,ii) = my_sig;
end

sig_norm = sig./norm(sig);

fprintf('pEPG calc. complete\n')

end

