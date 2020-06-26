function dM = bloch_eq(t,M,rf,f0,t_rf)
% function dM = bloch_eq(t,M,b1,f0)
% See Eqn 3.4 of Bernstein


gamma_rads_per_T = 2*pi*42.577479e6;
delta_rads = f0*2*pi;

b1 = interp1(t_rf,rf,t);

dM(1) = delta_rads * M(2);
dM(2) = gamma_rads_per_T * b1 * M(3) - delta_rads * M(1);
dM(3) = -gamma_rads_per_T * b1 * M(2);

dM = dM';