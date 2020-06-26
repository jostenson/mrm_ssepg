function T_m = rfTransEPG_array(alpha,phi)
%rfTransEPG returns the tranistion matrix describing an RF pulse of flip 
%angle alpha (degrees) with a phase of phi (radians); phi is relative to
%the x-axis; T_m is 3x3xN matrix where N is the number of elements in alpha

N = numel( alpha );

T_m = zeros( 9, N );

alpha = alpha(:)';

T_m(1,:) = cosd(alpha/2).^2;
T_m(2,:) = exp(-2i*phi)*sind(alpha/2).^2;
T_m(3,:) = -1i/2*exp(-1i*phi)*sind(alpha);
T_m(4,:) = exp(2i*phi)*sind(alpha/2).^2;
T_m(5,:) = cosd(alpha/2).^2;
T_m(6,:) = 1i/2*exp(1i*phi)*sind(alpha);
T_m(7,:) = -1i*exp(1i*phi)*sind(alpha);
T_m(8,:) = 1i*exp(-1i*phi)*sind(alpha);
T_m(9,:) = cosd(alpha);

T_m = reshape( T_m, [3 3 numel(alpha)] );



end

