% Nodal Period Calculator
% Author: V. Swaminathan
% Purpose: To calculate the Nodal Period of an Orbit

%% Function

function P_n = nodalPeriod(P_k, J2, R_a, p, e, i)
    
    temp = sqrt(1 - e^2)*(2 - 3*(sind(i)^2)) + (4 - 5*(sind(i)^2));
    P_n = P_k/(1 + ((3*J2*R_a)/(4*p))*temp);

end