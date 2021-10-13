% AN Nodal Regression Rate Calculator
% Author: V. Swaminathan
% Purpose: To calculate the Nodal Regression of an Orbit's AN taking into
%          account J2, J4, J6 zonal harmonics

%% Function

function d_AN = regressionRateAN(n, J2, R_a, p, i, e)
    
    term_J2 = -(3/2)*n*J2*((R_a/p)^2)*cosd(i);
    
    temp = 12 - 4*(e^2) - (80+5*(e^2))*(sind(i)^2);
    term_J4J6 = (3/32)*n*(J2^2)*((R_a/p)^4)*cosd(i)*temp;
    
    d_AN = term_J2 + term_J4J6;

end