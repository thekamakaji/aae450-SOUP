% SSO Calculator
% Author: V. Swaminathan
% Purpose: Calculates SSO inclination given orbital altitude
%          for circular orbits with a < 12352 km

%% Function
function i_SSO = makeSSO(r_a)
    
    i_SSO = acosd(-1* ((r_a/12352)^(7/2)));

end