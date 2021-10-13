% Orbital Period Calculator
% Author: V. Swaminathan
% Purpose: To calculate the Keplerian Orbital Period

%% Function

function P_k = keplerianPeriod(a, mu)

    P_k = 2*pi*sqrt((a^3)/mu);

end