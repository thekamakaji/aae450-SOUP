% Half Ground Range Angle Calculator
% Author: V. Swaminathan
% Purpose: To calculate ground range angle seen by given sensor and orbit
%          info for satellite.

%% Function

function theta = halfGRA(rs, phi, R_geo)

    gamma = asind(rs*sind(phi)/R_geo); % Intermediate Angle
    rho = R_geo*cosd(gamma) + rs*cosd(phi); % Slant Range to edge of sensor coverage area
    
    theta = asind(rho*sind(phi)/R_geo); % HGRA
    
end