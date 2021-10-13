% Geodetic Radius Calculator
% Author: V. Swaminathan
% Purpose: To calculate the geodetic radius of Earth at a target latitude

%% Function

function RLat = geodeticR(R_A, R_B, Lat)
    
    num = (((R_A^2)*cosd(Lat))^2) + (((R_B^2)*sind(Lat))^2);
    den = ((R_A*cosd(Lat))^2) + ((R_B*sind(Lat))^2);
    
    RLat = sqrt(num/den);

end