% Surface Longitude Calculator
% Author: V. Swaminathan
% Purpose: To calculate the Surface Longitude of a pass at target latitude

%% Function

function surfLong = surfaceLong(w, v, AN, i, dLong)

    num1 = cosd(w + v)*sind(AN) + sind(w + v)*cosd(AN)*cosd(i);
    den1 = cosd(w + v)*cosd(AN) - sind(w + v)*sind(AN)*cosd(i);
    
    surfLong = atand(num1/den1) + (v/(2*pi))*dLong;

end