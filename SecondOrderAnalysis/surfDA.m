% Surface Dihedral Angle
% Author: V. Swaminathan
% Purpose: To calculate the surface dihedral angle (analogous to the
%          coverage in longitude of sensor at target latitude) using
%          spherical trigonometry

%% Function

function Aphi = surfDA(HGRA, Lat)

    Aphi = 2*acosd((cosd(HGRA)-(sind(Lat)^2))/(cosd(Lat)^2)); % Surf. Dihedral Angle [deg.]

end