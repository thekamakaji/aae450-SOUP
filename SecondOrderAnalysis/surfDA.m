% Surface Dihedral Angle
% Author: V. Swaminathan
% Purpose: To calculate the surface dihedral angle (analogous to the
%          coverage in longitude of sensor at target latitude) using
%          spherical trigonometry

%% Function

function Aphi = surfDA(HGRA, Lat)

    Aphi = 2*acosd((cos(HGRA)-(sin(Lat)^2))/(cos(Lat)^2)); % Surf. Dihedral Angle [deg.]

end