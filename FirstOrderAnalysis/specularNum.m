% Calculates number of specular points needed for given resolution and area
function [n] = specularNum(d_swath, resolution)
  A = pi*(d_swath^2)/4; % Area covered by ground swath
  n = A / resolution; % Specular grid points needed
end
