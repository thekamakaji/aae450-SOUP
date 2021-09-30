% Calculates swath diameter of given instrument converage at user orbit
function [d] = swathDiameter(r_sat)
  in_ang = 60; % incidence angle [deg.]
  h = r_sat - 6378.14; % alt. of sat
  d = 2*(tand(in_ang)*h);
end
