% Calculates period of orbit [days] from orbital radius [km]
function [p] = periodFromRadius(r_sat, mu)
  p = (2*pi*sqrt(r_sat^3 / mu)) / (3600*24);
endfunction
