% Calculate the orbit shift in degrees longitude --> p_sat in days
function [d_lambda] = orbitShiftAN(p_sat)
  d_lambda = (-0.2507)*p_sat*(24*60); % degrees of shift per rev
end
