% Calculate transmitter-facing swath diameter for instrument
function [d] = swathTx(r_sat, r_Tx)
  in_ang = 60; % incidence angle [deg.]
  h = r_Tx - r_sat; % Instrument requires Rx lower orbit than Tx
  d = 2*(tand(in_ang)*h);
endfunction
