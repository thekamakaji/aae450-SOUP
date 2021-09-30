% Orbit Determination for SOoP Satellite
% V. Swaminathan - Team SOUP

% Universal Constants
R_E = 6378.14; % Radius of Earth [km]
C_E = 40002; % Circumference of Earth [km]
mu_E = 3.986e5; % Grav. Parameter of Earth [km^3 / s^2]
J2_E = 0.00108; % J2 parameter of Earth

% Keplerian Orbit of Transmitter (Tx) Satellite - HARDCODE
orb_Tx_a = 20180 + R_E; % Orbit radius [km]
orb_Tx_e = 1; % Eccentricity
orb_Tx_i = 86.4; % Inclination [deg.]
orb_Tx_RAAN = 0; % Right Ascension of AN [deg.]
orb_Tx_P = 1; % Orbital Period [days]

% Keplerian Orbit of Receiver (Rx) Satellite - HARDCODE
orb_Rx_a = 750 + R_E; % Orbit radius [km]
orb_Rx_e = 1; % Eccentricity
orb_Rx_i = 86.4; % Inclination [deg.]
orb_Rx_RAAN = 0; % Right Ascension of AN [deg.]

% Keplerian Orbit of Satellites - Calculated Values (Update vars.)
orb_Tx_P = periodFromRadius(orb_Tx_a, mu_E); % Orbital Period [days]
orb_Rx_P = periodFromRadius(orb_Rx_a, mu_E); % Orbital Period [days]

g_shiftRx = orbitShiftAN(orb_Rx_P) % Shift westward per orbit [deg.] --> works
g_shiftRx_km = g_shiftRx*R_E*(pi/180) % Shift westward per orbit [km]
g_swath_km = swathDiameter(orb_Rx_a) % Swath diameter of Rx instrument [km]
RxTx_swath_km = swathTx(orb_Rx_a, orb_Tx_a) % Swath diameter of Rx facing Tx [km]
spec_points_SSM = specularNum(g_swath_km, 100) % Specular points needed for Soil Surface Moisture (10km x 10km resolution)
spec_points_RZSM = specularNum(g_swath_km, 1600) % Specular points needed for Root Zone Soil Moisture (40km x 40km resolution)
numSatsNeeded_SC1 = C_E / (2*g_swath_km) % Number of Satellites needed to cover all of earth, every orbit (High value)
numSatsNeeded_SC2 = C_E / (2*abs(g_shiftRx_km)) % Number of Satellites needed to cover all of earth (Med value)

% Specific Calculations
d_bw_satsInPlane = ((orb_Tx_a + R_E)*2*pi) / 11 % Distance b/w IridiumNEXT satellites within each orbital plane [km]
target_RxTx_swath = d_bw_satsInPlane / 2 % Target diameter of RxTx swath to ensure continuous signal reception [km]
paramRxTx_swath = target_RxTx_swath / RxTx_swath_km % Get to < 1 in order to have continuous coverage

% Synodic Period Calculation (2-day revisit time minimum)
if orb_Rx_P < orb_Tx_P
  P_syn = 1 / ((1/orb_Rx_P)-(1/orb_Tx_P)) % Synodic Period [days]
end
if orb_Rx_P > orb_Tx_P
  P_syn = 1 / ((1/orb_Tx_P)-(1/orb_Rx_P)) % Synodic Period [days]
end
