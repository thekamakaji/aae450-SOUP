% Propulsion Sizing Code
% Author: V. Swaminathan
% Version: 11/18/2021 1805 EST
% Purpose: To help select Team SoUP's satellite propulsion system and ideal
%          propellant mass
% 
% To Be Provided By Code:
% - Propellant cost of any given orbital changes (% of mass)
% - Total Expected DeltaV requirement
%

%% Initialize Workspace

clear;
close all;

%% Constant Parameters
R_A = 6378.137; % Earth Equatorial Radius [km]
R_B = 6356.752; % Earth Polar Radius [km]
MU = 3.986e5; % Gravitational Param. for Earth [km^3/s^2]
J2 = 1082.63e-6; % J2 param. for Earth
SideDay_E = 86164; % Sidereal Day Length [s]
SolarDay_E = 86400; % Solar Day Length [s]
W_E = (-0.2507)/60; % Rate of Earth Rotation [deg./s]
g = 9.81; % Gravitational Acceleration Constant on Earth [m/s^2]

%% Common Factors

m_i = 40; % Satellite Inert Mass [kg]

% Planned Mission Duration / Satellite Lifetime
mission_time = 5; % [years]

% Deorbit Altitude --> bring satellite into high air drag zone
r_d = 200 + R_A; % 200 km altitude selected --> high enough air drag

%% Main Plane Delivery Calculations

% Target Orbital Conditions
Rx_a = 415 + R_A; % Orbit radius [km] (semimajor axis (circular))
Rx_e = 0; % Eccentricity
Rx_i = makeSSO(Rx_a); % Inclination [deg.]
Rx_RAAN = 0; % Right Ascension of AN [deg.]
Rx_w = 0; % Argument of Perigee [deg.]
Rx_p = Rx_a*(1 - Rx_e^2); % Semilatus Rectum p [km]

% Initial Delivery Orbit from Launch Vehicle
r_1 = 400 + R_A; % All sats delivered to 400 km altitude orbit

% Hohmann Transfer to Correct Altitude - Timed to space by True Anomaly
dV_1_1 = changeApogee(r_1, Rx_a, MU); % [km/s]
dV_1_2 = circularize(r_1, Rx_a, MU); % [km/s]

% Station Keeping Burns (countering drag per year) - update from simulation
dV_1_s = 0.02; % [km/s]

% Deorbit Burn --> Drop Perigee to 200 km
dV_1_d = changePerigee(Rx_a, r_d, MU);

% Total dV/Propellant Expenditure
I_sp = 60:10:350; % Range of specific impulses [s]
dV_1_total = dV_1_1 + dV_1_2 + dV_1_d + mission_time*dV_1_s;

for n = 1:1:length(I_sp)
    mp_1_total(n) = propExpended(dV_1_total, I_sp(n), m_i, g);
end

figure(1)
plot(I_sp, mp_1_total)
hold on
grid on
xlabel("I_s_p [s]")
ylabel("Propellant Expenditure Over Mission Lifetime [kg]")
title("Propellant Use vs. Thruster I_s_p - Main Plane (24)")

%% Secondary Plane Delivery Calculations - One Launch

% Target Orbital Conditions
Rx_a = 523 + R_A; % Orbit radius [km] (semimajor axis (circular))
Rx_e = 0; % Eccentricity
Rx_i = makeSSO(Rx_a); % Inclination [deg.]
Rx_RAAN = 0; % Right Ascension of AN [deg.]
Rx_w = 0; % Argument of Perigee [deg.]
Rx_p = Rx_a*(1 - Rx_e^2); % Semilatus Rectum p [km]

% Initial Delivery Orbit from Launch Vehicle
r_1 = 400 + R_A; % All sats delivered to 400 km altitude orbit
i_1 = makeSSO(415 + R_A); % All sats delivered into main plane at start

% Hohmann Transfer to Correct Altitude - Timed to space by True Anomaly
dV_2_1 = changeApogee(r_1, Rx_a, MU); % [km/s]
dV_2_2 = circularize(r_1, Rx_a, MU); % [km/s]

% Inclination Plane Change
v_c_2 = velocityCircular(Rx_a, MU); % [km/s]
dV_2_3 = changeInclination((Rx_i - i_1), v_c_2); % [km/s]

% Station Keeping Burns (countering drag per year) - update from simulation
dV_2_s = 0.02; % [km/s]

% Deorbit Burn --> Drop Perigee to 200 km
dV_2_d = changePerigee(Rx_a, r_d, MU);

% Total dV/Propellant Expenditure
dV_2_total = dV_2_1 + dV_2_2 + dV_2_3 + dV_2_d + mission_time*dV_2_s;

for n = 1:1:length(I_sp)
    mp_2_total(n) = propExpended(dV_2_total, I_sp(n), m_i, g);
end

figure(2)
plot(I_sp, mp_2_total)
hold on
grid on
xlabel("I_s_p [s]")
ylabel("Propellant Expenditure Over Mission Lifetime [kg]")
title("Propellant Use vs. Thruster I_s_p - Secondary Plane (12) - One Launch")

%% Secondary Plane Delivery Calculations - Separate Launch

% Target Orbital Conditions
Rx_a = 523 + R_A; % Orbit radius [km] (semimajor axis (circular))
Rx_e = 0; % Eccentricity
Rx_i = makeSSO(Rx_a); % Inclination [deg.]
Rx_RAAN = 0; % Right Ascension of AN [deg.]
Rx_w = 0; % Argument of Perigee [deg.]
Rx_p = Rx_a*(1 - Rx_e^2); % Semilatus Rectum p [km]

% Initial Delivery Orbit from Launch Vehicle
r_1 = 500 + R_A; % Secondary plane sats delivered to 500 km orbit

% Hohmann Transfer to Correct Altitude - Timed to space by True Anomaly
dV_3_1 = changeApogee(r_1, Rx_a, MU); % [km/s]
dV_3_2 = circularize(r_1, Rx_a, MU); % [km/s]

% Inclination Plane Change --> No longer necessary

% Station Keeping Burns (countering drag per year) - update from simulation
dV_3_s = 0.02; % [km/s]

% Deorbit Burn --> Drop Perigee to 200 km
dV_3_d = changePerigee(Rx_a, r_d, MU);

% Total dV/Propellant Expenditure
dV_3_total = dV_3_1 + dV_3_2 + dV_3_d + mission_time*dV_3_s;

for n = 1:1:length(I_sp)
    mp_3_total(n) = propExpended(dV_3_total, I_sp(n), m_i, g);
end

figure(3)
plot(I_sp, mp_3_total)
hold on
grid on
xlabel("I_s_p [s]")
ylabel("Propellant Expenditure Over Mission Lifetime [kg]")
title("Propellant Use vs. Thruster I_s_p - Secondary Plane (12) - Separate Launch")

%% Console Outputs

fprintf("------------------------------------------------------------------------")
fprintf("\nSelected Mission Lifetime: %0.0f years\n", mission_time)
fprintf("Satellite Inert Mass: %0.2f kg\n", m_i)
fprintf("Deorbit Altitude: %0.1f km\n", r_d - R_A)

fprintf("------------------------------------------------------------------------")
fprintf("\n\nMain Plane deltaV per Satellite:\n\n")
fprintf("Hohmann Transfer: %0.3f [m/s]\n", (dV_1_1 + dV_1_2)*1000)
fprintf("Station Keeping: %0.3f [m/s]\n", mission_time*dV_1_s*1000)
fprintf("Deorbit Burn: %0.3f [m/s]\n", dV_1_d*1000)
fprintf("\nTotal dV: %0.3f [m/s]\n\n", dV_1_total*1000)

fprintf("------------------------------------------------------------------------")
fprintf("\n\nSecondary Plane deltaV per Satellite: One Launch\n\n")
fprintf("Hohmann Transfer: %0.3f [m/s]\n", (dV_2_1 + dV_2_2)*1000)
fprintf("Inclination Change: %0.3f [m/s]\n", (dV_2_1 + dV_2_2)*1000)
fprintf("Station Keeping: %0.3f [m/s]\n", mission_time*dV_2_s*1000)
fprintf("Deorbit Burn: %0.3f [m/s]\n", dV_2_d*1000)
fprintf("\nTotal dV: %0.3f [m/s]\n\n", dV_2_total*1000)

fprintf("------------------------------------------------------------------------")
fprintf("\n\nSecondary Plane deltaV per Satellite: Two Launches\n\n")
fprintf("Hohmann Transfer: %0.3f [m/s]\n", (dV_3_1 + dV_3_2)*1000)
fprintf("Station Keeping: %0.3f [m/s]\n", mission_time*dV_3_s*1000)
fprintf("Deorbit Burn: %0.3f [m/s]\n", dV_3_d*1000)
fprintf("\nTotal dV: %0.3f [m/s]\n\n", dV_3_total*1000)

fprintf("------------------------------------------------------------------------\n")

%% Maneuver Calculation Functions

% Altitude dV Changes

% Change Apogee
% r_1 = starting circular orbit radius [km]
% r_2 = final apogee radius [km]
% MU = Gravitational Parameter for Body being Orbited [km^3/s^2]
function dV = changeApogee(r_1, r_2, MU)

    a_minus = r_1;
    a_plus = (r_1 + r_2) / 2;

    v_minus = sqrt(MU / a_minus);
    v_plus = sqrt((2*(sqrt(MU / a_minus))^2) - (MU / a_plus));

    dV = v_plus - v_minus;

end

% Circularize
% r_1 = starting perigee orbit radius [km]
% r_2 = final circular orbit radius [km]
% MU = Gravitational Parameter for Body being Orbited [km^3/s^2]
function dV = circularize(r_1, r_2, MU)

    a_minus = (r_1 + r_2) / 2;
    a_plus = r_2;

    v_plus = sqrt(MU / a_plus);
    v_minus = sqrt((2*(sqrt(MU / a_plus))^2) - (MU / a_minus));

    dV = v_plus - v_minus;

end

% Change Perigee
% r_1 = starting circular orbit radius [km]
% r_2 = final perigee radius [km]
% MU = Gravitational Parameter for Body being Orbited [km^3/s^2]
function dV = changePerigee(r_1, r_2, MU)

    a_minus = r_1;
    a_plus = (r_1 + r_2) / 2;

    v_minus = sqrt(MU / a_minus);
    v_plus = sqrt((2*(sqrt(MU / a_minus))^2) - (MU / a_plus));

    dV = abs(v_plus - v_minus);

end

% Change Inclination
% delta_i = change in inclination [deg.]
% v_current = current circular orbit velocity [km]
function dV = changeInclination(delta_i, v_current)

    dV = 2 * v_current * sind(delta_i / 2);

end

% Propellant Mass Expended
% deltaV = total deltaV needed for mission [km/s]
% I_sp = specific impulse of s/c thruster [s]
% m_i = inert mass of s/c [kg]
% g = gravitational acceleration of body [m/s^2]
function m_p = propExpended(deltaV, I_sp, m_i, g)

    m_p = (1 - exp((-deltaV*1000)/(g*I_sp))) * m_i;

end

% Circular Orbital Velocity
% r = orbital radius [km]
% MU = Gravitational Parameter for Body being Orbited [km^3/s^2]
function v_c = velocityCircular(r, MU)

    v_c = sqrt(MU / r);

end

