% 2nd Order Analysis - Numerical Solution
% Author: V. Swaminathan
% Version: 11/09/2021 1054 EST
% Purpose: To help verify that Team SoUP's satellite constellation meets
%          requirements of AAE450 Project (2nd order analysis)
% 
% To Be Provided By Code:
% - Location/Time Histories of Tx/Rx satellites
% - Relative position to Tx sats from Rx sats
% - "Coverage" requirement verification (soft, mainly for SWE)

%% Initialize Workspace
clear all;
close all;

%% Constant Parameters
R_A = 6378.137; % Earth Equatorial Radius [km]
R_B = 6356.752; % Earth Polar Radius [km]
MU = 3.986e5; % Gravitational Param. for Earth [km^3/s^2]
J2 = 1082.63e-6; % J2 param. for Earth
SideDay_E = 86164; % Sidereal Day Length [s]
SolarDay_E = 86400; % Solar Day Length [s]
W_E = (-0.2507)/60; % Rate of Earth Rotation [deg./s]

%% Satellite Scenario Definition
p_days = 0.25; % Days of simulation propogation

startTime = datetime(2020,5,11,12,35,38); % Start Epoch
stopTime = startTime + days(p_days); % End Epoch
sampleTime = 60; % Sample Time [s]

sc_main = satelliteScenario(startTime, stopTime, sampleTime);

%% Tx Satellite Import

% MUOS

% MERIDIAN

% SKYNET

% ORBCOMM

% IRIDIUM

% GPS

% GLONASS

% BEIDOU

% GALLILEO

% SWARM


%% Rx Constellation Definitions

% Number of Satellites in Constellation
numsats_main = 24;
numsats_second = 12;

% COMMON FACTORS - Main Plane
Rx_a = (415 + R_A)*1000; % Orbit radius [km] (semimajor axis (circular))
Rx_e = 0; % Eccentricity
Rx_i = makeSSO_NS(Rx_a); % Inclination [deg.]
Rx_RAAN = 0; % Right Ascension of AN [deg.]
Rx_w = 0; % Argument of Perigee [deg.]

Rx_TA = linspace(0,360,numsats_main); % True Anomaly Values [deg.]

% PRIMARY ORBITAL PLANE (24)
Rx_1_1 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(1), "OrbitPropagator","sgp4");
Rx_1_2 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(2), "OrbitPropagator","sgp4");
Rx_1_3 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(3), "OrbitPropagator","sgp4");
Rx_1_4 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(4), "OrbitPropagator","sgp4");
Rx_1_5 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(5), "OrbitPropagator","sgp4");
Rx_1_6 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(6), "OrbitPropagator","sgp4");
Rx_1_7 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(7), "OrbitPropagator","sgp4");
Rx_1_8 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(8), "OrbitPropagator","sgp4");
Rx_1_9 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(9), "OrbitPropagator","sgp4");
Rx_1_10 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(10), "OrbitPropagator","sgp4");
Rx_1_11 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(11), "OrbitPropagator","sgp4");
Rx_1_12 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(12), "OrbitPropagator","sgp4");
Rx_1_13 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(13), "OrbitPropagator","sgp4");
Rx_1_14 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(14), "OrbitPropagator","sgp4");
Rx_1_15 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(15), "OrbitPropagator","sgp4");
Rx_1_16 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(16), "OrbitPropagator","sgp4");
Rx_1_17 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(17), "OrbitPropagator","sgp4");
Rx_1_18 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(18), "OrbitPropagator","sgp4");
Rx_1_19 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(19), "OrbitPropagator","sgp4");
Rx_1_20 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(20), "OrbitPropagator","sgp4");
Rx_1_21 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(21), "OrbitPropagator","sgp4");
Rx_1_22 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(22), "OrbitPropagator","sgp4");
Rx_1_23 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(23), "OrbitPropagator","sgp4");
Rx_1_24 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(24), "OrbitPropagator","sgp4");

% COMMON FACTORS - Second Plane
Rx_a = (523 + R_A)*1000; % Orbit radius [km] (semimajor axis (circular))
Rx_e = 0; % Eccentricity
Rx_i = makeSSO_NS(Rx_a); % Inclination [deg.]
Rx_RAAN = 0; % Right Ascension of AN [deg.]
Rx_w = 0; % Argument of Perigee [deg.]

Rx_TA = linspace(0,360,numsats_second); % True Anomaly Values [deg.]

% SECONDARY ORBITAL PLANE (12)
Rx_2_1 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(1), "OrbitPropagator","sgp4");
Rx_2_2 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(2), "OrbitPropagator","sgp4");
Rx_2_3 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(3), "OrbitPropagator","sgp4");
Rx_2_4 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(4), "OrbitPropagator","sgp4");
Rx_2_5 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(5), "OrbitPropagator","sgp4");
Rx_2_6 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(6), "OrbitPropagator","sgp4");
Rx_2_7 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(7), "OrbitPropagator","sgp4");
Rx_2_8 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(8), "OrbitPropagator","sgp4");
Rx_2_9 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(9), "OrbitPropagator","sgp4");
Rx_2_10 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(10), "OrbitPropagator","sgp4");
Rx_2_11 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(11), "OrbitPropagator","sgp4");
Rx_2_12 = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(12), "OrbitPropagator","sgp4");

%% Orbital Simulation/Ground Tracks

% Viewer Options - Colors (Change Secondary Plane)
Rx_2_1.MarkerColor = [0 1 0];
Rx_2_1.Orbit.LineColor = [0 1 0];
Rx_2_1.LabelFontColor = [0 1 0];

Rx_2_2.MarkerColor = [0 1 0];
Rx_2_2.Orbit.LineColor = [0 1 0];
Rx_2_2.LabelFontColor = [0 1 0];

Rx_2_3.MarkerColor = [0 1 0];
Rx_2_3.Orbit.LineColor = [0 1 0];
Rx_2_3.LabelFontColor = [0 1 0];

Rx_2_4.MarkerColor = [0 1 0];
Rx_2_4.Orbit.LineColor = [0 1 0];
Rx_2_4.LabelFontColor = [0 1 0];

Rx_2_5.MarkerColor = [0 1 0];
Rx_2_5.Orbit.LineColor = [0 1 0];
Rx_2_5.LabelFontColor = [0 1 0];

Rx_2_6.MarkerColor = [0 1 0];
Rx_2_6.Orbit.LineColor = [0 1 0];
Rx_2_6.LabelFontColor = [0 1 0];

Rx_2_7.MarkerColor = [0 1 0];
Rx_2_7.Orbit.LineColor = [0 1 0];
Rx_2_7.LabelFontColor = [0 1 0];

Rx_2_8.MarkerColor = [0 1 0];
Rx_2_8.Orbit.LineColor = [0 1 0];
Rx_2_8.LabelFontColor = [0 1 0];

Rx_2_9.MarkerColor = [0 1 0];
Rx_2_9.Orbit.LineColor = [0 1 0];
Rx_2_9.LabelFontColor = [0 1 0];

Rx_2_10.MarkerColor = [0 1 0];
Rx_2_10.Orbit.LineColor = [0 1 0];
Rx_2_10.LabelFontColor = [0 1 0];

Rx_2_11.MarkerColor = [0 1 0];
Rx_2_11.Orbit.LineColor = [0 1 0];
Rx_2_11.LabelFontColor = [0 1 0];

Rx_2_12.MarkerColor = [0 1 0];
Rx_2_12.Orbit.LineColor = [0 1 0];
Rx_2_12.LabelFontColor = [0 1 0];

% Initialize Orbital Scenario Viewer
View_Orb = satelliteScenarioViewer(sc_main, "Dimension","2D");

% Initialize Ground Tracks
groundTrack(Rx_1_1, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_2, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_3, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_4, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_5, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_6, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_7, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_8, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_9, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_10, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_11, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_12, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_13, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_14, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_15, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_16, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_17, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_18, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_19, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_20, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_21, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_22, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_23, "LeadTime", p_days*24*60*60);
groundTrack(Rx_1_24, "LeadTime", p_days*24*60*60);

groundTrack(Rx_2_1, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_2, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_3, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_4, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_5, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_6, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_7, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_8, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_9, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_10, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_11, "LeadTime", p_days*24*60*60);
groundTrack(Rx_2_12, "LeadTime", p_days*24*60*60);

% Visual of Orbital Propagation
play(sc_main);

%% Data Processing/Export

% Tx Location/Time History

% Rx Location/Time History - MAIN PLANE
[pos_Rx_1_1, vel_Rx_1_1, time_vec] = states(Rx_1_1);
[pos_Rx_1_2, vel_Rx_1_2] = states(Rx_1_2);
[pos_Rx_1_3, vel_Rx_1_3] = states(Rx_1_3);
[pos_Rx_1_4, vel_Rx_1_4] = states(Rx_1_4);
[pos_Rx_1_5, vel_Rx_1_5] = states(Rx_1_5);
[pos_Rx_1_6, vel_Rx_1_6] = states(Rx_1_6);
[pos_Rx_1_7, vel_Rx_1_7] = states(Rx_1_7);
[pos_Rx_1_8, vel_Rx_1_8] = states(Rx_1_8);
[pos_Rx_1_9, vel_Rx_1_9] = states(Rx_1_9);
[pos_Rx_1_10, vel_Rx_1_10] = states(Rx_1_10);
[pos_Rx_1_11, vel_Rx_1_11] = states(Rx_1_11);
[pos_Rx_1_12, vel_Rx_1_12] = states(Rx_1_12);
[pos_Rx_1_13, vel_Rx_1_13] = states(Rx_1_13);
[pos_Rx_1_14, vel_Rx_1_14] = states(Rx_1_14);
[pos_Rx_1_15, vel_Rx_1_15] = states(Rx_1_15);
[pos_Rx_1_16, vel_Rx_1_16] = states(Rx_1_16);
[pos_Rx_1_17, vel_Rx_1_17] = states(Rx_1_17);
[pos_Rx_1_18, vel_Rx_1_18] = states(Rx_1_18);
[pos_Rx_1_19, vel_Rx_1_19] = states(Rx_1_19);
[pos_Rx_1_20, vel_Rx_1_20] = states(Rx_1_20);
[pos_Rx_1_21, vel_Rx_1_21] = states(Rx_1_21);
[pos_Rx_1_22, vel_Rx_1_22] = states(Rx_1_22);
[pos_Rx_1_23, vel_Rx_1_23] = states(Rx_1_23);
[pos_Rx_1_24, vel_Rx_1_24] = states(Rx_1_24);

% Rx Location/Time History - SECONDARY PLANE
[pos_Rx_2_1, vel_Rx_2_1, time] = states(Rx_2_1);
[pos_Rx_2_2, vel_Rx_2_2] = states(Rx_2_2);
[pos_Rx_2_3, vel_Rx_2_3] = states(Rx_2_3);
[pos_Rx_2_4, vel_Rx_2_4] = states(Rx_2_4);
[pos_Rx_2_5, vel_Rx_2_5] = states(Rx_2_5);
[pos_Rx_2_6, vel_Rx_2_6] = states(Rx_2_6);
[pos_Rx_2_7, vel_Rx_2_7] = states(Rx_2_7);
[pos_Rx_2_8, vel_Rx_2_8] = states(Rx_2_8);
[pos_Rx_2_9, vel_Rx_2_9] = states(Rx_2_9);
[pos_Rx_2_10, vel_Rx_2_10] = states(Rx_2_10);
[pos_Rx_2_11, vel_Rx_2_11] = states(Rx_2_11);
[pos_Rx_2_12, vel_Rx_2_12] = states(Rx_2_12);

% Tx/Rx Relative Position History




