% 2nd Order Analysis - Numerical Solution
% Author: V. Swaminathan
% Version: 11/15/2021 1912 EST
% Purpose: To help verify that Team SoUP's satellite constellation meets
%          requirements of AAE450 Project (2nd order analysis)
% 
% To Be Provided By Code:
% - Location/Time Histories of Tx/Rx satellites
% - Relative position to Tx sats from Rx sats
% - "Coverage" requirement verification (soft, mainly for SWE)

%% Initialize Workspace

fprintf("Clearing Workspace...")

clear;
close all;

fprintf("Done!\n")

tic; % Start Run Timer

%% Constant Parameters
R_A = 6378.137; % Earth Equatorial Radius [km]
R_B = 6356.752; % Earth Polar Radius [km]
MU = 3.986e5; % Gravitational Param. for Earth [km^3/s^2]
J2 = 1082.63e-6; % J2 param. for Earth
SideDay_E = 86164; % Sidereal Day Length [s]
SolarDay_E = 86400; % Solar Day Length [s]
W_E = (-0.2507)/60; % Rate of Earth Rotation [deg./s]

%% Satellite Scenario Definition
p_days = 1; % Days of simulation propogation
viewOpt_3D = 0; % 1 = show 3D viewer, 0 = 2D only

startTime = datetime(2020,5,11,12,35,38); % Start Epoch
stopTime = startTime + days(p_days); % End Epoch
sampleTime = 60; % Sample Time [s]

fprintf("Initializing orbital simulation...")

sc_main = satelliteScenario(startTime, stopTime, sampleTime);

fprintf("Done!\n")

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

fprintf("Defining Constellation - Main Plane (24)...")

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

fprintf("Done!\n")

% COMMON FACTORS - Second Plane
Rx_a = (523 + R_A)*1000; % Orbit radius [km] (semimajor axis (circular))
Rx_e = 0; % Eccentricity
Rx_i = makeSSO_NS(Rx_a); % Inclination [deg.]
Rx_RAAN = 0; % Right Ascension of AN [deg.]
Rx_w = 0; % Argument of Perigee [deg.]

Rx_TA = linspace(0,360,numsats_second); % True Anomaly Values [deg.]

% SECONDARY ORBITAL PLANE (12)

fprintf("Defining Constellation - Secondary Plane (12)...")

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

fprintf("Done!\n")

%% Orbital Simulation/Ground Tracks in 3D

if viewOpt_3D == 1
    
    fprintf("Preparing 3D Viewer...")

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
    View_Orb = satelliteScenarioViewer(sc_main);
    
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
    
    fprintf("Done!\n")

    % Visual of Orbital Propagation
    play(sc_main);

else
    fprintf("3D VIEWER DISABLED - ENABLE IN SCENARIO DEFINITION SECTION\n")
end

%% Data Processing/Export

fprintf("Processing Data from Orbit Propagation...")

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

fprintf("Done!\n")

%% Converting to LLA --> Need help speeding up?

fprintf("Converting position data from ECI to LLA Format...\n")
poolobj = parpool(4); % Parallel Pool w/ 4 Threads

fprintf('Progress Bar: %.0f Iterations needed\n', length(pos_Rx_1_1));
fprintf([repmat('.',1,length(pos_Rx_1_1)) '\n\n']);

% Data Conversion - Multithreaded
parfor n = 1:length(pos_Rx_1_1)

    % Main Plane - 24
    lla_1_1(n,:) = eci2lla(pos_Rx_1_1(:,n)', datevec(time_vec(n)));    
    lla_1_2(n,:) = eci2lla(pos_Rx_1_2(:,n)', datevec(time_vec(n)));
    lla_1_3(n,:) = eci2lla(pos_Rx_1_3(:,n)', datevec(time_vec(n)));
    lla_1_4(n,:) = eci2lla(pos_Rx_1_4(:,n)', datevec(time_vec(n)));
    lla_1_5(n,:) = eci2lla(pos_Rx_1_5(:,n)', datevec(time_vec(n)));
    lla_1_6(n,:) = eci2lla(pos_Rx_1_6(:,n)', datevec(time_vec(n)));
    lla_1_7(n,:) = eci2lla(pos_Rx_1_7(:,n)', datevec(time_vec(n)));
    lla_1_8(n,:) = eci2lla(pos_Rx_1_8(:,n)', datevec(time_vec(n)));
    lla_1_9(n,:) = eci2lla(pos_Rx_1_9(:,n)', datevec(time_vec(n)));
    lla_1_10(n,:) = eci2lla(pos_Rx_1_10(:,n)', datevec(time_vec(n)));
    lla_1_11(n,:) = eci2lla(pos_Rx_1_11(:,n)', datevec(time_vec(n)));
    lla_1_12(n,:) = eci2lla(pos_Rx_1_12(:,n)', datevec(time_vec(n)));
    lla_1_13(n,:) = eci2lla(pos_Rx_1_13(:,n)', datevec(time_vec(n)));
    lla_1_14(n,:) = eci2lla(pos_Rx_1_14(:,n)', datevec(time_vec(n)));
    lla_1_15(n,:) = eci2lla(pos_Rx_1_15(:,n)', datevec(time_vec(n)));
    lla_1_16(n,:) = eci2lla(pos_Rx_1_16(:,n)', datevec(time_vec(n)));
    lla_1_17(n,:) = eci2lla(pos_Rx_1_17(:,n)', datevec(time_vec(n)));
    lla_1_18(n,:) = eci2lla(pos_Rx_1_18(:,n)', datevec(time_vec(n)));
    lla_1_19(n,:) = eci2lla(pos_Rx_1_19(:,n)', datevec(time_vec(n)));
    lla_1_20(n,:) = eci2lla(pos_Rx_1_20(:,n)', datevec(time_vec(n)));
    lla_1_21(n,:) = eci2lla(pos_Rx_1_21(:,n)', datevec(time_vec(n)));
    lla_1_22(n,:) = eci2lla(pos_Rx_1_22(:,n)', datevec(time_vec(n)));
    lla_1_23(n,:) = eci2lla(pos_Rx_1_23(:,n)', datevec(time_vec(n)));
    lla_1_24(n,:) = eci2lla(pos_Rx_1_24(:,n)', datevec(time_vec(n)));
    
    % Secondary Plane - 12
    lla_2_1(n,:) = eci2lla(pos_Rx_2_1(:,n)', datevec(time_vec(n)));    
    lla_2_2(n,:) = eci2lla(pos_Rx_2_2(:,n)', datevec(time_vec(n)));
    lla_2_3(n,:) = eci2lla(pos_Rx_2_3(:,n)', datevec(time_vec(n)));
    lla_2_4(n,:) = eci2lla(pos_Rx_2_4(:,n)', datevec(time_vec(n)));
    lla_2_5(n,:) = eci2lla(pos_Rx_2_5(:,n)', datevec(time_vec(n)));
    lla_2_6(n,:) = eci2lla(pos_Rx_2_6(:,n)', datevec(time_vec(n)));
    lla_2_7(n,:) = eci2lla(pos_Rx_2_7(:,n)', datevec(time_vec(n)));
    lla_2_8(n,:) = eci2lla(pos_Rx_2_8(:,n)', datevec(time_vec(n)));
    lla_2_9(n,:) = eci2lla(pos_Rx_2_9(:,n)', datevec(time_vec(n)));
    lla_2_10(n,:) = eci2lla(pos_Rx_2_10(:,n)', datevec(time_vec(n)));
    lla_2_11(n,:) = eci2lla(pos_Rx_2_11(:,n)', datevec(time_vec(n)));
    lla_2_12(n,:) = eci2lla(pos_Rx_2_12(:,n)', datevec(time_vec(n)));

    % Loop Visualization - Out-of-order Execution
    %fprintf("Iteration %.0f / %.0f completed\n", n, length(pos_Rx_1_1));
    fprintf('\b|\n');

end
delete(poolobj); % Close parallel processing pool

%% Coverage of the Entire World

fprintf("\nPlotting World Map Coverage...")
figure(2)
worldmap world
load coastlines
[latcells, loncells] = polysplit(coastlat, coastlon);
plotm(coastlat, coastlon, 'black')
hold on
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
geoshow('worldlakes.shp', 'FaceColor', 'cyan')
geoshow('worldrivers.shp', 'Color', 'blue')
title("Worldwide Coverage")

% Main Plane - 24
plotm(lla_1_1(:,1), lla_1_1(:,2), 'red')
plotm(lla_1_2(:,1), lla_1_2(:,2), 'red')
plotm(lla_1_3(:,1), lla_1_3(:,2), 'red')
plotm(lla_1_4(:,1), lla_1_4(:,2), 'red')
plotm(lla_1_5(:,1), lla_1_5(:,2), 'red')
plotm(lla_1_6(:,1), lla_1_6(:,2), 'red')
plotm(lla_1_7(:,1), lla_1_7(:,2), 'red')
plotm(lla_1_8(:,1), lla_1_8(:,2), 'red')
plotm(lla_1_9(:,1), lla_1_9(:,2), 'red')
plotm(lla_1_10(:,1), lla_1_10(:,2), 'red')
plotm(lla_1_11(:,1), lla_1_11(:,2), 'red')
plotm(lla_1_12(:,1), lla_1_12(:,2), 'red')
plotm(lla_1_13(:,1), lla_1_13(:,2), 'red')
plotm(lla_1_14(:,1), lla_1_14(:,2), 'red')
plotm(lla_1_15(:,1), lla_1_15(:,2), 'red')
plotm(lla_1_16(:,1), lla_1_16(:,2), 'red')
plotm(lla_1_17(:,1), lla_1_17(:,2), 'red')
plotm(lla_1_18(:,1), lla_1_18(:,2), 'red')
plotm(lla_1_19(:,1), lla_1_19(:,2), 'red')
plotm(lla_1_20(:,1), lla_1_20(:,2), 'red')
plotm(lla_1_21(:,1), lla_1_21(:,2), 'red')
plotm(lla_1_22(:,1), lla_1_22(:,2), 'red')
plotm(lla_1_23(:,1), lla_1_23(:,2), 'red')
plotm(lla_1_24(:,1), lla_1_24(:,2), 'red')

% Secondary Plane - 12
plotm(lla_2_1(:,1), lla_2_1(:,2), 'blue')
plotm(lla_2_2(:,1), lla_2_2(:,2), 'blue')
plotm(lla_2_3(:,1), lla_2_3(:,2), 'blue')
plotm(lla_2_4(:,1), lla_2_4(:,2), 'blue')
plotm(lla_2_5(:,1), lla_2_5(:,2), 'blue')
plotm(lla_2_6(:,1), lla_2_6(:,2), 'blue')
plotm(lla_2_7(:,1), lla_2_7(:,2), 'blue')
plotm(lla_2_8(:,1), lla_2_8(:,2), 'blue')
plotm(lla_2_9(:,1), lla_2_9(:,2), 'blue')
plotm(lla_2_10(:,1), lla_2_10(:,2), 'blue')
plotm(lla_2_11(:,1), lla_2_11(:,2), 'blue')
plotm(lla_2_12(:,1), lla_2_12(:,2), 'blue')

fprintf("Done!\n")

%% Coverage of Europe

fprintf("Plotting Europe Map Coverage...")
figure(3)
h_1 = worldmap('Europe');
getm(h_1,"MapProjection");
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
geoshow('worldlakes.shp', 'FaceColor', 'cyan')
geoshow('worldrivers.shp', 'Color', 'blue')
geoshow('worldcities.shp', 'Marker', '.',...
                           'MarkerEdgeColor', 'magenta')
title("Coverage of Europe")

% Main Plane - 24
plotm(lla_1_1(:,1), lla_1_1(:,2), 'red')
plotm(lla_1_2(:,1), lla_1_2(:,2), 'red')
plotm(lla_1_3(:,1), lla_1_3(:,2), 'red')
plotm(lla_1_4(:,1), lla_1_4(:,2), 'red')
plotm(lla_1_5(:,1), lla_1_5(:,2), 'red')
plotm(lla_1_6(:,1), lla_1_6(:,2), 'red')
plotm(lla_1_7(:,1), lla_1_7(:,2), 'red')
plotm(lla_1_8(:,1), lla_1_8(:,2), 'red')
plotm(lla_1_9(:,1), lla_1_9(:,2), 'red')
plotm(lla_1_10(:,1), lla_1_10(:,2), 'red')
plotm(lla_1_11(:,1), lla_1_11(:,2), 'red')
plotm(lla_1_12(:,1), lla_1_12(:,2), 'red')
plotm(lla_1_13(:,1), lla_1_13(:,2), 'red')
plotm(lla_1_14(:,1), lla_1_14(:,2), 'red')
plotm(lla_1_15(:,1), lla_1_15(:,2), 'red')
plotm(lla_1_16(:,1), lla_1_16(:,2), 'red')
plotm(lla_1_17(:,1), lla_1_17(:,2), 'red')
plotm(lla_1_18(:,1), lla_1_18(:,2), 'red')
plotm(lla_1_19(:,1), lla_1_19(:,2), 'red')
plotm(lla_1_20(:,1), lla_1_20(:,2), 'red')
plotm(lla_1_21(:,1), lla_1_21(:,2), 'red')
plotm(lla_1_22(:,1), lla_1_22(:,2), 'red')
plotm(lla_1_23(:,1), lla_1_23(:,2), 'red')
plotm(lla_1_24(:,1), lla_1_24(:,2), 'red')

% Secondary Plane - 12
plotm(lla_2_1(:,1), lla_2_1(:,2), 'blue')
plotm(lla_2_2(:,1), lla_2_2(:,2), 'blue')
plotm(lla_2_3(:,1), lla_2_3(:,2), 'blue')
plotm(lla_2_4(:,1), lla_2_4(:,2), 'blue')
plotm(lla_2_5(:,1), lla_2_5(:,2), 'blue')
plotm(lla_2_6(:,1), lla_2_6(:,2), 'blue')
plotm(lla_2_7(:,1), lla_2_7(:,2), 'blue')
plotm(lla_2_8(:,1), lla_2_8(:,2), 'blue')
plotm(lla_2_9(:,1), lla_2_9(:,2), 'blue')
plotm(lla_2_10(:,1), lla_2_10(:,2), 'blue')
plotm(lla_2_11(:,1), lla_2_11(:,2), 'blue')
plotm(lla_2_12(:,1), lla_2_12(:,2), 'blue')

fprintf("Done!\n")

%% Coverage of USA

fprintf("Plotting USA Map Coverage...")
figure(4)
h_2 = worldmap('USA');
getm(h_2,"MapProjection");
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
geoshow('worldlakes.shp', 'FaceColor', 'cyan')
geoshow('worldrivers.shp', 'Color', 'blue')
geoshow('worldcities.shp', 'Marker', '.',...
                           'MarkerEdgeColor', 'magenta')
title("Coverage of USA")

% Main Plane - 24
plotm(lla_1_1(:,1), lla_1_1(:,2), 'red')
plotm(lla_1_2(:,1), lla_1_2(:,2), 'red')
plotm(lla_1_3(:,1), lla_1_3(:,2), 'red')
plotm(lla_1_4(:,1), lla_1_4(:,2), 'red')
plotm(lla_1_5(:,1), lla_1_5(:,2), 'red')
plotm(lla_1_6(:,1), lla_1_6(:,2), 'red')
plotm(lla_1_7(:,1), lla_1_7(:,2), 'red')
plotm(lla_1_8(:,1), lla_1_8(:,2), 'red')
plotm(lla_1_9(:,1), lla_1_9(:,2), 'red')
plotm(lla_1_10(:,1), lla_1_10(:,2), 'red')
plotm(lla_1_11(:,1), lla_1_11(:,2), 'red')
plotm(lla_1_12(:,1), lla_1_12(:,2), 'red')
plotm(lla_1_13(:,1), lla_1_13(:,2), 'red')
plotm(lla_1_14(:,1), lla_1_14(:,2), 'red')
plotm(lla_1_15(:,1), lla_1_15(:,2), 'red')
plotm(lla_1_16(:,1), lla_1_16(:,2), 'red')
plotm(lla_1_17(:,1), lla_1_17(:,2), 'red')
plotm(lla_1_18(:,1), lla_1_18(:,2), 'red')
plotm(lla_1_19(:,1), lla_1_19(:,2), 'red')
plotm(lla_1_20(:,1), lla_1_20(:,2), 'red')
plotm(lla_1_21(:,1), lla_1_21(:,2), 'red')
plotm(lla_1_22(:,1), lla_1_22(:,2), 'red')
plotm(lla_1_23(:,1), lla_1_23(:,2), 'red')
plotm(lla_1_24(:,1), lla_1_24(:,2), 'red')

% Secondary Plane - 12
plotm(lla_2_1(:,1), lla_2_1(:,2), 'blue')
plotm(lla_2_2(:,1), lla_2_2(:,2), 'blue')
plotm(lla_2_3(:,1), lla_2_3(:,2), 'blue')
plotm(lla_2_4(:,1), lla_2_4(:,2), 'blue')
plotm(lla_2_5(:,1), lla_2_5(:,2), 'blue')
plotm(lla_2_6(:,1), lla_2_6(:,2), 'blue')
plotm(lla_2_7(:,1), lla_2_7(:,2), 'blue')
plotm(lla_2_8(:,1), lla_2_8(:,2), 'blue')
plotm(lla_2_9(:,1), lla_2_9(:,2), 'blue')
plotm(lla_2_10(:,1), lla_2_10(:,2), 'blue')
plotm(lla_2_11(:,1), lla_2_11(:,2), 'blue')
plotm(lla_2_12(:,1), lla_2_12(:,2), 'blue')

fprintf("Done!\n")

%% Coverage of 100km x 100km grid at Equator

fprintf("Plotting 100x100km @ Equator Coverage...")
figure(5)
h_3 = worldmap([-.45 0.45], [-.45 0.45]);
getm(h_3,"MapProjection");
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
geoshow('worldlakes.shp', 'FaceColor', 'cyan')
geoshow('worldrivers.shp', 'Color', 'blue')
geoshow('worldcities.shp', 'Marker', '.',...
                           'MarkerEdgeColor', 'magenta')
title("100x100 km at Equator")

% Main Plane - 24
plotm(lla_1_1(:,1), lla_1_1(:,2), 'red')
plotm(lla_1_2(:,1), lla_1_2(:,2), 'red')
plotm(lla_1_3(:,1), lla_1_3(:,2), 'red')
plotm(lla_1_4(:,1), lla_1_4(:,2), 'red')
plotm(lla_1_5(:,1), lla_1_5(:,2), 'red')
plotm(lla_1_6(:,1), lla_1_6(:,2), 'red')
plotm(lla_1_7(:,1), lla_1_7(:,2), 'red')
plotm(lla_1_8(:,1), lla_1_8(:,2), 'red')
plotm(lla_1_9(:,1), lla_1_9(:,2), 'red')
plotm(lla_1_10(:,1), lla_1_10(:,2), 'red')
plotm(lla_1_11(:,1), lla_1_11(:,2), 'red')
plotm(lla_1_12(:,1), lla_1_12(:,2), 'red')
plotm(lla_1_13(:,1), lla_1_13(:,2), 'red')
plotm(lla_1_14(:,1), lla_1_14(:,2), 'red')
plotm(lla_1_15(:,1), lla_1_15(:,2), 'red')
plotm(lla_1_16(:,1), lla_1_16(:,2), 'red')
plotm(lla_1_17(:,1), lla_1_17(:,2), 'red')
plotm(lla_1_18(:,1), lla_1_18(:,2), 'red')
plotm(lla_1_19(:,1), lla_1_19(:,2), 'red')
plotm(lla_1_20(:,1), lla_1_20(:,2), 'red')
plotm(lla_1_21(:,1), lla_1_21(:,2), 'red')
plotm(lla_1_22(:,1), lla_1_22(:,2), 'red')
plotm(lla_1_23(:,1), lla_1_23(:,2), 'red')
plotm(lla_1_24(:,1), lla_1_24(:,2), 'red')

% Secondary Plane - 12
plotm(lla_2_1(:,1), lla_2_1(:,2), 'blue')
plotm(lla_2_2(:,1), lla_2_2(:,2), 'blue')
plotm(lla_2_3(:,1), lla_2_3(:,2), 'blue')
plotm(lla_2_4(:,1), lla_2_4(:,2), 'blue')
plotm(lla_2_5(:,1), lla_2_5(:,2), 'blue')
plotm(lla_2_6(:,1), lla_2_6(:,2), 'blue')
plotm(lla_2_7(:,1), lla_2_7(:,2), 'blue')
plotm(lla_2_8(:,1), lla_2_8(:,2), 'blue')
plotm(lla_2_9(:,1), lla_2_9(:,2), 'blue')
plotm(lla_2_10(:,1), lla_2_10(:,2), 'blue')
plotm(lla_2_11(:,1), lla_2_11(:,2), 'blue')
plotm(lla_2_12(:,1), lla_2_12(:,2), 'blue')

fprintf("Done!\n")

%% Coverage of 100km x 100km grid at Latitude = 70 deg.

fprintf("Plotting 100x100km @ Lat. = 70 deg. Coverage...")
figure(6)
h_4 = worldmap([69.55 70.45], [-.45 0.45]);
getm(h_4,"MapProjection");
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
geoshow('worldlakes.shp', 'FaceColor', 'cyan')
geoshow('worldrivers.shp', 'Color', 'blue')
geoshow('worldcities.shp', 'Marker', '.',...
                           'MarkerEdgeColor', 'magenta')
title("100x100 km at Lat = 70 deg.")

% Main Plane - 24
plotm(lla_1_1(:,1), lla_1_1(:,2), 'red')
plotm(lla_1_2(:,1), lla_1_2(:,2), 'red')
plotm(lla_1_3(:,1), lla_1_3(:,2), 'red')
plotm(lla_1_4(:,1), lla_1_4(:,2), 'red')
plotm(lla_1_5(:,1), lla_1_5(:,2), 'red')
plotm(lla_1_6(:,1), lla_1_6(:,2), 'red')
plotm(lla_1_7(:,1), lla_1_7(:,2), 'red')
plotm(lla_1_8(:,1), lla_1_8(:,2), 'red')
plotm(lla_1_9(:,1), lla_1_9(:,2), 'red')
plotm(lla_1_10(:,1), lla_1_10(:,2), 'red')
plotm(lla_1_11(:,1), lla_1_11(:,2), 'red')
plotm(lla_1_12(:,1), lla_1_12(:,2), 'red')
plotm(lla_1_13(:,1), lla_1_13(:,2), 'red')
plotm(lla_1_14(:,1), lla_1_14(:,2), 'red')
plotm(lla_1_15(:,1), lla_1_15(:,2), 'red')
plotm(lla_1_16(:,1), lla_1_16(:,2), 'red')
plotm(lla_1_17(:,1), lla_1_17(:,2), 'red')
plotm(lla_1_18(:,1), lla_1_18(:,2), 'red')
plotm(lla_1_19(:,1), lla_1_19(:,2), 'red')
plotm(lla_1_20(:,1), lla_1_20(:,2), 'red')
plotm(lla_1_21(:,1), lla_1_21(:,2), 'red')
plotm(lla_1_22(:,1), lla_1_22(:,2), 'red')
plotm(lla_1_23(:,1), lla_1_23(:,2), 'red')
plotm(lla_1_24(:,1), lla_1_24(:,2), 'red')

% Secondary Plane - 12
plotm(lla_2_1(:,1), lla_2_1(:,2), 'blue')
plotm(lla_2_2(:,1), lla_2_2(:,2), 'blue')
plotm(lla_2_3(:,1), lla_2_3(:,2), 'blue')
plotm(lla_2_4(:,1), lla_2_4(:,2), 'blue')
plotm(lla_2_5(:,1), lla_2_5(:,2), 'blue')
plotm(lla_2_6(:,1), lla_2_6(:,2), 'blue')
plotm(lla_2_7(:,1), lla_2_7(:,2), 'blue')
plotm(lla_2_8(:,1), lla_2_8(:,2), 'blue')
plotm(lla_2_9(:,1), lla_2_9(:,2), 'blue')
plotm(lla_2_10(:,1), lla_2_10(:,2), 'blue')
plotm(lla_2_11(:,1), lla_2_11(:,2), 'blue')
plotm(lla_2_12(:,1), lla_2_12(:,2), 'blue')

fprintf("Done!\n")

toc % Display total runtime of NumericalSolver_2OA
