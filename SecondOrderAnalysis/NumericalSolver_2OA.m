% 2nd Order Analysis - Numerical Solution
% Author: V. Swaminathan
% Version: 11/19/2021 1413 EST
% Purpose: To help verify that Team SoUP's satellite constellation meets
%          requirements of AAE450 Project (2nd order analysis)
% 
% To Be Provided By Code:
% - Location/Time Histories of Rx satellites
% - Ground Station Contact Times of Rx satellites 
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

p_days = 15; % Days of simulation propogation

% Options:
% - on --> Show 3D simulation
% - off --> Skip 3D viewer
viewOpt_3D = "off"; % Select 3D Viewer option

% Available Providers:
% - none --> Skip Ground Station Analysis Entirely
% - SSC --> Swedish Space Corporation
% - NEN --> Near Earth Network (NASA) --> NOT IMPLEMENTED YET
% - AWS --> Amazon Web Services --> NOT IMPLEMENTED YET
% - KSAT --> Kongsberg Satellite Services --> NOT IMPLEMENTED YET
ground_network = "none"; % Select Ground Station Provider

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

%% Ground Station Networks

% No Network Selected --> Run w/o Ground Station Analysis
if ground_network == "none"
    fprintf("GROUND STATIONS DISABLED - ENABLE IN SCENARIO DEFINITION SECTION\n")
end

% Near Earth Network (NEN)

% Swedish Space Corporation (SSC)
if ground_network == "SSC"

    fprintf("Defining Ground Station Network...")

    GS_SSC(1) = groundStation(sc_main, 19.0166, -155.6667, "Name", "SSC - South Point Satellite Station", "MinElevationAngle", 10);
    GS_SSC(2) = groundStation(sc_main, 26.7333, -81.00333, "Name", "SSC - Clewiston Satellite Station", "MinElevationAngle", 10);
    GS_SSC(3) = groundStation(sc_main, 64.8000, -147.6500, "Name", "SSC - North Pole Satellite Station", "MinElevationAngle", 10);
    GS_SSC(4) = groundStation(sc_main, 68.4000, -133.5000, "Name", "SSC - Inuvik Satellite Station", "MinElevationAngle", 10);
    GS_SSC(5) = groundStation(sc_main, -33.0133, -70.6666, "Name", "SSC - Santiago Satellite Station", "MinElevationAngle", 10);
    GS_SSC(6) = groundStation(sc_main, -52.9333, -70.8500, "Name", "SSC - Punta Arenas Satellite Station", "MinElevationAngle", 10);
    GS_SSC(7) = groundStation(sc_main, 67.8833, 21.00666, "Name", "SSC - Esrange Space Center", "MinElevationAngle", 10);
    GS_SSC(8) = groundStation(sc_main, 59.2000, 18.00833, "Name", "SSC - Stockholm Teleport", "MinElevationAngle", 10);
    GS_SSC(9) = groundStation(sc_main, 13.1000, 100.9167, "Name", "SSC - Siracha Satellite Station", "MinElevationAngle", 10);
    GS_SSC(10) = groundStation(sc_main, -29.0083, 135.5833, "Name", "SSC - WASC", "MinElevationAngle", 10);
   
    fprintf("Done!\n")

end

% Amazon Web Services (AWS)

% Kongsberg Satellite Services (KSAT)


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

for n = 1:1:24
    name_temp = "SoUP-1-" + n;
    Rx_1(n) = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(n), "OrbitPropagator","sgp4", "Name", name_temp);
end

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

for n = 1:1:12
    name_temp = "SoUP-2-" + n;
    Rx_2(n) = satellite(sc_main,Rx_a,Rx_e,Rx_i,Rx_RAAN,Rx_w,Rx_TA(n), "OrbitPropagator","sgp4", "Name", name_temp);
end

fprintf("Done!\n")

%% Ground Station Contact History

if ground_network == "SSC"

    fprintf("Calculating Ground Station Contacts w/ SSC - Main Plane (24)...")
    
    % Main Plane Access - SSC
    for m = 1:1:length(GS_SSC)
        for n = 1:1:24
            AC_1(n, m) = access(Rx_1(n), GS_SSC(m));
        end
    end
    
    fprintf("Done!\n")
    fprintf("Calculating Ground Station Contacts w/ SSC - Second Plane (12)...")
    
    % Secondary Plane Access - SSC
    for m = 1:1:length(GS_SSC)
        for n = 1:1:12
            AC_2(n, m) = access(Rx_2(n), GS_SSC(m));
        end
    end
    
    fprintf("Done!\n")
    fprintf("Calculating Total Time in Contact for Each Satellite...")
    
    CT_1_temp = zeros(24);
    CT_2_temp = zeros(12);
    for m = 1:1:length(GS_SSC)
        for n = 1:1:24
            CT_1_temp(n) = CT_1_temp(n) + sum(accessIntervals(AC_1(n,m)).Duration(:));
        end
        for n = 1:1:12
            CT_2_temp(n) = CT_2_temp(n) + sum(accessIntervals(AC_2(n,m)).Duration(:));
        end
    end
    CT_1 = CT_1_temp(:,1)';
    CT_2 = CT_2_temp(:,1)';
    
    fprintf("Done!\n")

end

%% Orbital Simulation/Ground Tracks in 3D

if viewOpt_3D == "on"
    
    fprintf("Preparing 3D Viewer...")

    % Viewer Options - Colors (Change Secondary Plane)
    for n = 1:1:12
        Rx_2(n).MarkerColor = [0 1 0];
        Rx_2(n).Orbit.LineColor = [0 1 0];
        Rx_2(n).LabelFontColor = [0 1 0];
    end
    
    % Initialize Orbital Scenario Viewer
    View_Orb = satelliteScenarioViewer(sc_main, "CameraReferenceFrame","Inertial");
    
    % Initialize Ground Tracks
    for n = 1:1:24
        groundTrack(Rx_1(n), "LeadTime", p_days*24*3600);
    end

    for n = 1:1:12
        groundTrack(Rx_2(n), "LeadTime", p_days*24*3600);
    end
    
    fprintf("Done!\n")

    % Visual of Orbital Propagation
    play(sc_main);

else
    fprintf("3D VIEWER DISABLED - ENABLE IN SCENARIO DEFINITION SECTION\n")
end

%% Position/Time/Velocity Data Processing

fprintf("Processing Data from Orbit Propagation...")

% Tx Location/Time History

% Rx Location/Time History - MAIN PLANE
[pos_Rx_1_1, vel_Rx_1_1, time_vec] = states(Rx_1(1));
[pos_Rx_1_2, vel_Rx_1_2] = states(Rx_1(2));
[pos_Rx_1_3, vel_Rx_1_3] = states(Rx_1(3));
[pos_Rx_1_4, vel_Rx_1_4] = states(Rx_1(4));
[pos_Rx_1_5, vel_Rx_1_5] = states(Rx_1(5));
[pos_Rx_1_6, vel_Rx_1_6] = states(Rx_1(6));
[pos_Rx_1_7, vel_Rx_1_7] = states(Rx_1(7));
[pos_Rx_1_8, vel_Rx_1_8] = states(Rx_1(8));
[pos_Rx_1_9, vel_Rx_1_9] = states(Rx_1(9));
[pos_Rx_1_10, vel_Rx_1_10] = states(Rx_1(10));
[pos_Rx_1_11, vel_Rx_1_11] = states(Rx_1(11));
[pos_Rx_1_12, vel_Rx_1_12] = states(Rx_1(12));
[pos_Rx_1_13, vel_Rx_1_13] = states(Rx_1(13));
[pos_Rx_1_14, vel_Rx_1_14] = states(Rx_1(14));
[pos_Rx_1_15, vel_Rx_1_15] = states(Rx_1(15));
[pos_Rx_1_16, vel_Rx_1_16] = states(Rx_1(16));
[pos_Rx_1_17, vel_Rx_1_17] = states(Rx_1(17));
[pos_Rx_1_18, vel_Rx_1_18] = states(Rx_1(18));
[pos_Rx_1_19, vel_Rx_1_19] = states(Rx_1(19));
[pos_Rx_1_20, vel_Rx_1_20] = states(Rx_1(20));
[pos_Rx_1_21, vel_Rx_1_21] = states(Rx_1(21));
[pos_Rx_1_22, vel_Rx_1_22] = states(Rx_1(22));
[pos_Rx_1_23, vel_Rx_1_23] = states(Rx_1(23));
[pos_Rx_1_24, vel_Rx_1_24] = states(Rx_1(24));


% Rx Location/Time History - SECONDARY PLANE
[pos_Rx_2_1, vel_Rx_2_1] = states(Rx_2(1));
[pos_Rx_2_2, vel_Rx_2_2] = states(Rx_2(2));
[pos_Rx_2_3, vel_Rx_2_3] = states(Rx_2(3));
[pos_Rx_2_4, vel_Rx_2_4] = states(Rx_2(4));
[pos_Rx_2_5, vel_Rx_2_5] = states(Rx_2(5));
[pos_Rx_2_6, vel_Rx_2_6] = states(Rx_2(6));
[pos_Rx_2_7, vel_Rx_2_7] = states(Rx_2(7));
[pos_Rx_2_8, vel_Rx_2_8] = states(Rx_2(8));
[pos_Rx_2_9, vel_Rx_2_9] = states(Rx_2(9));
[pos_Rx_2_10, vel_Rx_2_10] = states(Rx_2(10));
[pos_Rx_2_11, vel_Rx_2_11] = states(Rx_2(11));
[pos_Rx_2_12, vel_Rx_2_12] = states(Rx_2(12));

% Tx/Rx Relative Position History

fprintf("Done!\n")

%% Converting to LLA

poolobj = parpool(4); % Parallel Pool w/ 4 Threads
fprintf("Converting position data from ECI to LLA Format...\n")

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
h_3 = worldmap("Ecuador");
getm(h_3,"MapProjection");
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
geoshow('worldlakes.shp', 'FaceColor', 'cyan')
geoshow('worldrivers.shp', 'Color', 'blue')
geoshow('worldcities.shp', 'Marker', '.',...
                           'MarkerEdgeColor', 'magenta')
title("100x100 km near Equator")

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

scaleruler on
scaleruler("units", "km")

fprintf("Done!\n")

%% Coverage of 100km x 100km grid at Latitude ~ 70 deg.

fprintf("Plotting 100x100km @ Lat. ~ 70 deg. Coverage...")
figure(6)
h_4 = worldmap("Iceland");
getm(h_4,"MapProjection");
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
geoshow('worldlakes.shp', 'FaceColor', 'cyan')
geoshow('worldrivers.shp', 'Color', 'blue')
geoshow('worldcities.shp', 'Marker', '.',...
                           'MarkerEdgeColor', 'magenta')
title("100x100 km near Lat = 70 deg.")

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

scaleruler on
scaleruler("units", "km")

fprintf("Done!\n")

toc % Display total runtime of NumericalSolver_2OA
