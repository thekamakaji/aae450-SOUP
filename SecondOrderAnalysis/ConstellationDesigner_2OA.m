% 2nd Order Analysis - Constellation Designer
% Author: V. Swaminathan
% Version: 10/13/2021 0834 EST
% Purpose: To help develop Team SoUP's satellite constellation to meet
%          requirements of AAE450 Project (2nd order analysis)
% 
% To Be Provided By Code:
% - Maximum Revisit Times for Target Latitudes
% - Design spaces for inclination and altitude
% - Satellite constellation parameters

%% Initialize Workspace
clear all;
close all;

%% Constants
R_A = 6378.137; % Earth Equatorial Radius [km]
R_B = 6356.752; % Earth Polar Radius [km]
MU = 3.986e5; % Gravitational Param. for Earth [km^3/s^2]
J2 = 1082.63e-6; % J2 param. for Earth
SideDay_E = 86164; % Sidereal Day Length [s]
SolarDay_E = 86400; % Solar Day Length [s]
W_E = (-0.2507)/60; % Rate of Earth Rotation [deg./s]

%% Surface Target/Propogation Definitions

Lat = 30; % Target Latitude(s) for Observation [deg.] --> Positive only (assume mirrored performance)
startLong = 0; % Starting Longitude of Reciever in Ref. Plane [deg.]
timeTotal = 3*SolarDay_E; % Time to propogate through [s] --> Multiplier is Days (24 hrs.)

% Print Surface Target/Propogation Definitions
fprintf("-----------------------------------------------------------------\n");
fprintf("Target Latitude: %.2f degrees\n", Lat);
fprintf("Starting Longitude of Ref. Orbit: %.2f degrees\n", startLong);
fprintf("Propogation Time: %.2f days\n", timeTotal/SolarDay_E);

%% Satellite/Constellation Definitions

% Define Walker Constellation Parameters to Test
C_totalSats = 1; % Total Number of Rx satellites in Constellation
C_planes = 1; % Number of Equally spaced orbital planes
C_spacing = 0; % Spacing b/w satellites in adjacent planes (true anomaly "slots")

% Sensor/Instrumentation Parameters
S_phi = 60; % Half-cone boresight angle of Antennae [deg.] (Incd. Ang.)

% Orbital Parameters for Reference Receiver (Rx)
Rx_a = 750 + R_A; % Orbit radius [km] (semimajor axis (circular))
Rx_e = 0; % Eccentricity
Rx_i = 85; % Inclination [deg.]
Rx_RAAN = 0; % Right Ascension of AN [deg.]
Rx_w = 0; % Argument of Perigee [deg.]
Rx_p = Rx_a*(1 - Rx_e^2); % Semilatus Rectum p [km]
Rx_v = asind(sind(Lat)/sind(Rx_i)) - Rx_w; % True Anomaly @ target Latitude [deg.]
Rx_rs = Rx_p / (1 + (Rx_e*cosd(Rx_v))); % Orbital Radius @ Target Latitude [km]
Rx_hs = Rx_rs - geodeticR(R_A, R_B, Lat); % Orbital Altitude @ Target Latitude [km]

% Print out parameters from satellite definitions
fprintf("-----------------------------------------------------------------\n");
fprintf("Number of Satellites in Constellation: %.0f\n", C_totalSats);
fprintf("Number of Equally Spaced Orbital Planes: %.0f\n", C_planes);
fprintf("Relative Spacing of Sats. in Adjacent Planes: %.0f\n", C_spacing);
fprintf("-----------------------------------------------------------------\n");
fprintf("Ref. Orbit Altitude: %.2f km\n", Rx_a - R_A);
fprintf("Ref. Orbit Inclination: %.2f degrees\n", Rx_i);
fprintf("Ref. Orbit Eccentricity: %.2f\n", Rx_e);
fprintf("Ref. Orbit RAAN: %.2f degrees\n", Rx_RAAN);
fprintf("Ref. Orbit Arg. of Perigee: %.2f degrees\n", Rx_w);
fprintf("Orbital Radius @ Target Lat.: %.2f km\n", Rx_rs);
fprintf("Orbital Altitude @ Target Lat.: %.2f km\n", Rx_hs);

%% Calculation - Longitude of Successive Passes of Reference Receiver Sat.

% Basic Orbital Parameters
Rx_Pk = keplerianPeriod(Rx_a, MU); % Keplerian Orbit Period [s]
Rx_Pn = nodalPeriod(Rx_Pk, J2, R_A, Rx_p, Rx_e, Rx_i); % Nodal Period [s]
Rx_n = 360/Rx_Pk; % Mean motion [deg./s]

% Ground Track Drift Parameters
Rx_dANfromJ = regressionRateAN(Rx_n, J2, R_A, Rx_p, Rx_i, Rx_e); % Shift in AN caused by J2, J3, J4 zonal harmonics [deg./s]
Rx_dLong = Rx_Pn * (W_E + Rx_dANfromJ); % Shift in Surface Longitude per orbit [deg.]
Rx_surfLongPass1 = startLong + surfaceLong(Rx_w, Rx_v, Rx_RAAN, Rx_i, Rx_dLong); % Surface Longitude of 1st pass over target latitude [deg.]

j_n = (timeTotal / Rx_Pn); % Number of orbits within propogation time

% Surface Longitudes @ Consecutive Passes [deg.]
% !!! NEED TO ADD DESCENDING PASS COVERAGE !!!
Rx_surfLong(1) = Rx_surfLongPass1;
for x = 2:1:j_n
    Rx_surfLong(x) = Rx_surfLongPass1 + (x-1)*Rx_dLong;
end

%% Calculation - Satellite Constellation Passes Over Target Latitude
% !!! NEED TO ADD DESCENDING PASS COVERAGE !!!

C_satsPerPlane = C_totalSats/C_planes; % Number of Satellites in Each Orbital Plane

% Longitude Passes of First Satellite in Each of Multiple Planes @ Target Lat. [deg.]
C_surfLongFirst(1,:) = Rx_surfLong;
for m = 1:1:(C_planes - 1)
    C_surfLongFirst(m+1,:) = Rx_surfLong + 2*pi*m*((1/C_planes)+(C_spacing/C_totalSats));
end

% Add Longitude Passes of Multiple Satellites in Same Planes @ Target Lat. [deg.]
% !!! NEEDS TO BE CHECKED AND VERIFIED !!!
for x = 1:1:C_planes
    tempVec = C_surfLongFirst(x,:);
    for L = 1:1:(C_satsPerPlane - 1)
        temp = C_surfLongFirst(x,:) + (L/C_satsPerPlane)*Rx_dLong;
        tempVec = [tempVec, temp];
    end
    C_surfLongAll(x,:) = tempVec;
end
% !!! VERIFY ABOVE CODE SECTION !!!

% Total Set of Passes by Longitude in Given Propogation Time
C_allPasses(1,:) = C_surfLongAll(1,:);
for x = 2:1:C_planes
    C_allPasses = [C_allPasses C_surfLongAll(x,:)];
end

% Normalize Passes to 360 deg. Earth
C_allPasses = wrapTo180(C_allPasses);
C_allPasses = unique(C_allPasses);

% Print out "gaps" in Longitude passes (not accounting for sensor capabilities
fprintf("-----------------------------------------------------------------\n");
fprintf("Number of Orbits by Each Sat. within Propogation Time: %.0f\n", j_n);
fprintf("Largest Gap in Longitude of Passes: %.3f degrees\n", max(diff(C_allPasses)));
fprintf("Largest Gap in Longitude of Passes: %.3f km\n", (max(diff(C_allPasses))/360)*(2*pi*R_A));

%% Calculation - Determination of Visible Longitudes

% !!! ISSUES FOR FOLLOWING LAT. VALUES !!! --> DEBUG NEEDED
% ~30 +- 2, ~60 +- 7, ~11, ~1 - ~8, ~20 - ~21, ~50

% Basic Sensor Parameters
% !!! DEBUG FOR LOWER LAT. VALUES !!!
R_geo = geodeticR(R_A, R_B, Lat); % Geodetic Radius @ Target Lat. [km]
HGRA = halfGRA(Rx_rs, S_phi, R_geo); % Half Ground Range Angle [deg.]
SDA = surfDA(HGRA, Lat); % Surface Dihedral Angle (longitude coverage of sensor) [deg.]

% Print out Longitude Coverage Range of Sensor
fprintf("-----------------------------------------------------------------\n");
fprintf("Coverage Range of Sensor in Surface Longitude: %.3f degrees\n", SDA);
fprintf("Coverage Range of Sensor in Surface Longitude: %.3f km\n", (SDA/360)*(2*pi*R_A));

% Longitude Coverage Range
% !!! IMPLEMENTATION NEEDED !!!

%% Plots of MRT/ART and Design Spaces
% !!! IMPLEMENTATION NEEDED !!!





