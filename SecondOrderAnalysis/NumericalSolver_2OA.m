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
p_days = 2; % Days of simulation propogation

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

% PRIMARY ORBITAL PLANE (24)

% SECONDARY ORBITAL PLANE (12)


%% Orbital Simulation


%% Data Processing/Export

% Tx Location/Time History

% Rx Location/Time History

% Tx/Rx Relative Position History


