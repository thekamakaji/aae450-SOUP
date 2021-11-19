%% Spectator Point Inclination Calculator

incT = 64.8; % Transmitter sat inclination [degrees]
incR = 85; % Receiver sat inclination [degrees]
rE = 6378; % Radius of Earth [km]
RecAlt = 750; % Altitude of Receiver Sat [km]
TrAlt = 25508; % Transmitter sat (GLONASS) altitude [km]

% Law of Cosine
syms a
aSol = solve(a^2 + rE * a + rE^2 - (rE + RecAlt)^2 == 0, a);
SpecToRec = double(max(aSol)); % distance between Spectator point and Receiver sat [km]
recAngle = acosd((rE^2 + (rE + RecAlt)^2 - SpecToRec^2)/(2*rE*(rE+RecAlt)));
% angle between spectator point and receiver sat [degrees]

SpecInc = incR - recAngle % Inclination of Spectator Point [degrees]