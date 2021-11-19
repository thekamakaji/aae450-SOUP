% Specular Point Calculator
% Author: Judy Park
% Purpose: To calculate the latitude and longitude of the spectator point
% on spherical Earth given the latitudes, longitudes, and altitudes of the
% receiver and transmitter satellites
% 
% To Be Provided By Code:
% - Visibility of the spectator point from satellites
% - Latitude, Longitude, and Altitude of spectator point on spherical Earth

function [visible, latSpec, longSpec] = specularPoint(latR, longR, hR, latT, longT, hT)

thetaR = (90 - abs(latR)) * sign(latR);
thetaT = (90 - abs(latT)) * sign(latT);

rE = 6378.1363; % radius of Earth [km]
dist = sqrt(2) * rE * sqrt(1 - sind(thetaR).*sind(thetaT).*cosd(longR-longT) ...
    - cosd(thetaR).*cosd(thetaT));
angleRT = acosd(1 - dist.^2/(2*rE^2));

% Solve for distances and angles using Trigonometry
syms dR dT incidence x
eqn1 = incidence == acosd((rE^2 + dR^2 - (rE + hR)^2)/(-2*dR*rE));
eqn2 = incidence == acosd((rE^2 + dT^2 - (rE + hT)^2)/(-2*dT*rE));
eqn3 = x == acosd((rE^2 + (rE + hR)^2 - dR^2)/(2*rE*(rE + hR)));
eqn4 = x == angleRT - acosd((rE^2 + (rE + hT)^2 - dT^2)/(2*rE*(rE + hT)));
eqns = [eqn1 eqn2 eqn3 eqn4];
S = vpasolve(eqns, [dR dT incidence x]);
incAngle = S.incidence;
specAngle = S.x;

% Find out whether the satellites can see the spectator point or not and
% find the latitude and longitude of the spectator point.
% Assume that both latitudes of receiver and transmitter satellites cannot
% be +-90 degrees 
if incAngle <= 60
    
    if abs(latR) == 90 && abs(latT) ~= 90
        longR = longT;
        [latSpec, longSpec, visible] = calcSpec(latR, longR, latT, longT, specAngle, angleRT);  
    elseif abs(latT) == 90 && abs(latR) ~= 90
        longT = longR;
        [latSpec, longSpec, visible] = calcSpec(latR, longR, latT, longT, specAngle, angleRT);
    elseif abs(latR) == 90 && abs(latT) == 90 && latR * latT > 0
        fprintf('Error: The satellites are at the same location.');
        visible = 0;
    else
    [latSpec, longSpec, visible] = calcSpec(latR, longR, latT, longT, specAngle, angleRT);
    end
else
    fprintf('The satellites cannot see the spectator point.');
    visible = 0;
end
end