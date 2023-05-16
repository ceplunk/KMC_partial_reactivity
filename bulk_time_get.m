%% bulk_time_get.m
%% Function to sample time for a particle diffusing in 3D to diffuse to a
%  sphere of radius d
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

% relies on tau_PS_B.mat from calculating_PS_bernoff.m
% this is the cdf from Bernoff et al 2018, equations 3.8 and 3.10, using
% linear interpolation and the relationship t = d^2*tau/Dt/pi^2

function [tnew] = bulk_time_get(d, tau_PS_B, Dt)
randT = rand;

for tauIndex=1:115000
    if tau_PS_B(tauIndex,2)>=randT, break, end
end

if tauIndex == 1
    tauU = tau_PS_B(tauIndex,1);
    tauL = 0;
    PSU = tau_PS_B(tauIndex,2);
    PSL = 0;
    fprintf('need smaller t \n')
else
    tauU = tau_PS_B(tauIndex,1);
    tauL = tau_PS_B(tauIndex-1,1);
    PSU = tau_PS_B(tauIndex,2);
    PSL = tau_PS_B(tauIndex-1,2);
end
tau = (tauL*(PSU-randT) + tauU*(randT - PSL))/( PSU - PSL );
tnew = d^2*tau/Dt/pi^2;