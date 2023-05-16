%% cyl_time_get.m
%% Function for calculating the time to diffuse to a cylinder of radius d
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

function [tnew] = cyl_time_get(d, cyl_cdf, Dt)
% relies on cyl_cdf.mat from calculating_exit_cyl_unit.m
% this is from Theorem 2 from Ciesielski and Taylor 1962, which gives the
% tail distribution, and then using linear interpolation and the
% relationship t = d^2*tau/Dt
if d == 0
    d = 1e-20;
end

randT = rand;

for tauIndex=1:200
    if cyl_cdf(tauIndex,2)>=randT, break, end
end


if tauIndex == 1 % this might break things
    tauU = cyl_cdf(tauIndex,1);
    tauL = 0;
    PSU = cyl_cdf(tauIndex,2);
    PSL = 0;
    fprintf('need smaller t \n')
else
    tauU = cyl_cdf(tauIndex,1);
    tauL = cyl_cdf(tauIndex-1,1);
    PSU = cyl_cdf(tauIndex,2);
    PSL = cyl_cdf(tauIndex-1,2);
end
tau = (tauL*(PSU-randT) + tauU*(randT - PSL))/( PSU - PSL );
tau_test = tau >0;
if ~tau_test
    tau
    tauIndex
    randT
    [1 1]*[1 1 1]
end
tnew = d^2*tau/Dt;