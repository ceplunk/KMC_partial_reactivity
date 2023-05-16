%% calculating_PS_bernoff.m
%% Function to calculate the cdf for time to exit a sphere of radius 1
%  Following Bernoff et al 2018, equations 3.8 and 3.10
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

tau1 = [1:0.0001:10 10.001:0.001:25];
PS1 = 1 - 2.*exp(-tau1) + 2.*exp(-4*tau1) - 2.*exp(-9*tau1) ...
    + 2*exp(-16*tau1) - 2*exp(-25*tau1);

tau2 = 0.0001:0.0001:0.9999;
PS2 = 2*sqrt(pi./tau2).*( exp(-pi^2*(1/2)^2./tau2) ...
    + exp(-pi^2*(3/2)^2./tau2) + exp(-pi^2*(5/2)^2./tau2) ...
    + exp(-pi^2*(7/2)^2./tau2) + exp(-pi^2*(9/2)^2./tau2) );

tau = [tau2 tau1];
PS = [PS2 PS1];
tau_PS_B = [tau' PS'];
save('tau_PS_B.mat', 'tau_PS_B')