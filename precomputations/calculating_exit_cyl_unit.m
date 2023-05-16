%% calculating_exit_cyl_unit.m
%% Function to calculate the cdf for time to exit a cylinder of radius 1
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

% Theorem 2 from Ciesielski and Taylor 1962 gives tail distribution 
% Relies on the function besselzero.m

L = 1;
range = logspace(-2,1.1,200);

delta = 2;
nu = delta/2-1;
Linv = 1./L;
terms = 20000000; % 2e6
bzeros = besselzero(nu, terms); % first terms zeros of besselj(nu,k)
inputRL = range' * Linv;
tailDist = zeros(size(inputRL));

for k = 1:terms
    tailDist = tailDist + 1/(2^(nu-1)*gamma(nu+1)) * ...
        bzeros(k)^(nu-1) / besselj(nu+1,bzeros(k)) * exp(- bzeros(k)^2/2 * inputRL );
end

mycdf = 1 - tailDist;
cyl_cdf = abs([range' mycdf]);
save('cyl_cdf.mat', 'cyl_cdf')