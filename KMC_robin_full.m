%% KMC_robin_full.m
%% Function for kinetic Monte Carlo simulations of survival probability
%  for a single particle with isotropic diffusion above a reflective plane
%  with infinite Robin disks of radius eps and reactivity kappa/eps on a
%  unit square grid
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

% we are looking at the PDE
% \partial_t S = \Delta S,    z > 0, t > 0
% \partial_z S = kappa S,     z = 0, (x,y) in patches of radius eps
% \partial_z S = 0,           z = 0, (x,y) not in patches of radius eps

function [finalTime, successes] = KMC_robin_full(kappa,eps,z0,totTrials,cyl_cdf,tau_PS_B)

reactivity = kappa/eps;
L = 1;

% numerical inversion is based on the following data
tau_PS_B = tau_PS_B(33:end,:);

successes = 0;
finalTime = zeros(1,totTrials);

parfor trial = 1:totTrials
    timeElapsed = 0;
    % initialize location of particles
    x = rand;
    y = rand;
    z = z0;

    % after initializing, simulate until success or wanders too far away
    while 1<2

        % simulate to z=0 plane
        % use a uniform random variable with the correct transformation
        % to simulate how long it took to reach the plane
        t1 = 0.25*(z/erfcinv(rand))^2;
        timeElapsed = timeElapsed + t1;
        % using that time, use normal random variables to simulate new 
        % location
        x = x + sqrt(2*t1)*randn;
        y = y + sqrt(2*t1)*randn;
        z = 0;

        % calculate the distance to BC patch
        % bound: positive or negative to determine if outside or inside
        % dist: actually calculate distance; always positive
        dist = sqrt((mod(x,1/sqrt(L))-0.5/sqrt(L))^2 + ...
            (mod(y,1/sqrt(L))-0.5/sqrt(L))^2) - eps;
        
        if dist <= 0
            % calculate time tau for (x,y) to diffuse to x^2 + y^2 > dist^2
            tau = cyl_time_get(-dist,cyl_cdf,1);
            
            tau_test = tau > 0;
            if ~tau_test
                [1 1]*[1 1 1]
            end

            % for all time 0 < t < tau, z sees a kappa BC
            % can just use equation 4.4 to determine if absorbed or where
            % it goes next

            [z,t3] = absorb_time_or_loc(tau, reactivity);

            if isnan(z)
                successes = successes + 1;
                timeElapsed = timeElapsed + t3;
                finalTime(trial) = timeElapsed;
                break

            % distribute (x,y) on x^2 + y^2 = dist^2
            else
                timeElapsed = timeElapsed + tau;
                angle = 2*pi*rand;
                x = x + dist*cos(angle);
                y = y + dist*sin(angle);
            end

        else
            % if BC is not met, determine where the particle goes next
            locvec = randn(1,3);
            const=dist/sqrt(locvec(1)^2 + locvec(2)^2 + locvec(3)^2);
            z = const*abs(locvec(1));
            x = const*locvec(2)+x;
            y = const*locvec(3)+y;
            t2 = bulk_time_get(dist, tau_PS_B, 1);
            timeElapsed = timeElapsed + t2;

        end

    end

end

