%% KMC_local_patch.m
%% Function for kinetic Monte Carlo simulations of local time accumulated
%  by a single circular patch with given radius, for isotropic diffusion
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

% Diffusion equation we are considering:
% \partial_t p = \Delta p,  z > 0, (x,y) \in \R^2, t > 0
% \partial_z p = 1,         z = 0, (x,y) in patch
% \partial_z p = 0,      	z = 0, (x,y) not in patch

function times = KMC_local_patch(R,Rinit,Rmax,totTrials,cyl_cdf)

% use a small initial radius Rinit that is outside the boundary conditions 
% for z = 0

times = nan(1,totTrials);

parfor trial = 1:totTrials
    % initialize location of particles
    % start with 3 iid standard normal random variables
    locvec = randn(1,3);
    % adjust by this constant to make distance to origin equal to Rinit
    const = Rinit/sqrt(locvec(1)^2 + locvec(2)^2 + locvec(3)^2);
    x = const*locvec(1);
    y = const*locvec(2);
    z = const*abs(locvec(3));
    totaldist = sqrt(x^2+y^2+z^2);
    localtime = 0;

    % after initializing, simulate until success or wanders too far away
    while 1 < 2

        % simulate to z=0 plane
        % use a uniform random variable with the correct transformation
        % to simulate how long it took to reach the plane
        tt = 0.25*(z/erfcinv(rand))^2;
        % using that time, use normal random variables to simulate new 
        % location
        x = x + sqrt(2*tt)*randn;
        y = y + sqrt(2*tt)*randn;

        % calculate the distance to BC patch
        dist = sqrt( x^2 + y^2 ) - R;
        % check if BC is met
        if dist <= 0
            % calculate time tau for (x,y) to diffuse to x^2 + y^2 > dist^2
            % Deaconu et al. 2017
            tau = cyl_time_get(-dist,cyl_cdf,1)

            % for all time 0 < t < tau, calculate how much local time is
            % accumulated and where z exits
            % sample change in local time and then z location
            dlt =  2*sqrt(tau)*erfinv(rand);
            localtime = localtime + dlt;
            z = sqrt(dlt^2-4*tau*log(rand))-dlt;
            
            % distribute (x,y) on x^2 + y^2 = dist^2
            angle = 2*pi*rand;
            x = x + dist*cos(angle);
            y = y + dist*sin(angle);
            % check distance from origin. If > Rmax, exit
            totaldist = sqrt(x^2 + y^2 + z^2);

        else
            % check if too far away
            totaldist = sqrt(x^2 + y^2);
            if totaldist > Rmax
                times(trial) = localtime;
                break
            end

            % if BC is not met, determine where the particle goes next
            locvec = randn(1,3);
            const=dist/sqrt(locvec(1)^2 + locvec(2)^2 + locvec(3)^2);
            z = const*abs(locvec(1));
            x = const*locvec(2)+x;
            y = const*locvec(3)+y;

            % check distance from origin. If > Rmax, exit
            totaldist = sqrt(x^2 + y^2 + z^2);

        end

    end

end