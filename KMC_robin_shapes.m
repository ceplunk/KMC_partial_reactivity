%% KMC_robin_shapes.m
%% Function for kinetic Monte Carlo simulations of splitting probability
%  for a single particle with isotropic diffusion above a reflective plane
%  with one Robin patch with a specified shape and reactivity kappa
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

% we are looking at the PDE:
% \partial_t p = \Delta p,  z > 0, (x,y) \in \R^2, t > 0
% \partial_z p = kappa p,	z = 0, (x,y) in patch
% \partial_z p = 0,      	z = 0, (x,y) not in patch

% uses Residuals_ellipse.m for ellipses

function capacitance = KMC_robin_shapes(shapeID,kappa,Rpatch1,Rpatch2,Rmax,totTrials,cyl_cdf)

% use a small initial radius that is outside the boundary conditions for 
% z = 0
Rinit = 1.01*max(Rpatch1,Rpatch2);

successes = 0;

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

    % after initializing, simulate until success or wanders too far away
    while totaldist < Rmax

        % simulate to z=0 plane
        % use a uniform random variable with the correct transformation
        % to simulate how long it took to reach the plane
        tt = 0.25*(z/erfcinv(rand))^2;
        % using that time, use normal random variables to simulate new 
        % location
        x = x + sqrt(2*tt)*randn;
        y = y + sqrt(2*tt)*randn;

        % calculate the distance to BC patch
        % bound: positive or negative to determine if outside or inside
        % dist: actually calculate distance; always positive
        if shapeID == 1 % circle
            bound = sqrt( x^2 + y^2 )- Rpatch1;
            dist = abs(bound);
            
        elseif shapeID == 2 % ellipse
            bound = Rpatch2^2*x^2/Rpatch1^2 + y^2 - Rpatch2^2;
            if bound <= 5
                dist = Residuals_ellipse( [x,y], [0,0,Rpatch1,Rpatch2,pi]);
                dist = 0.98*sqrt(dist);
                % a little extra safe since this algorithm has a tolerance 
                % of 1e-9 for speed
            else
                dist = sqrt(x^2 + y^2) - 1; % assuming can fit ellipse into
                % circle of radius 1
            end
        elseif shapeID == 3 % rectangle
            d1 = abs(x) - Rpatch1;
            d2 = abs(y) - Rpatch2;
            bound = max(d1,d2);
            dist = abs(bound);
        else
            fprintf('error: shapeID not matched')
        end
        % check if BC is met
        if bound <= 0
            % calculate time tau for (x,y) to diffuse to x^2 + y^2 > dist^2
            % Deaconu et al. 2017
            tau = cyl_time_get(dist,cyl_cdf,1);
            
            tau_test = tau > 0;
            if ~tau_test
                [1 1]*[1 1 1]
            end

            % for all time 0 < t < tau, z sees a kappa BC
            % can just use equation 4.4 to determine if absorbed or where
            % it goes next

            z = absorb_or_loc(tau, kappa);

            if isnan(z)
                successes = successes + 1;
                break

            % distribute (x,y) on x^2 + y^2 = dist^2
            else
                angle = 2*pi*rand;
                x = x + dist*cos(angle);
                y = y + dist*sin(angle);
                % check distance from origin. If > Rmax, exit
                totaldist = sqrt(x^2 + y^2 + z^2);
            end

        else

            % check if too far away
            totaldist = sqrt(x^2 + y^2);
            if totaldist > Rmax
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

% calculate capacitance by approximate probability of success * Rinit
capacitance = successes*Rinit/totTrials;