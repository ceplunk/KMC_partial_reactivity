%% absorb_time_or_loc.m
%% Function to sample time of absorption or exit location after tau time
%  for 1D diffusion above a partially reactive boundary with fixed
%  reactivity
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

function [zOut,time] = absorb_time_or_loc(tau, reactivity)
% maximum tau will be 4.3438... , so zrange = linspace(0,10) is sufficient
zlen = 250;
zrange = logspace(-10,1,zlen);

pSpace = erfcx(reactivity*sqrt(tau)) - erfcx( (2*reactivity*tau + zrange) / ...
    (2 * sqrt(tau) ) ) .* exp( - zrange.^2 / (4*tau));

myRand = rand;
idx = [];
for testIndex=1:zlen
    if pSpace(testIndex)>myRand
        idx = testIndex;
        break
    end
end

if isempty(idx) == 0
    if idx == 1
        zOut = myRand/pSpace(idx)*zrange(idx);
    else
        zOut = (pSpace(idx) - myRand) / (pSpace(idx) - pSpace(idx-1)) * zrange(idx-1) ...
            + (myRand - pSpace(idx - 1)) / (pSpace(idx) - pSpace(idx-1)) * zrange(idx);
    end
    time = NaN;
else
    zOut = NaN;
    newRand = rand;
    trange = logspace(-20,log10(tau));
    if reactivity > 1e3
        trange = logspace(-30,log10(tau));
    end
    calcul = 1 - erfcx( (2*reactivity.*trange)./sqrt(4*trange));
    calcul = calcul/calcul(end);
    myidx = 50;
    for testIndex=1:length(trange)
        if calcul(testIndex) >= newRand
            myidx = testIndex;
            break
        end
    end
    if myidx == 1
        fprintf(num2str(myidx))
        time = trange(myidx)*newRand/calcul(myidx);
    else
        time = (trange(myidx-1)*(calcul(myidx)-newRand) + ...
            trange(myidx)*(newRand-calcul(myidx-1)))/...
            (calcul(myidx)-calcul(myidx-1));
    end
end