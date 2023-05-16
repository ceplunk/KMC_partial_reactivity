%% absorb_or_loc.m
%% Function to sample if absorbed or exit location after tau time
%  for 1D diffusion above a partially reactive boundary with fixed
%  reactivity
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

function zOut = absorb_or_loc(tau, reactivity)
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

if isempty(idx) == 1
    zOut = NaN;
elseif idx == 1
    zOut = myRand/pSpace(idx)*zrange(idx);
else
    zOut = (pSpace(idx) - myRand) / (pSpace(idx) - pSpace(idx-1)) * zrange(idx-1) ...
        + (myRand - pSpace(idx - 1)) / (pSpace(idx) - pSpace(idx-1)) * zrange(idx);
end