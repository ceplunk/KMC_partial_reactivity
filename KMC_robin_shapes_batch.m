%% KMC_robin_shapes_batch.m
%% Script to run a batch of the simulations in KMC_robin_shapes.m
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

rng shuffle

filename_end = '_0523';

parpool

loading = load('cyl_cdf.mat');
cyl_cdf = loading.cyl_cdf;

Rpatch1 = 1;
Rmax = 1e18;
tic

totTrials = 1e2;
kappaRange = logspace(-5,5,50);
RpatchRange = 1;
shapeNum = 1;

capacitances = zeros(length(kappaRange),length(RpatchRange));

for ii = 1:length(kappaRange)
    for jj = 1:length(RpatchRange)
        capacitances(ii,jj) = KMC_robin_shapes(shapeNum,kappaRange(ii), ...
            Rpatch1, RpatchRange(jj),  Rmax, totTrials, cyl_cdf);
    end
end

fprintf('Saving capacitances...\n')
save(strcat('capacitances_inner',filename_end,'.mat'),'capacitances', '-v7.3');

fprintf('Done.\n')
toc