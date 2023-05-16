%% KMC_local_patch_batch.m
%% Script to run a batch of the simulations in KMC_local_patch.m
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

rng shuffle

filename_end = '_0516';

parpool % change to parpool(32) for 32 cores, etc

loading = load('cyl_cdf.mat');
cyl_cdf = loading.cyl_cdf;

Rmax = 1e16;
tic

totTrials = 1e8;
R = 1;
Rinit = [1.1 2 3];

times = zeros(length(Rinit),totTrials);

for ii = 1:length(Rinit)
    times(ii,:) = KMC_local_patch(R,Rinit(ii),Rmax,totTrials,cyl_cdf);
end

fprintf('Saving times...\n')
save(strcat('times',filename_end,'.mat'),'times', '-v7.3');

fprintf('Calculating constants...\n')
fbar = mean(times,2);
K1 = mean(fbar.*Rinit');

fprintf(strcat('K1 = ', num2str(K1), '\n'))

save(strcat('K1',filename_end,'.mat'),'K1');

fprintf('Done.\n')
toc