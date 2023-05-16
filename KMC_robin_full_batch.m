%% KMC_local_patch_batch.m
%% Script to run a batch of the simulations in KMC_robin_full.m
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

rng shuffle

filename_end = '_0523';

parpool(32)

loading = load('tau_PS_B.mat');
tau_PS_B = loading.tau_PS_B;
loading = load('cyl_cdf.mat');
cyl_cdf = loading.cyl_cdf;
tic

totTrials = 1e5;
eps = logspace(-2,-1,5);
kappa = [1e-2 1e-1];
z0 = 1;

finalTimes = zeros(length(eps),length(kappa),totTrials);
kappaOutList = zeros(length(eps),length(kappa));

for ii = 1:length(eps)
    for jj = 1:length(kappa)
        finalTimes(ii,jj,:) = KMC_robin_full(kappa(jj),eps(ii),z0,totTrials,cyl_cdf,tau_PS_B);
        save(strcat('finalTimes', filename_end, '.mat'),'finalTimes', '-v7.3');
        string = strcat('Finished with eps = ', num2str(eps(ii)), 'kappa = ', num2str(kappa(jj)), '\n');
        fprintf(string)
    end
end


fprintf('Finished saving times. Now optimizing kappa values.\n')

for ii = 1:length(eps)
    for jj = 1:length(kappa)
        c0 = 2*kappa(jj)/(pi*(kappa(jj)+1)); % c0 is just a function of kappa
        f = @(kappaOut)ks_distance(kappaOut,sort(reshape(finalTimes(ii,jj,:),1,[])),z0,1);
        kappaTh = 2*pi*eps(ii)*c0;
        kappaOutList(ii,jj) = fminsearch(f,kappaTh);
    end
end

fprintf('Saving kappas...\n')
save(strcat('kappaOutList',filename_end,'.mat'),'kappaOutList', '-v7.3');

fprintf('Done.\n')

toc

