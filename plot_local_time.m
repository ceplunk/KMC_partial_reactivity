%% plot_local_time.m
%% Script to make plots of results of calculating the capacitance of a
%  patch with fixed flux -1
%  Claire Plunkett and Sean Lawley
%  University of Utah
%  May 2023

loading = load('K1_0213_more.mat');
K1 = loading.K1;
loading = load('mean_local_times_0213_more.mat');
mtimes = loading.mtimes;

rho = [1.1 1.5:0.5:5];

f = figure(1);
f.Position = [0 0 900 400];
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(rho,mtimes,'o')
hold on
plot(rho,K1./rho)
hold off
set(gca, 'FontName', 'Century Schoolbook','FontSize',12)
xlabel('Starting Radius')
ylabel('Expected Value of Local Time')
legend('KMC Results','Estimated Value')
ttl = title('A','FontSize',20);
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.14;
ttl.Position(2) = 0.99;
ttl.HorizontalAlignment = 'left';  

nexttile
plot(rho, abs(mtimes-K1./rho)./mtimes,'o')
set(gca, 'FontName', 'Century Schoolbook','FontSize',12)
xlabel('Starting Radius')
ylabel('Relative Error')
ttl = title('B','FontSize',20);
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.14;
ttl.Position(2) = 0.99;
ttl.HorizontalAlignment = 'left';  
