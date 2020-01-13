function plot_coherent_feedforwardloop

% This is to plot all the data

close all 
clear all

% Load data
load coherent_ffl

% Plot steady-state or delay ratio
steady_state(o);

delay_fun(o);

%%
function steady_state(o)

% plot steady state, 1x1
M = 100; % number of random samplings of parameter space
N = 100; % numer of random temperture changes about chosen parameter point

ratio = [];
ratio_alpha_x = [];
ratio_gamma_x = [];
ratio_alpha_y = [];
ratio_gamma_y = [];
ratio_K_x = [];
for i = 1:M
    ratio = [ratio, o.default(i).r_cffl_ss_f/o.default(i).cffl_ss_f];
    ratio_alpha_x = [ratio_alpha_x, o.default(i).Q_alphaX];
    ratio_gamma_x = [ratio_gamma_x, o.default(i).Q_gammaX];
    ratio_alpha_y = [ratio_alpha_y, o.default(i).Q_alphaY];
    ratio_gamma_y = [ratio_gamma_y, o.default(i).Q_gammaY];
    ratio_K_x = [ratio_K_x, o.default(i).Q_K];
end

% mean(ratio)
% std(ratio)
% figure;
% [nelements, centers] = hist(ratio, 0:.1:3);
% bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
% hold on;
% x = linspace(0,1);
% plot(1*ones(size(x)), x, 'k', 0.66*ones(size(x)), x, 'k--', 1.5*ones(size(x)), x, 'k--');
% hold on;
% axis([-.1 2.5 0 1.1]);


figure;
[nelements, centers] = hist(ratio_alpha_x, 2:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(2*ones(size(x)), x, 'k--', 3*ones(size(x)), x, 'k--');
hold on;
axis([1.9 3.1 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

figure;
[nelements, centers] = hist(ratio_gamma_x, 2:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(2*ones(size(x)), x, 'k--', 3*ones(size(x)), x, 'k--');
hold on;
axis([1.9 3.1 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


figure;
[nelements, centers] = hist(ratio_alpha_y, 2:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(2*ones(size(x)), x, 'k--', 3*ones(size(x)), x, 'k--');
hold on;
axis([1.9 3.1 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

figure;
[nelements, centers] = hist(ratio_gamma_y, 2:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(2*ones(size(x)), x, 'k--', 3*ones(size(x)), x, 'k--');
hold on;
axis([1.9 3.1 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


figure;
[nelements, centers] = hist(ratio_K_x, 0.5:.1:1.4);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(0.66*ones(size(x)), x, 'k--', 1.5*ones(size(x)), x, 'k--');
hold on;
axis([0.4 1.6 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

%%
function delay_fun(o)

% plot steady state, 1x1
M = 100; % number of random samplings of parameter space
N = 100; % numer of random temperture changes about chosen parameter point

ratio = [];
for i = 1:M
    ratio = [ratio, o.default(i).r_delay/o.default(i).delay];
end

mean(ratio)
std(ratio)
figure;
[nelements, centers] = hist(ratio, 0:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(1*ones(size(x)), x, 'k');
x = linspace(0,1);
plot(1*ones(size(x)), x, 'k', 0.5*ones(size(x)), x, 'k--', .33*ones(size(x)), x, 'k--');
hold on;
hold on;
 axis([-.1 1.5 0 1.1]);

