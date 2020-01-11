function plot_negativefeedback

% This is to plot the data for negative feedback for the loaded data

close all 
clear all

% load negative_feedback
% load negative_feedback_gamma10
% load negative_feedback_gammap1

% Uncomment to plot steady_state, transient response, identify parameters,
% or plot the negative feedback regimes

% steady_state(o);

% transient_response(o);

% identify_parameters(o);

% negativefeedbackregime(o);

%%
function negativefeedbackregime(o)
%%
X = logspace(-1,5);
K = o.K/10;
y1 = o.alpha./(1+(X/K));
y2 = o.gamma*X;
figure;
loglog(X,y1,'k', X, y2, 'k', 'LineWidth', 4);
axis([1e-1 1e5 1e-1 1e3]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

K = o.K*100;
y1 = o.alpha./(1+(X/K));
y2 = o.gamma*X;
figure;
loglog(X,y1,'k', X, y2, 'k', 'LineWidth', 4);
axis([1e-1 1e5 1e-1 1e3]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
function identify_parameters(o)

M = 100;
N = 100;

t = o.time;

alpha = o.alpha*10.^o.set(1,:);
gamma = o.gamma*10.^o.set(2,:);
K = o.K*10.^o.set(3,:);
% three parameters, so two ratios
p1 = K;
p2 = alpha./(gamma);

distribution = []; % for those that are ratio(i) > 0.5

for i=1:M
    
    y = o.default(i).transient_response/o.default(i).steady_state;
    time_indices = find(y > 0.4 & y < 0.6);

    slowerthantwofold(i) = 0;
    fasterthantwofold(i) = 0;

    for j=1:N
        y_r = o.default(i).r_transient(j,:);
        
        region = y_r(time_indices)' < 2*y(time_indices) - y(time_indices).*y(time_indices);
        slowerthantwofold(i) = slowerthantwofold(i) + length(find(region == 1));
        fasterthantwofold(i) = fasterthantwofold(i) + length(find(region == 0));
        
    end
    ratio(i) = slowerthantwofold(i)/(slowerthantwofold(i) + fasterthantwofold(i));
    xss(i) = o.default(i).steady_state;
    meanQ(i) = mean(o.default(i).r_steady_state/o.default(i).steady_state);
    stdQ(i) = std(o.default(i).r_steady_state/o.default(i).steady_state);

    figure(1000);
    loglog(p1(i), p2(i), 'ko', 'MarkerFaceColor', (1-ratio(i))*[1 1 1], 'MarkerSize', 15);
    hold on;


end

figure(1000);
x = logspace(0, 2);
loglog(x, x, 'k');

figure;
semilogx(p2./p1, meanQ, 'ks', 'MarkerFaceColor', [1 1 1], 'MarkerSize', 6);
hold on;
semilogx(p2./p1, stdQ, 'kd', 'MarkerFaceColor', [1 1 1], 'MarkerSize', 6);
hold on;
semilogx(p2./p1, xss./p2, 'ko', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 6);
hold on;

% last curve can be calculated analytically
% xss = (-K + sqrt(K^2 + 4*K*alpha/gamma))/2;
% p1 = K;
% p2 = alpha/(gamma);
% xss/p2 = (-(p1/p2) + sqrt((p1/p2)^2 + 4*p1/p2))/2

xx = logspace(-2,4);
yy = (-1./xx + sqrt((1./xx).^2 + 4./xx))/2;

semilogx(xx, yy,'k');
axis([1e-2 1e4 0 1.2]);


%%
function transient_response(o)

M = 100;
N = 100;

t = o.time;

x_bins = 0:0.02:1;
y_bins = 0:0.02:1;
y_bins_count = 0*ones(length(x_bins), length(y_bins));


for i=1:M
    
% find the time bins corresponding to thedefault trajectory for this parameter set
    y = o.default(i).transient_response/o.default(i).steady_state;
    for k = 1:length(x_bins)-1
        time_bin(k).indices = find(y > y_bins(k) & y < y_bins(k+1));
    end
    time_bin(length(x_bins)).indices =  find(y > y_bins(length(x_bins)-1));
    
    
    for j=1:N
        y_r = o.default(i).r_transient(j,:);
%         for each trajectory, populate the bin count
        for k = 1:length(x_bins)
            n = hist(y_r(time_bin(k).indices), y_bins);
            binary_n = n>0;
            y_bins_count(k,:) = y_bins_count(k,:) + binary_n;
        end
        

    end
end



figure(1);
imagesc((1-y_bins_count'/N*M));
colormap('gray');
set(gca,'YDir','normal'); % to keep the y-axis from being inverted
hold on;
% the grid
t = 0:0.01:10;
exponent = [0.5, 1, 2, 3];
for i = 1:length(exponent)
    figure(1);
    plot(length(x_bins)*(1-exp(-t)), length(y_bins)*(1-exp(-exponent(i)*(t))),'k');
    hold on;
end
axis([0 1.01*length(x_bins) 0 1.01*length(y_bins)]);

%%
function steady_state(o)

% plot steady state, 1x1
M = 100; % number of random samplings of parameter space
N = 100; % number of random temperture changes about chosen parameter point

ratio = [];
ratio_alpha = [];
ratio_gamma = [];
ratio_K = [];
for i = 1:M
    ratio = [ratio, o.default(i).r_steady_state/o.default(i).steady_state];
    ratio_alpha = [ratio_alpha, o.default(i).r_alpha];
    ratio_gamma = [ratio_gamma, o.default(i).r_gamma];
    ratio_K = [ratio_K, o.default(i).r_K];    
end

mean(ratio)
std(ratio)
figure;
[nelements, centers] = hist(ratio, 0:.1:2);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(1*ones(size(x)), x, 'k', 0.66*ones(size(x)), x, 'k--', 1.5*ones(size(x)), x, 'k--');
hold on;
axis([-.1 2.1 0 1.1]);

figure;
[nelements, centers] = hist(ratio_alpha, 2:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(2*ones(size(x)), x, 'k--', 3*ones(size(x)), x, 'k--');
hold on;
axis([1.9 3.1 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

figure;
[nelements, centers] = hist(ratio_gamma, 2:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(2*ones(size(x)), x, 'k--', 3*ones(size(x)), x, 'k--');
hold on;
axis([1.9 3.1 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


figure;
[nelements, centers] = hist(ratio_K, 0.5:.1:1.4);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(0.66*ones(size(x)), x, 'k--', 1.5*ones(size(x)), x, 'k--');
hold on;
axis([0.4 1.6 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
