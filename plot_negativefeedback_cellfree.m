function plot_negativefeedback_cellfree

% This is to plot the data for negative feedback for cell-free context

close all 
clear all

% Load the data
% load negative_feedback_cellfree_resource_limited

% uncomment to plot parameter distributions, transient response

% parameter_distributions(o);

% transient_response(o);


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
    y = o.default(i).transient_response;
    y = y/y(end);
%     figure(i);
%     plot(t, y/y(end), 'b');
%     hold on;
    for k = 1:length(x_bins)-1
        time_bin(k).indices = find(y > y_bins(k) & y < y_bins(k+1));
    end
    time_bin(length(x_bins)).indices =  find(y > y_bins(length(x_bins)-1));
    
    
    for j=1:N
        y_r = o.default(i).r_transient(j,:);
        y_r = y_r/y_r(end);
%         figure(i);
%     plot(t, y_r/y_r(end), 'r');
%     hold on;
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
function parameter_distributions(o)

% plot steady state, 1x1
M = 100; % number of random samplings of parameter space
N = 100; % number of random temperture changes about chosen parameter point

ratio_alpha = [];
ratio_gamma = [];
ratio_K = [];
for i = 1:M
    ratio_alpha = [ratio_alpha, o.default(i).r_alpha];
    ratio_gamma = [ratio_gamma, o.default(i).r_gamma];
    ratio_K = [ratio_K, o.default(i).r_K];    
end


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
