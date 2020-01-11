function plot_nofeedback

% This is to plot the data for no feedback circuit for loaded data

close all 
clear all

% load no_feedback
% load no_feedback_gamma10
% load no_feedback_gammap1


% Uncomment for steady state or transient response

% steady_state(o);

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
axis([0 1.01*length(x_bins) 0 1.01*length(y_bins)])



%%
function steady_state(o)

% plot steady state, 1x1
M = 100; % number of random samplings of parameter space
N = 100; % number of random temperture changes about chosen parameter point

ratio = [];
ratio_alpha = [];
ratio_gamma = [];
for i = 1:M
    ratio = [ratio, o.default(i).r_steady_state/o.default(i).steady_state];
    ratio_alpha = [ratio_alpha, o.default(i).r_alpha];
    ratio_gamma = [ratio_gamma, o.default(i).r_gamma];
    
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
plot(1*ones(size(x)), x, 'k', 2*ones(size(x)), x, 'k--', 3*ones(size(x)), x, 'k--');
hold on;
axis([1.9 3.1 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

figure;
[nelements, centers] = hist(ratio_gamma, 2:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(1*ones(size(x)), x, 'k', 2*ones(size(x)), x, 'k--', 3*ones(size(x)), x, 'k--');
hold on;
axis([1.9 3.1 0 1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
