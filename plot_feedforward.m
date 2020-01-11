function plot_feedforward

% This is to plot all the data of loaded feedforward circuits

close all 
clear all

% Load required data
% load feedforward_data
% load feedforward_data_gamma10
% load feedforward_data_gammap1

% Choose steady-state or transient to plot either

steady_state(o);

% transient_response(o);


%%
function transient_response(o)

M = 100;
N = 100;

t = o.time;

peak_ratio = [];
peak_time_ratio = [];
decay_time_ratio = [];
for i = 1:M
%     figure(1);
%     plot(t, o.default(i).transient_response,'r');
%     hold on;

    % This determines the maximum value of the pulse and the time that it
    % is reached
    [default_peak, default_peak_index] = max(o.default(i).transient_response);
    
    % This determines the decay time
    default_difference_from_final_value = abs(o.default(i).transient_response(end) - o.default(i).transient_response(default_peak_index:end));
    default_decay_time_index = min(find(default_difference_from_final_value < 0.05)); %5% settling time criterion
    default_decay_time = t(default_peak_index + default_decay_time_index);
    for j=1:N
%         figure(1);
%         plot(t, o.default(i).r_transient(j,:),'b');
%         hold on;
    
        % This determines the maximum value of the pulse and the time that it
        % is reached
        [r_peak, r_peak_index] = max(o.default(i).r_transient(j,:));
        peak_ratio = [peak_ratio, r_peak/default_peak];
        peak_time_ratio = [peak_time_ratio, t(r_peak_index)/t(default_peak_index)];
        
        % This determines the decay time
        r_difference_from_final_value = abs(o.default(i).r_transient(j,end) - o.default(i).r_transient(j,r_peak_index:end));
        r_decay_time_index = min(find(r_difference_from_final_value < 0.05)); %5% settling time criterion
        r_decay_time = t(r_peak_index + r_decay_time_index);
        decay_time_ratio = [decay_time_ratio, r_decay_time/default_decay_time];
        
        clear r_peak* r_d*
    end
    clear default_*
end

figure;
[nelements, centers] = hist(peak_ratio, 0:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(1*ones(size(x)), x, 'k', 0.66^3*ones(size(x)), x, 'k--', 1.5^3*ones(size(x)), x, 'k--');
hold on;
axis([-.1 3.5 0 1.1]);
mean(peak_ratio)
std(peak_ratio)

figure;
[nelements, centers] = hist(peak_time_ratio, 0:.05:1.4);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(1*ones(size(x)), x, 'k', 0.33*ones(size(x)), x, 'k--', 0.5*ones(size(x)), x, 'k--');
hold on;
axis([-.1 1.5 0 1.1]);
mean(peak_time_ratio)
std(peak_time_ratio)


figure;
[nelements, centers] = hist(decay_time_ratio, 0:.05:1.4);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(1*ones(size(x)), x, 'k', 0.33*ones(size(x)), x, 'k--', 0.5*ones(size(x)), x, 'k--');
hold on;
axis([-.1 1.5 0 1.1]);
mean(decay_time_ratio)
std(decay_time_ratio)

%%
function steady_state(o)

% plot steady state, 1x1
M = 100; % number of random samplings of parameter space
% N = 100; % number of random temperture changes about chosen parameter point

ratio = [];
ratio_alpha_x = [];
ratio_gamma_x = [];
ratio_alpha_y = [];
ratio_gamma_y = [];
ratio_K_x = [];
for i = 1:M
    ratio = [ratio, o.default(i).r_steady_state/o.default(i).steady_state];
    ratio_alpha_x = [ratio_alpha_x, o.default(i).r_alpha_x];
    ratio_gamma_x = [ratio_gamma_x, o.default(i).r_gamma_x];
    ratio_alpha_y = [ratio_alpha_y, o.default(i).r_alpha_y];
    ratio_gamma_y = [ratio_gamma_y, o.default(i).r_gamma_y];
    ratio_K_x = [ratio_K_x, o.default(i).r_K_x];
end

mean(ratio)
std(ratio)
figure;
[nelements, centers] = hist(ratio, 0:.1:3);
bar(centers, nelements/max(nelements), 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'BarWidth', 1);
hold on;
x = linspace(0,1);
plot(1*ones(size(x)), x, 'k', 0.66^3*ones(size(x)), x, 'k--', 1.5^3*ones(size(x)), x, 'k--');
hold on;
axis([-.1 3.5 0 1.1]);

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


