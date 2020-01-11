function nofeedback_gammap1

% Calculate temperature dependence using Q10 for no feedback circuit with
% degradation parameter gamma = 0.1/hr

close all 
clear all


% parameters default
o.alpha = 100; %nM/hr
o.gamma = 0.1; %/hr

% Change parameter sets
M = 100;
set = random('uniform', -1, 1, [2, M]);
time = 0:1e-1:100; % in hours

o.set = set;
o.time = time;

for j=1:M
    j
    tic
    % parameter regime point
    p.alpha = o.alpha*10^set(1,j); %nM/hr
    p.gamma = o.gamma*10^set(2,j); %/hr

    % calculate quantities
    steady_state = p.alpha/p.gamma;
    transient_response = steady_state*(1-exp(-p.gamma*time));

    o.default(j).steady_state = steady_state;
    o.default(j).transient_response = transient_response;
    
    % random
    N = 100;
    r_alpha = random('uniform', 2, 3, [1,N]);
    r_gamma = random('uniform', 2, 3, [1,N]);
    r_steady_state = (p.alpha*r_alpha)./(p.gamma*r_gamma);
    r_transient = 1 - exp(-p.gamma*diag(r_gamma)*repmat(time, N, 1));
    
    o.default(j).r_alpha = r_alpha;
    o.default(j).r_gamma = r_gamma;
    o.default(j).r_steady_state = r_steady_state;
    o.default(j).r_transient = r_transient;

toc
end

save no_feedback_gammap1 o