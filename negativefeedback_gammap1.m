function negativefeedback_gammap1

% Calculate temperature dependence using Q10 for negative feedback circuit
% with degradation rate parameter gamma = 0.1/hr


close all 
clear all

% parameters default
o.alpha = 100; %nM/hr
o.gamma = 0.1; %/hr
o.K = 10; %nM

% time
time = 0:1e-1:100; % in hours
o.time = time;
    
% Change parameter sets
M = 100;
set = random('uniform', -1, 1, [3, M]);

o.set = set;

for j=1:M
    j
    tic
    % parameter regime point
    p.alpha = o.alpha*10^set(1,j); %nM/hr
    p.gamma = o.gamma*10^set(2,j); %/hr
    p.K = o.K*10^set(3,j); %nM

    % calculate quantities
    steady_state = (-p.K + sqrt(p.K^2 + 4*p.K*p.alpha/p.gamma))/2;
    transient_response = simulate(p, time);
    
    o.default(j).steady_state = steady_state;
    o.default(j).transient_response = transient_response;

    % random
    N = 100;
    r_alpha = random('uniform', 2, 3, [1,N]);
    r_gamma = random('uniform', 2, 3, [1,N]);
    r_K = random('uniform', 2/3, 3/2, [1,N]);
    r_steady_state = (-p.K*r_K + sqrt((p.K*r_K).^2 + 4*p.K*r_K*p.alpha.*r_alpha./r_gamma/p.gamma))/2;
    
    o.default(j).r_alpha = r_alpha;
    o.default(j).r_gamma = r_gamma;
    o.default(j).r_K = r_K;
    o.default(j).r_steady_state = r_steady_state;
    
    for i = 1:N
        q.alpha = r_alpha(i)*p.alpha;
        q.gamma = r_gamma(i)*p.gamma;
        q.K = r_K(i)*p.K;
        y = simulate(q, time);
        r_transient(i,:) = y/r_steady_state(i);

    end
    o.default(j).r_transient = r_transient;
  
    toc
end

save negative_feedback_gammap1 o

%% Function simulates model
function yy = simulate(p, time)

y0 = 0;

T = time;
% solve the problem using ode23s
[t,y] = ode23s(@f,T,y0,[],p);

tt = t;
yy = y(:,1);

% figure(100);
% plot(tt, yy, 'b');
% hold on;
%-------------------------------------------------
function dydt = f(t,y,p)

X = y(1);


dydt = [p.alpha/(1+(X/p.K)) - p.gamma*X];