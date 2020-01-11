function negativefeedback_cellfree_resource_limited

% Calculate temperature dependence using Q10 for negative feedback circuit
% in cell-free context with resource limitation


close all 
clear all

% parameters default
o.alpha = 100; %nM/hr
o.gamma = 1; %/hr
o.K = 10; %nM

% time
time = 0:1e-2:10; % in hours
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
    transient_response = simulate(p, time);
    
    o.default(j).transient_response = transient_response;

    % random
    N = 100;
    r_alpha = random('uniform', 2, 3, [1,N]);
    r_gamma = random('uniform', 2, 3, [1,N]);
    r_K = random('uniform', 2/3, 3/2, [1,N]);
    
    o.default(j).r_alpha = r_alpha;
    o.default(j).r_gamma = r_gamma;
    o.default(j).r_K = r_K;
    
    for i = 1:N
        q.alpha = r_alpha(i)*p.alpha;
        q.gamma = r_gamma(i)*p.gamma;
        q.K = r_K(i)*p.K;
        y = simulate(q, time);
        r_transient(i,:) = y;

    end
    o.default(j).r_transient = r_transient;
  
    toc
end

save negative_feedback_cellfree_resource_limited o

%% Function simulates model
function yy = simulate(p, time)

y0 = [1, 0];

T = time;
% solve the problem using ode23s
[t,y] = ode23s(@f,T,y0,[],p);

tt = t;
yy = y(:,2);

% figure(100);
% plot(tt, yy, 'b');
% hold on;
%-------------------------------------------------
function dydt = f(t,y,p)

B = y(1);
X = y(2);


dydt = [-p.gamma*B
        p.alpha*B/(1+(X/p.K))];