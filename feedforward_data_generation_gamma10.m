function feedforward_data_generation_gamma10

% Calculate temperature dependence using Q10 for an incoerent feedforward
% loop circuit with degradation rate, gamma = 10/hr


close all 
clear all

% parameters default
o.alpha_x = 100; %nM/hr
o.gamma_x = 10; %/hr
o.K_x = 10; %nM
o.alpha_y = 100; %nM/hr
o.gamma_y = 10; %/hr
% u is from 1 to 2


% time
time = 0:1e-3:10; % in hours
o.time = time;
    
% Change parameter sets
M = 100;
set = random('uniform', -1, 1, [5, M]);

o.set = set;

for j=1:M
    j
    tic
    % parameter regime point
    p.alpha_x = o.alpha_x*10^set(1,j); %nM/hr
    p.gamma_x = o.gamma_x*10^set(2,j); %/hr
    p.K_x = o.K_x*10^set(3,j); %nM
    p.alpha_y = o.alpha_y*10^set(4,j); %nM/hr
    p.gamma_y = o.gamma_y*10^set(5,j); %/hr
			
    % calculate quantities
    steady_state = p.K_x*p.alpha_y*p.gamma_x/(p.alpha_x*p.gamma_y);
    transient_response = simulate(p, time);
    
    o.default(j).steady_state = steady_state;
    o.default(j).transient_response = transient_response;

    % random
    N = 100;
    r_alpha_x = random('uniform', 2, 3, [1,N]);
    r_gamma_x = random('uniform', 2, 3, [1,N]);
    r_K_x = random('uniform', 2/3, 3/2, [1,N]);
    r_alpha_y = random('uniform', 2, 3, [1,N]);
    r_gamma_y = random('uniform', 2, 3, [1,N]);
	
    r_steady_state = r_K_x.*p.K_x.*r_alpha_y.*p.alpha_y.*r_gamma_x.*p.gamma_x./(r_alpha_x.*p.alpha_x.*r_gamma_y.*p.gamma_y);
    
    o.default(j).r_alpha_x = r_alpha_x;
    o.default(j).r_gamma_x = r_gamma_x;
    o.default(j).r_K_x = r_K_x;
    o.default(j).r_alpha_y = r_alpha_y;
    o.default(j).r_gamma_y = r_gamma_y;
    o.default(j).r_steady_state = r_steady_state; % so this is not so much a ratio as the steasy-state after temperature change
    
    for i = 1:N
        q.alpha_x = r_alpha_x(i)*p.alpha_x;
        q.gamma_x = r_gamma_x(i)*p.gamma_x;
        q.K_x = r_K_x(i)*p.K_x;
        q.alpha_y = r_alpha_y(i)*p.alpha_y;
        q.gamma_y = r_gamma_y(i)*p.gamma_y;
		
        y = simulate(q, time);
        r_transient(i,:) = y;

        %         figure(3);
%         plot(transient_response/steady_state, r_transient(i,:),'c');
%         hold on;
%         
% 		figure(1);
% 	  	plot(time, r_transient(i,:),'c');
% 		hold on;
    end
    o.default(j).r_transient = r_transient;

    toc
end

save feedforward_data_gamma10 o

%% Function simulates model
function yy = simulate(p, time)

% initial condition is steady-state
y0 = [p.alpha_x/p.gamma_x, p.K_x*p.alpha_y*p.gamma_x/p.alpha_x/p.gamma_y];

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

X = y(1);
Y = y(2);


dydt = [p.alpha_x*2 - p.gamma_x*X
		p.alpha_y*2*p.K_x/X - p.gamma_y*Y];