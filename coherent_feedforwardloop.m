function coherent_feedforwardloop

% In this code, we estimate temperature dependence of the coherent feedforward loop

close all 
clear all

% a default parameter set
o.alphaX = 100;
o.gammaX = 1; 
o.alphaY = 100; 
o.gammaY = 1;
o.K = 100;

% time
time = 0:1e-2:100; % in hours
o.time = time;
    
% Change parameter sets
M = 100;
set = random('uniform', -1, 1, [5, M]);

o.set = set;

for j=1:M
    j
    tic
    % parameter regime point (for each j there is a parameter set scaled by% 0.1 to 10)
    p.alphaX = o.alphaX*10^set(1,j); %nM/hr
    p.gammaX = o.gammaX*10^set(2,j); %/hr
    p.alphaY = o.alphaY*10^set(3,j); %nM/hr
    p.gammaY = o.gammaY*10^set(4,j); %/hr
    p.K = o.K*10^set(5,j); %nM

    % calculate quantities CFFL
    y_cffl = simulate(p, time); % nominal trajectory
    cffl_ss_0 = (p.alphaY*p.alphaX)/(p.gammaY*p.K*p.gammaX); % initial nominal steady state;
    cffl_ss_f = 4*(p.alphaY*p.alphaX)/(p.gammaY*p.K*p.gammaX); % final nominal steady state;
    cffl_mid = cffl_ss_0+(cffl_ss_f-cffl_ss_0)/2;
    [~, cffl_midx] = min(abs(y_cffl-cffl_mid));
    
    % calculate quantities nominal no feedforward
    y_ref = simulate1(p, time); % nominal trajectory
    ref_ss_0 = (p.alphaY)/(p.gammaY); % initial steady state;
    ref_ss_f = (p.alphaY*2)/(p.gammaY); % final steady state;
    ref_mid = ref_ss_0+(ref_ss_f-ref_ss_0)/2;
    [~,ref_midx] = min(abs(y_ref-ref_mid));

    
%     figure
%     plot(time,y_cffl)
%     hold on
%     plot(time,y_ref)
    
    delay=(cffl_midx-ref_midx)*1e-2; %1e-2 is the grid size
   
    
    o.default(j).cffl_ss_f = cffl_ss_f; % nominal steady state stored to j position
    o.default(j).cffl_response = y_cffl; % nominal transient response stored to j position
    o.default(j).delay = delay; 

    
    
    % random temperature changes about a nominal parameter set
    N = 100; % Number of Q10 scaling
    Q_alphaX = random('uniform', 2, 3, [1,N]); % 100 scaling factor for alphaX
    Q_gammaX = random('uniform', 2, 3, [1,N]); % 100 scaling factor for gammaX
    Q_alphaY = random('uniform', 2, 3, [1,N]); % 100 scaling factor for alphaY
    Q_gammaY = random('uniform', 2, 3, [1,N]); % 100 scaling factor for gammaY
    Q_K = random('uniform', 2/3, 3/2, [1,N]);  % 100 scaling factor for K
    
    o.default(j).Q_alphaX = Q_alphaX;
    o.default(j).Q_gammaX = Q_gammaX;
    o.default(j).Q_alphaY = Q_alphaY;
    o.default(j).Q_gammaY = Q_gammaY;
    o.default(j).Q_K = Q_K;
    
    
    % for CFFL
    Q_cffl_ss=(Q_alphaY.*Q_alphaX)./(Q_K.*Q_gammaX.*Q_gammaY);
    r_cffl_ss_0 = cffl_ss_0*Q_cffl_ss; 
    o.default(j).r_cffl_ss_0 = r_cffl_ss_0;
    r_cffl_ss_f = cffl_ss_f*Q_cffl_ss; 
    o.default(j).r_cffl_ss_f = r_cffl_ss_f; 
    r_cffl_mid = r_cffl_ss_0+(r_cffl_ss_f-r_cffl_ss_0)/2; 
    o.default(j).r_cffl_mid = r_cffl_mid; 
    
    % for Open Loop (no feedforward)
    Q_ref_ss=(Q_alphaY)./(Q_gammaY);
    r_ref_ss_0 = ref_ss_0*Q_ref_ss; 
    o.default(j).r_ref_ss_0 = r_ref_ss_0;
    r_ref_ss_f = ref_ss_f*Q_ref_ss; 
    o.default(j).r_ref_ss_f = r_ref_ss_f; 
    r_ref_mid = r_ref_ss_0+(r_ref_ss_f-r_ref_ss_0)/2; 
    o.default(j).ref_mid = r_ref_mid; 
    
   
    
    for i = 1:N
        q.alphaX = Q_alphaX(i)*p.alphaX;
        q.gammaX = Q_gammaX(i)*p.gammaX;
        q.alphaY = Q_alphaY(i)*p.alphaY;
        q.gammaY = Q_gammaY(i)*p.gammaY;
        q.K = Q_K(i)*p.K;
        
        
        
        r_y_cffl = simulate(q, time);
        [~,r_cffl_midx(i)]=min(abs(r_y_cffl-r_cffl_mid(i)));
        r_cffl_response(:,i)=r_y_cffl;
        
        r_y_ref=simulate1(q,time);
        [~,r_ref_midx(i)]=min(abs(r_y_ref-r_ref_mid(i)));
        r_ref_response(:,i)=r_y_ref;
        
        r_delay(i)=((r_cffl_midx(i)-r_ref_midx(i)))*1e-2;

        
        

    end
    
    
    o.default(j).r_delay=r_delay;
    o.default(j).r_cffl_response=r_cffl_response;
    o.default(j).r_ref_response=r_ref_response;


    
    toc
end

save coherent_ffl o

%% Function simulates model
function yy = simulate(p, time)
u=1;
y0 = (p.alphaY*p.alphaX)/(p.gammaY*p.gammaX*p.K);
x0 = (p.alphaX*u/p.gammaX);
T = time;
% solve the problem using ode23s
[t,y] = ode23s(@f,T,[x0 y0],[],p);

tt = t;
yy = y(:,2);

% figure(100);
% plot(tt, yy, 'b');
% hold on;
%-------------------------------------------------
function dxydt = f(t,xy,p)

x = xy(1);
y = xy(2);
u=2;

dxydt = [p.alphaX*u - p.gammaX*x;...
         p.alphaY*u*x/p.K - p.gammaY*y];
     
%% Function open loop
function yy = simulate1(p, time)
u=1;
x0 = (p.alphaY*u/p.gammaY);
T = time;
% solve the problem using ode23s
[t,x] = ode23s(@f_ol,T,x0,[],p);

tt = t;
yy = x;

% figure(100);
% plot(tt, yy, 'b');
% hold on;
%-------------------------------------------------
function dxydt = f_ol(t,x,p)

u=2;

dxydt = p.alphaY*u - p.gammaY*x;
         