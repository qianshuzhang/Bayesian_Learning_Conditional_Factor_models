clc;
clear;

%% Data simulation or data loading

sim = 1;       % 1, using simulated data

if sim == 1
   
    T = 5000;
    % annualized mean is 10%
    params.mu     = 0.10/252; 
    % jump size parameters
    params.theta  = -0.02;     params.sigma  = 0.010; 
    % diffusion volatility
    params.omega1 = 1e-10;     params.beta1  = 0.950;   params.alpha1 = 2.2e-3;
    params.delta1 = 0.030;     params.gamma1 = 0.01;
    % jump intensity
    params.omega2 = 1e-8;      params.beta2  = 0.950;   params.alpha2 = 5.0e-4;
    params.delta2 = 2.8e-3;    params.gamma2 = 0.100;
    
    V0 = [0.025/252; 5/252];
    
    [Ret, V, J] = ModelSim(params, V0, T);
    
%     figure(1)
%     subplot(3, 1, 1), plot(Ret)
%     subplot(3, 1, 2), plot(J.y)
%     subplot(3, 1, 3), plot(J.n)
%     
%     figure(2)
%     subplot(2, 1, 1), plot(sqrt(V.h1))
%     subplot(2, 1, 2), plot(V.h2)
%     
else
    
    load('Y');
    T = size(Y, 2);
    
end

%% ------------implement Online Learning-----------------------------------
% -------------------------------------------------------------------------
%particle numbers
Nparam = 100;

% Parameters need to be estimated:
% [delta,theta,alpha0z,alpha1z,beta1z,gammaz,
% alpha0y,alpha1y,beta1y,gammay,lambdaz,lambday];
params = [params.mu; params.theta; params.sigma; params.omega1; params.beta1;...
          params.alpha1; params.delta1; params.gamma1; params.omega2; params.beta2;...
          params.alpha2; params.delta2; params.gamma2]';   

X = repmat(params, Nparam, 1);

States.h1 = 0.025/252 * ones(Nparam, 1);
States.h2 = 5/252 * ones(Nparam, 1);

States.e  = sqrt(States.h1).*randn(Nparam, 1);
States.y  = zeros(Nparam, 1);

%%
h1 = [];
h2 = [];
for t = 1:T
    
    t
    
    [l(:,t), States] = filtering_llh(Ret(t), X, States);
    
    h1 = [h1; mean(States.h1)];
    h2 = [h2; mean(States.h2)];
    
end

figure(1), plot(sum(l')/T)
figure(2), plot(h1), hold on, plot(V.h1, 'r')
figure(3), plot(h2), hold on, plot(V.h2, 'r')