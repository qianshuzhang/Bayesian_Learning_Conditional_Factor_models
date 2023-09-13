function [Ret, V] = ModelSim(params, V0, T,n_factors)


% simulation of non-affine GARCH jump model as specified in Learning
% Comovements of Tails, Volatility and Returns

T0 = 2000;

% number of factors

h = V0;


Ret  = zeros(T, n_factors);
V.h = zeros(T, n_factors);

%xi = exp(params.theta + params.sigma^2/2) - 1;

z = mvnrnd(zeros(1, size(params.A*diag(h)*params.A', 1)), params.A*diag(h)*params.A');
y = 0;
%ztilde = randn;
%ytilde = 0;

Ret(1,:) = params.mu' + z;

V.h(1,:) = h;
for t = 2:T+T0
   
   h = max(params.omega + params.beta.*h + params.alpha.*(Ret(t-1,:)'-params.mu).^2 ...
        , 1e-20);
   

   %ytilde = (y - params.theta*h2)/sqrt((params.theta^2 + params.sigma^2)*h2);
   %ztilde = randn;
   
   z  = mvnrnd(zeros(1, size(params.A*diag(h)*params.A', 1)), params.A*diag(h)*params.A');
   
   Ret(t,:)  = params.mu' + z;
   V.h(t,:) = h';
   
end

Ret  = Ret(T0+1:end,:);
V.h = V.h(T0+1:end,:);

